"""
Steady-state part-load analysis of a single-zone VAV air-cooling system equipped
with a DX air-cooling coil and optionally an air-heating coil.
"""
from collections.abc import Iterator
import warnings
from dataclasses import dataclass
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid, CoolPropWarning
from hvac.cooling_load_calc import ConditionedZone
from hvac.air_conditioning import AirConditioningProcess, AirStream, AdiabaticMixing
from hvac.charts.psychrometric_chart import PsychrometricChart, StatePoint
from .cooling_design import Output
from hvac.logging import ModuleLogger


warnings.filterwarnings('ignore', category=CoolPropWarning)


Q_ = Quantity
logger = ModuleLogger.get_logger(__name__)
logger.setLevel(ModuleLogger.ERROR)


# Define standard (NTP) air:
Air = Fluid('Air')
standard_air = Air(
    T=Q_(20, 'degC'),
    P=Q_(101_325, 'Pa')
)


@dataclass
class DesignData:
    """Groups the design data of the single-zone VAV air-cooling system.

    Attributes
    ----------
    Q_zone_sen: Quantity
        Sensible zone load at design conditions.
    Q_zone_lat: Quantity
        Latent zone load at design conditions.
    V_dot_vent_ntp:
        Minimum required volume flow rate of outdoor ventilation air referred
        to NTP-standard air (dry air at 20 °C and 101.325 kPa).
    V_dot_supply_ntp:
        Required supply air volume flow rate at design conditions referred
        to NTP-standard air.
    outdoor_air:
        State of outdoor air at design conditions.
    zone_air:
        Required state of zone air.
    supply_air:
        Required state of supply air at design conditions.
    ----------------------------------------------------------------------------
    m_dot_vent:
        Outdoor ventilation air mass flow rate. Internally calculated.
    m_dot_supply:
        Supply air mass flow rate. Internally calculated.
    m_dot_recir:
        Recirculation air mass flow rate. Internally calculated.
    V_dot_recir_ntp:
        Recirculation air volume flow rate. Internally calculated.
    """
    Q_zone_sen: Quantity
    Q_zone_lat: Quantity
    V_dot_vent_ntp: Quantity
    V_dot_supply_ntp: Quantity
    outdoor_air: HumidAir
    zone_air: HumidAir
    supply_air: HumidAir

    def __post_init__(self):
        self.m_dot_vent = self.V_dot_vent_ntp * standard_air.rho
        self.m_dot_supply = self.V_dot_supply_ntp * standard_air.rho
        self.m_dot_recir = self.m_dot_supply - self.m_dot_vent
        self.V_dot_recir_ntp = self.m_dot_recir / standard_air.rho


class VAVSingleZoneAirCoolingSystem:
    """
    Model for analyzing the steady-state part-load operation of a single-zone
    VAV air-cooling system equipped with an air-cooling coil and optionally an
    air-heating coil.

    The aim is to maintain the zone air setpoint temperature under the given
    set of external conditions (outdoor air state, sensible and latent zone
    load). To accomplish this, the mass flow rate of supply air is modulated,
    while the cooled air leaving the cooling coil remains at a fixed state
    (it is assumed that both temperature and humidity remain fixed).
    However, the supply air mass flow rate cannot drop below a certain minimum
    limit (to ensure proper mixing of supply air with zone air), and it also
    cannot become greater than the maximum capacity of the supply fan.
    If the minimum limit of the supply air mass flow is reached, the cooled air
    needs to be heated if a heating coil is present or the setpoint temperature
    of the cooled air needs to be increased. Otherwise, the zone air temperature
    will decrease below its setpoint.
    If the maximum limit of the supply air mass flow rate is reached, the zone
    air temperature will increase beyond its setpoint, unless the setpoint
    temperature of the cooled air can be decreased. Otherwise, the zone air
    temperature will increase above its setpoint.

    Mixing control (economizer operation) may allow reaching the cooling air
    temperature setpoint without mechanical cooling if outdoor air temperature
    is lower than this cooling air temperature setpoint.
    If outdoor air temperature is between the cooling air temperature setpoint
    and the zone air setpoint temperature, the temperature of the mixed air will
    be higher than the outdoor air temperature, and it will cost more cooling
    energy to cool this air to the cooling air temperature setpoint. Therefore,
    recirculation of return air is set to a minimum (zero), and all supply air
    to the zone will come from outdoors.
    If outdoor air temperature is above the zone air temperature setpoint, the
    outdoor air flow rate is restricted to what's required for zone ventilation
    only and recirculation of return air is now at its maximum.
    """
    def __init__(
        self,
        design_data: DesignData,
        T_set_cool: Quantity | None = None,
        W_set_cool: Quantity | None = None,
        heating_coil_present: bool = True,
        variable_cooling_setpoint: bool = False,
        rel_m_dot_supply_min: Quantity = Q_(60, 'pct'),
        units: dict[str, tuple[str, int]] | None = None
    ) -> None:
        """
        Creates a `VAVSingleZoneAirCoolingSystem` object.

        Parameters
        ----------
        design_data:
            `DesignData` instance that groups the design data of the AC system
            (see class `DesignData` in this module).
        T_set_cool: optional
            Setpoint temperature of the air leaving the cooling coil. If None,
            the supply air temperature in `design_data` will be used.
        W_set_cool: optional
            The humidity ratio of the air leaving the cooling coil that belongs
            with the setpoint temperature of the cooled air. If None, the supply
            air humidity ratio in `design_data` will be used.
        heating_coil_present: default True
            Indicates if a heating coil is present downstream of the cooling
            coil or not.
        variable_cooling_setpoint: default False
            Indicates if the setpoint of the cooling air at the cooling coil
            outlet can be shifted in case no heating coil is present.
        rel_m_dot_supply_min: default 60 %
            The minimum mass flow rate of supply air to ensure proper mixing
            with zone air expressed as a fraction of the supply air mass flow
            rate determined for summer peak design conditions.
        units: optional
            Units to be used for displaying the results. Keys are:
            'm_dot' for mass flow rate, 'V_dot' for volume flow rate, 'Q_dot'
            for heat transfer rate, 'T' for temperature, 'W' for absolute
            humidity, 'RH' for relative humidity, and 'SHR' for the sensible heat
            ratio. Values are 2-tuples: the first element is the unit
            (a string), the second element is the number of decimals to be
            displayed (an int).
        """
        self.design_data = design_data
        # Set the normal cooling coil setpoint:
        self.T_set_cool = (
            T_set_cool if T_set_cool is not None
            else self.design_data.supply_air.Tdb
        )
        self.W_set_cool = (
            W_set_cool if W_set_cool is not None
            else self.design_data.supply_air.W
        )
        # Indicate the control configuration of the air-cooling system:
        self.heating_coil_present = heating_coil_present
        if not self.heating_coil_present and variable_cooling_setpoint:
            self.variable_cooling_setpoint = True
        else:
            # If a heating coil is present, the cooling coil setpoint is always
            # fixed, even if the user should have indicated otherwise.
            self.variable_cooling_setpoint = False
        # Set the allowable minimum flow rate of the supply air:
        self.rel_m_dot_supply_min = rel_m_dot_supply_min
        # The minimum supply air mass flow rate must not be less than the
        # minimum ventilation air requirement and also must not be less than the
        # indicated fraction of the design value to ensure proper mixing of
        # supply air with zone air:
        self.m_dot_supply_min = max(
            self.design_data.m_dot_vent,
            rel_m_dot_supply_min * self.design_data.m_dot_supply
        )
        # The maximum supply air mass flow rate is limited to its design value,
        # which is considered to be the maximum flow rate that the supply fan can
        # displace:
        self.m_dot_supply_max = self.design_data.m_dot_supply
        # Determine the ADP of the cooling coil at design conditions.
        # It will be assumed that the ADP of the cooling coil is a constant.
        # (--> To what extent is this assumption correct?)
        mixing_chamber_des = AdiabaticMixing(
            in1=AirStream(
                state=self.design_data.outdoor_air,
                m_da=self.design_data.m_dot_vent
            ),
            in2=AirStream(
                state=self.design_data.zone_air,
                m_da=self.design_data.m_dot_recir
            ),
            out=AirStream(m_da=self.design_data.m_dot_supply)
        )
        mixed_air_des = mixing_chamber_des.stream_out.state
        cooling_coil_des = AirConditioningProcess(
            air_in=mixed_air_des,
            air_out=self.design_data.supply_air,
            m_da=self.design_data.m_dot_supply,
            Q_sen=self.design_data.Q_zone_sen,
            Q_lat=self.design_data.Q_zone_lat
        )
        self.ADP_cc = cooling_coil_des.ADP
        # Set the wanted measuring units of the output values:
        self.units = units or {}
        # Define internal instance variables:
        self.outdoor_air: HumidAir | None = None
        self.rel_Q_zone_sen: Quantity | None = None
        self.rel_Q_zone_lat: Quantity | None = None
        self.Q_dot_zone_sen: Quantity | None = None
        self.Q_dot_zone_lat: Quantity | None = None
        self.SHR_zone: Quantity | None = None
        self.T_zone_sp: Quantity | None = None
        self.m_dot_supply: Quantity | None = None
        self.m_dot_recir: Quantity | None = None
        self.m_dot_vent: Quantity | None = None
        self.cooled_air: HumidAir | None = None
        self.supply_air: HumidAir | None = None
        self.return_air: HumidAir | None = None
        self.mixed_air: HumidAir | None = None
        self.Q_dot_cc: Quantity | None = None
        self.SHR_cc: Quantity | None = None
        self.Q_dot_hc: Quantity | None = None

    def _control_supply_fan(
        self,
        supply_air: HumidAir,
    ) -> tuple[Quantity, HumidAir, HumidAir]:
        """Determines with the given state of supply air, the supply air mass
        flow rate to maintain the zone air setpoint temperature.

        However, the supply air mass flow rate cannot become smaller than the
        minimum limit, neither become greater than the maximum limit. If one
        of both limiting cases is met, the required state of supply air will
        also need to be recalculated.

        Returns
        -------
        m_dot_supply:
            The supply air mass flow rate needed to maintain the zone air
            setpoint temperature, or the minimum or maximum limit of the supply
            air mass flow rate.
        supply_air:
            The final state of supply air needed to maintain the zone air
            temperature (if possible, depending on the resulting supply air mass
            flow rate and the configuration of the air-cooling control system.)
        return_air:
            The state of the return air that depends on the sensible and
            latent cooling load of the zone, the supply air mass flow rate, and
            the final state of the supply air.
        """
        # Calculate the mass flow rate of supply air needed to maintain the
        # zone air setpoint temperature with the cooling coil being operated at
        # its normal setpoint.
        # The moisture content of the zone and return air will depend on the
        # moisture content of the supply air and on the latent load of the zone.
        zone = AirConditioningProcess(
            air_in=supply_air,
            T_ao=self.T_zone_sp,
            Q_sen=self.Q_dot_zone_sen,
            Q_lat=self.Q_dot_zone_lat,
            h_w=Q_(0, 'J / kg')
        )
        m_dot_supply = zone.m_da
        return_air = zone.air_out
        if m_dot_supply < self.m_dot_supply_min:
            # The mass flow rate of supply air cannot be smaller than the minimum
            # limit.
            m_dot_supply = self.m_dot_supply_min
            # We also need to recalculate the required state of supply air to
            # maintain the zone air temperature and the resulting state of the
            # return air from the zone (its humidity) depending on the
            # configuration of the air-cooling control system.
            if self.heating_coil_present:
                # To maintain the zone air setpoint temperature, the supply air
                # needs to be heated (see `_control_heating_coil`).
                # Heating the supply-air does not change its moisture content,
                # only its temperature is raised.
                zone = AirConditioningProcess(
                    W_ai=supply_air.W,
                    T_ao=self.T_zone_sp,
                    m_da=m_dot_supply,
                    Q_sen=self.Q_dot_zone_sen,
                    Q_lat=self.Q_dot_zone_lat,
                    h_w=Q_(0, 'J / kg')
                )
                supply_air = zone.air_in  # new required state of supply air
                return_air = zone.air_out
            elif self.variable_cooling_setpoint:
                # If no heating coil is present, but the cooling coil controller
                # is able to increase the cooling coil setpoint temperature to
                # maintain the zone air setpoint temperature.
                # (--> To what extent is this always possible?)
                zone = AirConditioningProcess(
                    W_ai=supply_air.W,
                    # Actually, `W_ai` will also change because the cooling
                    # process line is sloped --> see `_control_cooling_coil`.
                    T_ao=self.T_zone_sp,
                    m_da=m_dot_supply,
                    Q_sen=self.Q_dot_zone_sen,
                    Q_lat=self.Q_dot_zone_lat,
                    h_w=Q_(0, 'J / kg')
                )
                supply_air = zone.air_in  # new required state of supply air
                return_air = zone.air_out
            else:
                # The zone air setpoint temperature cannot be maintained.
                # Use outdoor air as supply air and turn off the cooling
                # coil --> see `_control_cooling_coil`.
                supply_air = self.outdoor_air  # new state of supply air
                zone = AirConditioningProcess(
                    air_in=supply_air,
                    m_da=m_dot_supply,
                    Q_sen=self.Q_dot_zone_sen,
                    Q_lat=self.Q_dot_zone_lat,
                    h_w=Q_(0, 'J / kg')
                )
                return_air = zone.air_out
                try:
                    _ = return_air.RH
                except ValueError:
                    return_air = HumidAir(
                        Tdb=return_air.Tdb,
                        RH=Q_(100, 'pct')
                    )
        # The mass flow rate of supply air cannot be greater than the minimum
        # limit:
        if m_dot_supply > self.m_dot_supply_max:
            m_dot_supply = self.m_dot_supply_max
            if self.variable_cooling_setpoint:
                # The cooling coil controller is able to decrease the cooling
                # coil setpoint temperature to maintain the zone air setpoint
                # temperature.(--> To what extent is this possible?)
                zone = AirConditioningProcess(
                    W_ai=supply_air.W,
                    # Actually, `W_ai` will also change because the cooling
                    # process line is sloped --> see `_control_cooling_coil`.
                    T_ao=self.T_zone_sp,
                    m_da=m_dot_supply,
                    Q_sen=self.Q_dot_zone_sen,
                    Q_lat=self.Q_dot_zone_lat,
                    h_w=Q_(0, 'J / kg')
                )
                supply_air = zone.air_in  # new required state of supply air
                return_air = zone.air_out
            else:
                # The zone air temperature will inevitably raise above its
                # setpoint.
                zone = AirConditioningProcess(
                    air_in=supply_air,
                    m_da=m_dot_supply,
                    Q_sen=max(Q_(1e-12, 'W'), self.Q_dot_zone_sen),
                    Q_lat=self.Q_dot_zone_lat,
                    h_w=Q_(0, 'J / kg')
                )
                return_air = zone.air_out
        return m_dot_supply, return_air, supply_air

    def _control_mixing_air(
        self,
        m_dot_supply: Quantity,
        return_air: HumidAir
    ) -> tuple[HumidAir, Quantity, Quantity]:
        """Determines the state of the mixed air, the mass flow rate of outdoor
        ventilation air and the mass flow rate of recirculated return air.

        There are three possible cases to consider:
        1. The outdoor air temperature is equal or higher than the zone air
           setpoint temperature. The mass flow rate of outdoor air is
           restricted to the minimum ventilation requirement. The mass flow rate
           of recirculation air on the other hand is 100 %.
        2. The outdoor air temperature is equal or lower than the cooling coil
           setpoint temperature. In that case, the cooling coil setpoint
           temperature can be attained by mixing outdoor air with return air and
           no mechanical cooling is needed.
        3. The outdoor air temperature is between the cooling coil setpoint
           temperature and the zone air setpoint temperature. The lowest
           possible temperature leaving the mixing chamber, i.e., closest to the
           cooling coil setpoint temperature, will be the temperature of the
           outdoor air. So, no mixing takes place: the mass flow rate of supply
           air comes directly from outdoors. The mass flow rate of recirculation
           air is 0 %.

        Returns
        -------
        mixed_air:
            The state of the air leaving the mixing chamber and entering the
            cooling coil.
        m_dot_vent:
            The mass flow rate of supply air taken in from outdoors.
        m_dot_recir:
            The mass flow rate of supply air recirculated from the zone.
        """
        if self.outdoor_air.Tdb >= self.T_zone_sp:
            # If the outdoor air temperature is equal or higher than the zone air
            # setpoint temperature: outdoor ventilation air flow rate is
            # restricted to the minimum ventilation rate required (determined
            # when designing the system) --> OA damper in its minimum open
            # position (NC) and RA damper in its maximum open position (NO).
            m_dot_vent = self.design_data.m_dot_vent
            m_dot_recir = m_dot_supply - m_dot_vent
            mixing_chamber = AdiabaticMixing(
                in1=AirStream(
                    state=return_air,
                    m_da=m_dot_recir
                ),
                in2=AirStream(
                    state=self.outdoor_air,
                    m_da=self.design_data.m_dot_vent
                ),
                out=AirStream(m_da=m_dot_supply)
            )
            mixed_air = mixing_chamber.stream_out.state
        elif self.outdoor_air.Tdb <= self.T_set_cool:
            # If outdoor air temperature is equal or less than the cooling coil
            # setpoint temperature, no mechanical cooling needed. Outdoor air is
            # mixed with return air to get at the cooling coil setpoint
            # temperature --> the mixing controller controls the position of the
            # OA and RA damper.
            mixing_chamber = AdiabaticMixing(
                in1=AirStream(
                    state=return_air
                ),
                in2=AirStream(
                    state=self.outdoor_air
                ),
                out=AirStream(Tdb=self.T_set_cool, m_da=m_dot_supply)
            )
            m_dot_vent = mixing_chamber.stream_in2.m_da
            m_dot_recir = mixing_chamber.stream_in1.m_da
            mixed_air = mixing_chamber.stream_out.state
        else:
            # If the outdoor air temperature is between the cooling coil setpoint
            # temperature and the zone air setpoint temperature, mixing cannot
            # attain the cooling coil setpoint temperature, because mixing
            # outdoor air with return air would inevitably raise the mixed air
            # temperature more above the cooling coil setpoint temperature.
            # In fact, it will cost less energy to cool the outdoor air to the
            # cooling coil temperature. So, no mixing takes place --> OA damper
            # in its maximum open and RA damper in its minimum open position.
            m_dot_vent = m_dot_supply
            m_dot_recir = Q_(0.0, 'kg / s')
            mixed_air = self.outdoor_air
        return mixed_air, m_dot_vent, m_dot_recir

    def _control_cooling_coil(
        self,
        m_dot_supply: Quantity,
        mixed_air: HumidAir,
        supply_air_req: HumidAir
    ) -> tuple[HumidAir, AirConditioningProcess | None]:
        """Determines if the cooling coil can be operational, and the state of
        air leaving the cooling coil.

        Parameters
        ----------
        m_dot_supply:
            The mass flow rate of supply air to the zone that flows through the
            cooling coil. This was determined in method `_control_supply fan()`.
        mixed_air:
            The state of air leaving the mixing chamber and entering the cooling
            coil. This was determined in method `_control_mixing_air()`.
        supply_air_req:
            The state of supply air required to maintain the zone air setpoint
            temperature, considered to be also the state of air leaving the
            cooling coil. This was determined in method `_control_supply fan()`.

        Returns
        -------
        cooled_air:
            The state of air leaving the cooling coil and entering the heating
            coil if it is present in the air-cooling system.
        cooling_coil:
            If mechanical cooling is possible, the cooling coil, being an
            `AirConditioningProcess` object is returned. Otherwise, None is
            returned.
        """
        if (m_dot_supply == self.m_dot_supply_min) and not self.heating_coil_present:
            # If the mass flow rate of supply air is at its minimum limit and
            # no heating coil is present, the zone air temperature will drop
            # below its setpoint if the cooling coil is active, unless the
            # cooling coil controller can increase the cooling coil setpoint
            # temperature to the required supply air temperature that's needed
            # to maintain the zone air setpoint temperature.
            if self.variable_cooling_setpoint:
                if mixed_air.Tdb > supply_air_req.Tdb:
                    # Only if the mixed air has a higher temperature than the
                    # required supply air temperature, the mixed air can be
                    # cooled to the required supply air temperature.
                    cooling_coil = AirConditioningProcess(
                        air_in=mixed_air,
                        m_da=m_dot_supply,
                        T_ao=supply_air_req.Tdb,
                        ADP=self.ADP_cc
                    )
                    cooled_air = cooling_coil.air_out
                else:
                    cooling_coil = None
                    cooled_air = mixed_air
            else:
                # The cooling coil is turned off when no heating coil is present
                # to prevent the zone air temperature of becoming too low,
                # whatever the temperature of the mixed air might be.
                cooling_coil = None
                cooled_air = mixed_air
        elif (m_dot_supply == self.m_dot_supply_max) and self.variable_cooling_setpoint:
            # If the mass flow rate of supply air is at its maximum limit, the
            # zone air temperature will raise above its setpoint, unless the
            # cooling coil controller can decrease the cooling coil setpoint
            # temperature to the required supply air temperature that's needed
            # to maintain the zone air setpoint temperature.
            if mixed_air.Tdb > supply_air_req.Tdb:
                cooling_coil = AirConditioningProcess(
                    air_in=mixed_air,
                    m_da=m_dot_supply,
                    T_ao=supply_air_req.Tdb,
                    ADP=self.ADP_cc
                )
                cooled_air = cooling_coil.air_out
            else:
                cooling_coil = None
                cooled_air = mixed_air
        elif mixed_air.Tdb > self.T_set_cool:
            # If the mass flow rate of supply air is between its minimum and
            # maximum limit, and the mixed air temperature is higher than the
            # cooling air setpoint temperature, the cooling coil controller
            # keeps the temperature of the cooled air fixed at its normal setpoint
            # temperature `T_set_cool`.
            cooling_coil = AirConditioningProcess(
                air_in=mixed_air,
                m_da=m_dot_supply,
                T_ao=self.T_set_cool,
                ADP=self.ADP_cc
            )
            cooled_air = cooling_coil.air_out
        else:
            # The cooling coil is turned off.
            cooling_coil = None
            cooled_air = mixed_air
        return cooled_air, cooling_coil

    def _control_heating_coil(
        self,
        m_dot_supply: Quantity,
        cooled_air: HumidAir,
        supply_air_req: HumidAir,
    ) -> tuple[HumidAir, AirConditioningProcess | None]:
        """If a heating coil is present, determines if it needs to be
        operational. Also returns the state of air leaving the heating coil, or,
        if no heating coil is present, that is leaving the cooling coil.

        Parameters
        ----------
        m_dot_supply:
            The mass flow rate of supply air to the zone that flows through the
            cooling coil. This was determined in method `_control_supply fan()`.
        cooled_air:
            The state of air leaving the cooling coil and entering the heating
            coil.
        supply_air_req:
            The state of supply air required to maintain the zone air setpoint
            temperature, considered to be also the state of air leaving the
            cooling coil. This was determined in method `_control_supply fan()`.

        Returns
        -------
        supply_air:
            The state of air leaving the heating coil if present. Otherwise,
            the state of air leaving the cooling coil is returned. This is also
            the state of air entering the zone.
        heating_coil:
            If a heating coil is present and the heating coil needs to be active
            to get at the required supply air temperature to maintain the zone
            air setpoint temperature, the heating coil, being an
            `AirConditioningProcess` object is returned. Otherwise, None is
            returned.
        """
        if self.heating_coil_present and (cooled_air.Tdb < supply_air_req.Tdb):
            # If the cooling coil outlet air temperature is less than the
            # required supply air temperature to maintain the setpoint zone air
            # temperature, supply air is heated to keep the zone air temperature
            # at its setpoint.
            supply_air = HumidAir(Tdb=supply_air_req.Tdb, W=cooled_air.W)
            heating_coil = AirConditioningProcess(
                air_in=cooled_air,
                air_out=supply_air,
                m_da=m_dot_supply
            )
        else:
            # If no heating coil is present or the cooled air leaving the cooling
            # coil doesn't need heating, return the state of air entering the
            # heating coil unaltered.
            supply_air = cooled_air
            heating_coil = None
        return supply_air, heating_coil

    def _zone_air(
        self,
        m_dot_supply: Quantity,
        supply_air: HumidAir
    ) -> HumidAir:
        """Determines the final state of the zone and return air in the space
        based on (1) the final state of supply air leaving the cooling or
        heating coil and that is entering the zone, (2) the sensible and latent
        cooling load of the zone, and (3) the mass flow rate of supply air.

        Parameters
        ----------
        m_dot_supply:
            The mass flow rate of supply air to the zone that flows through the
            cooling coil. This was determined in method `_control_supply fan()`.
        supply_air:
            The final state of air entering the zone.
        """
        zone = AirConditioningProcess(
            air_in=supply_air,
            m_da=m_dot_supply,
            Q_sen=self.Q_dot_zone_sen,
            Q_lat=self.Q_dot_zone_lat,
            h_w=Q_(0, 'J / kg')
        )
        return_air = zone.air_out
        try:
            _ = return_air.RH
        except ValueError:
            return_air = HumidAir(Tdb=return_air.Tdb, RH=Q_(100, 'pct'))
        return return_air

    def analyze(
        self,
        outdoor_air: HumidAir,
        rel_Q_zone_sen: Quantity,
        rel_Q_zone_lat: Quantity,
        T_zone_sp: Quantity | None = None
    ) -> Output:
        """Analyzes the operating state of the system for the external
        conditions given.

        Parameters
        ----------
        outdoor_air:
            The state of outdoor air at which the cooling load of the zone
            was calculated.
        rel_Q_zone_sen:
            Fraction of the sensible zone load at design conditions (either
            a fraction between 0 and 1, or a percentage between 0 and 100 %).
        rel_Q_zone_lat:
            Fraction of the latent zone load at design conditions (either
            a fraction between 0 and 1, or a percentage between 0 and 100 %).
        T_zone_sp: optional
            Setpoint of the zone air temperature at which the cooling load of
            the zone was calculated. If None, the zone air temperature in
            `design_data` will be used.

        Returns
        -------
        `Output` object (see class `Output` in module `cooling_design`).
        """
        self.outdoor_air = outdoor_air
        self.rel_Q_zone_sen = rel_Q_zone_sen
        self.rel_Q_zone_lat = rel_Q_zone_lat
        self.Q_dot_zone_sen = self.rel_Q_zone_sen * self.design_data.Q_zone_sen
        self.Q_dot_zone_lat = self.rel_Q_zone_lat * self.design_data.Q_zone_lat
        self.SHR_zone = self.Q_dot_zone_sen / (self.Q_dot_zone_sen + self.Q_dot_zone_lat)
        self.T_zone_sp = T_zone_sp if T_zone_sp is not None else self.design_data.zone_air.Tdb

        # We initially assume that the cooling coil is active and that the
        # cooling coil controller keeps the temperature of the air leaving
        # the cooling coil and being supplied to the zone at its normal setpoint.
        supply_air_ini = HumidAir(Tdb=self.T_set_cool, W=self.W_set_cool)
        logger.debug(
            "supply air (initial): "
            f"{supply_air_ini.Tdb.to('degC'):~P.2f} DB, "
            f"{supply_air_ini.W.to('g / kg'):~P.2f} AH "
            f"({supply_air_ini.RH.to('pct'):~P.0f} RH)"
        )
        # With the initial supply air state being selected, we determine (1) the
        # needed mass flow rate of supply air to maintain the zone air
        # temperature at its setpoint, (2) the resulting state of return air,
        # and (3) the required new state of supply air if the needed mass flow
        # rate of supply air cannot be attained with the initially assumed state
        # of supply air (i.e., if the supply air mass flow rate would be below
        # its minimum or above its maximum limit).
        m_dot_supply, return_air, supply_air_req = self._control_supply_fan(supply_air_ini)
        logger.debug(
            "supply air (required): "
            f"{supply_air_req.Tdb.to('degC'):~P.2f} DB, "
            f"{supply_air_req.W.to('g / kg'):~P.2f} AH "
            f"({supply_air_req.RH.to('pct'):~P.0f} RH)"
        )
        logger.debug(
            "return air (initial): "
            f"{return_air.Tdb.to('degC'):~P.2f} DB, "
            f"{return_air.W.to('g / kg'):~P.2f} AH "
            f"({return_air.RH.to('pct'):~P.0f} RH)"
        )
        # Mixing control:
        mixed_air, m_dot_vent, m_dot_recir = self._control_mixing_air(
            m_dot_supply,
            return_air
        )
        per_m_dot_recir = m_dot_recir.to('kg / hr') / m_dot_supply.to('kg / hr')
        per_m_dot_vent = m_dot_vent.to('kg / hr') / m_dot_supply.to('kg / hr')
        logger.debug(
            "percent mass flow rate recirculated air: "
            f"{per_m_dot_recir.to('pct'):~P.0f}"
        )
        logger.debug(
            "percent mass flow rate outdoor air: "
            f"{per_m_dot_vent.to('pct'):~P.0f}"
        )
        logger.debug(
            "mixed air: "
            f"{mixed_air.Tdb.to('degC'):~P.2f} DB, "
            f"{mixed_air.W.to('g / kg'):~P.2f} AH "
            f"({mixed_air.RH.to('pct'):~P.0f} RH)"
        )
        # Cooling coil control:
        cooled_air, cooling_coil = self._control_cooling_coil(
            m_dot_supply,
            mixed_air,
            supply_air_req
        )
        logger.debug(
            "cooled air: "
            f"{cooled_air.Tdb.to('degC'):~P.2f} DB, "
            f"{cooled_air.W.to('g / kg'):~P.2f} AH "
            f"({cooled_air.RH.to('pct'):~P.0f} RH)"
        )
        # Heating coil control:
        supply_air, heating_coil = self._control_heating_coil(
            m_dot_supply,
            cooled_air,
            supply_air_req
        )
        logger.debug(
            "supply air (final): "
            f"{supply_air.Tdb.to('degC'):~P.2f} DB, "
            f"{supply_air.W.to('g / kg'):~P.2f} AH "
            f"({supply_air.RH.to('pct'):~P.0f} RH)"
        )
        # Zone:
        return_air = self._zone_air(
            m_dot_supply,
            supply_air
        )
        logger.debug(
            "return air (final): "
            f"{return_air.Tdb.to('degC'):~P.2f} DB, "
            f"{return_air.W.to('g / kg'):~P.2f} AH "
            f"({return_air.RH.to('pct'):~P.0f} RH)"
        )
        # Assign the results to the instance variables of this object:
        self.m_dot_supply = m_dot_supply
        self.m_dot_vent = m_dot_vent
        self.m_dot_recir = m_dot_recir
        self.cooled_air = cooled_air
        self.supply_air = supply_air
        self.return_air = return_air
        self.mixed_air = mixed_air
        if cooling_coil is None:
            self.Q_dot_cc = Q_(0.0, 'W')
            self.SHR_cc = Q_(0.0, 'frac')
        else:
            self.Q_dot_cc = -cooling_coil.Q
            self.SHR_cc = cooling_coil.SHR
        if heating_coil is None:
            self.Q_dot_hc = Q_(0.0, 'W')
        else:
            self.Q_dot_hc = heating_coil.Q_sen
        # If both `cooling_coil` and `heating_coil` are None, the air-cooling
        # system is not operating. In that case, the actual state of the zone
        # air and the return air cannot be calculated from steady-state
        # equations, as the single-zone building will evolve towards a thermal
        # equilibrium with its outdoor environment.
        if cooling_coil is None and heating_coil is None:
            self.return_air = HumidAir(
                Tdb=self.outdoor_air.Tdb,
                W=self.return_air.W
            )
        return Output(
            self.m_dot_supply, self.m_dot_vent, self.m_dot_recir,
            self.outdoor_air, self.mixed_air, self.cooled_air,
            self.supply_air, self.return_air, self.design_data.zone_air,
            self.Q_dot_zone_sen, self.Q_dot_zone_lat, self.SHR_zone,
            self.Q_dot_cc, self.SHR_cc, self.Q_dot_hc, self.units
        )

    @property
    def psychrometric_chart(self) -> PsychrometricChart:
        """Returns a psychrometric chart with the AC processes drawn on it."""
        psy_chart = PsychrometricChart()
        psy_chart.plot_process(
            name='air mixing',
            start_point=StatePoint(
                self.outdoor_air.Tdb,
                self.outdoor_air.W
            ),
            end_point=StatePoint(
                self.return_air.Tdb,
                self.return_air.W
            ),
            mix_point=StatePoint(
                self.mixed_air.Tdb,
                self.mixed_air.W
            )
        )
        psy_chart.plot_process(
            name='air cooling',
            start_point=StatePoint(
                self.mixed_air.Tdb,
                self.mixed_air.W
            ),
            end_point=StatePoint(
                self.cooled_air.Tdb,
                self.cooled_air.W
            )
        )
        psy_chart.plot_process(
            name='air heating',
            start_point=StatePoint(
                self.cooled_air.Tdb,
                self.cooled_air.W
            ),
            end_point=StatePoint(
                self.supply_air.Tdb,
                self.supply_air.W
            )
        )
        psy_chart.plot_process(
            name='zone heating and humidification',
            start_point=StatePoint(
                self.supply_air.Tdb,
                self.supply_air.W
            ),
            end_point=StatePoint(
                self.return_air.Tdb,
                self.return_air.W
            )
        )
        return psy_chart


class CoolingSimData:
    """Helper class to prepare the simulation data needed to analyze the
    `VAVSingleZoneAirCoolingSystem` model on an hourly basis on a given day.
    """
    def __init__(
        self,
        zone: ConditionedZone,
        design_data: DesignData
    ) -> None:
        """Creates an instance of class `CoolingSimData`.

        Parameters
        ----------
        zone:
            `ConditionedZone` object that represents the thermal model of the
            single-zone building.
        design_data:
            `DesignData` object that contains data from the design
            calculations of the single-zone VAV air-cooling system, including
            the sensible and latent design cooling load.
        """
        # Get the outdoor air dry-bulb and wet-bulb temperatures for each hour
        # of the day (index 0 to 23):
        self.T_outdoor_db_rng = zone.weather_data.T_db_prof
        self.T_outdoor_wb_rng = zone.weather_data.T_wb_prof
        # Get the zone air setpoint temperature for each hour of the day
        # (index 0 to 23):
        self.T_set_zone_rng = [
            zone.T_zone(t * 3600)
            for t in range(len(self.T_outdoor_db_rng))
        ]
        # Get the dataframe with the zone loads (heat gains) for each hour of
        # the day:
        self.df_Q_zone = zone.solve(unit='kW')
        self.design_data = design_data
        # Combine outdoor air dry-bulb and wet-bulb temperatures in `HumidAir`
        # objects:
        self.outdoor_air_rng = [
            HumidAir(Tdb=T_outdoor_db, Twb=T_outdoor_wb)
            for T_outdoor_db, T_outdoor_wb
            in zip(self.T_outdoor_db_rng, self.T_outdoor_wb_rng)
        ]
        # 'Quantify' and put the sensible zone loads from the dataframe in a
        # list (index 0 to 23):
        self.Q_zone_sen_rng = [
            Q_(Q, 'kW')
            for Q in self.df_Q_zone['Q_dot_sen_zone'].values
        ]
        # 'Quantify' and put the latent zone loads from the dataframe in a
        # list (index 0 to 23):
        self.Q_zone_lat_rng = [
            Q_(Q, 'kW')
            for Q in self.df_Q_zone['Q_dot_lat_zone'].values
        ]
        # Determine the relative sensible zone loads (i.e., as a fraction of the
        # design load):
        self.rel_Q_zone_sen_rng = [
            Q_zone_sen.to('kW') / self.design_data.Q_zone_sen.to('kW')
            for Q_zone_sen in self.Q_zone_sen_rng
        ]
        # Determine the relative latent zone loads (i.e., as a fraction of the
        # design load):
        self.rel_Q_zone_lat_rng = [
            Q_zone_lat.to('kW') / self.design_data.Q_zone_lat.to('kW')
            for Q_zone_lat in self.Q_zone_lat_rng
        ]

    def __call__(self) -> Iterator[list[tuple[HumidAir, Quantity, Quantity]]]:
        """Returns an iterator over tuples which have three elements:
        -   the state of outdoor air
        -   the zone air setpoint temperature
        -   the corresponding relative sensible cooling load (ratio of the
            sensible part-load to the sensible design-load)
        -   the corresponding relative latent cooling load (ratio of the latent
            part-load to the latent design-load)

        These tuples can be used as arguments in the method `analyze` of class
        `VAVSingleZoneAirCoolingSystem`.
        """
        return zip(
            self.outdoor_air_rng,
            self.T_set_zone_rng,
            self.rel_Q_zone_sen_rng,
            self.rel_Q_zone_lat_rng
        )

    def get_sim_data(
        self,
        hour: int,
        is_relative: bool = True
    ) -> tuple[HumidAir, Quantity, Quantity, Quantity]:
        """Returns for the given hour (index 0 to 23) a tuple of three elements:
        -   the state of outdoor air
        -   the zone air setpoint temperature
        -   the corresponding relative sensible cooling load (ratio of the
            sensible part-load to the sensible design-load)
        -   the corresponding relative latent cooling load (ratio of the latent
            part-load to the latent design-load)
        """
        outdoor_air = self.outdoor_air_rng[hour]
        T_set_zone = self.T_set_zone_rng[hour]
        if is_relative:
            rel_Q_zone_sen = self.rel_Q_zone_sen_rng[hour]
            rel_Q_zone_lat = self.rel_Q_zone_lat_rng[hour]
            return outdoor_air, T_set_zone, rel_Q_zone_sen, rel_Q_zone_lat
        else:
            Q_zone_sen = self.Q_zone_sen_rng[hour]
            Q_zone_lat = self.Q_zone_lat_rng[hour]
            return outdoor_air, T_set_zone, Q_zone_sen, Q_zone_lat
