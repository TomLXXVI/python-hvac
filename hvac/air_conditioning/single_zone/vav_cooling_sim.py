"""
Steady-state part-load analysis of a single-zone VAV air-cooling system that can
both cool and heat the supply air to keep the zone air temperature fixed at its
setpoint value at all times.
"""
from collections.abc import Iterator
import warnings
from dataclasses import dataclass
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid, CoolPropWarning
from hvac.cooling_load_calc_old import FixedTemperatureZone
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
    Q_cc:
        Cooling coil load at design conditions.
    SHR_cc:
        Sensible heat ratio of cooling coil load at design conditions.
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
    Q_cc: Quantity
    SHR_cc: Quantity

    def __post_init__(self):
        self.m_dot_vent = self.V_dot_vent_ntp * standard_air.rho
        self.m_dot_supply = self.V_dot_supply_ntp * standard_air.rho
        self.m_dot_recir = self.m_dot_supply - self.m_dot_vent
        self.V_dot_recir_ntp = self.m_dot_recir / standard_air.rho


class VAVSingleZoneAirCoolingSystem:
    """
    Implements the control strategy of a single-zone VAV air-cooling system:
    - supply fan control
    - mixing control
    - cooling control
    - heating control
    """
    def __init__(
        self,
        design_data: DesignData,
        T_set_cool: Quantity | None = None,
        W_set_cool: Quantity | None = None,
        rel_m_dot_supply_min: Quantity = Q_(60, 'pct'),
        has_rev_heatpump: bool = False,
        has_var_cooling_setpoint: bool = False,
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
            (Nominal) setpoint temperature of the air leaving the cooling coil. 
            If None, the supply air temperature in `design_data` will be used.
        W_set_cool: optional
            The humidity ratio of the air leaving the cooling coil that belongs
            with the setpoint temperature of the cooled air. If None, the supply
            air humidity ratio in `design_data` will be used.
        rel_m_dot_supply_min: default 60 %
            The minimum mass flow rate of supply air to ensure proper mixing
            with zone air expressed as a fraction of the supply air mass flow
            rate determined for summer peak design conditions.
        has_rev_heatpump: default False
            Indicates if the airco-system has a reversible heat pump instead
            of a separate cooling coil and a separate heating coil.
        has_var_cooling_setpoint: default False
            Indicates if the cooling coil/heatpump controller has a variable 
            setpoint for the air leaving the cooling coil, or not.
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
        self.has_rev_heatpump = has_rev_heatpump
        self.has_var_cooling_setpoint = has_var_cooling_setpoint

        # Set the nominal cooling coil operating point:
        self.T_set_cool = (
            T_set_cool if T_set_cool is not None
            else self.design_data.supply_air.Tdb
        )
        self.W_set_cool = (
            W_set_cool if W_set_cool is not None
            else self.design_data.supply_air.W
        )

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

        # Set the maximum flow rate of the supply air:
        self.m_dot_supply_max = self.design_data.m_dot_supply
        # The maximum supply air mass flow rate is limited to its design value,
        # which is considered to be the maximum flow rate that the supply fan can
        # displace.

        # Determine the ADP of the cooling coil at design conditions.
        # It will be assumed that the ADP of the cooling coil is constant
        # (unless the absolute humidity ratio of this ADP is higher than the
        # absolute humidity ratio of the inlet air --> see method
        # `_cooling_coil_control`).
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
            Q=self.design_data.Q_cc,
            SHR=self.design_data.SHR_cc
        )
        self.ADP_cc = cooling_coil_des.ADP

        # Set the desired measuring units of the output values:
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
        """Determines from the given state of supply air, the required supply
        air mass flow rate to maintain the zone air setpoint temperature.
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
            The state of supply air needed to maintain the zone air
            temperature.
        return_air:
            The state of the return air that depends on the sensible and
            latent cooling load of the zone, the supply air mass flow rate, and
            the state of the supply air.
        """
        # Calculate the mass flow rate of supply air needed to maintain the
        # zone air setpoint temperature with the cooling coil being operated at
        # its normal setpoint. The moisture content of the zone and return air
        # will depend on the moisture content of the supply air and on the
        # latent load of the zone.
        zone = AirConditioningProcess(
            air_in=supply_air,
            T_ao=self.T_zone_sp,
            Q_sen=self.Q_dot_zone_sen,
            Q_lat=self.Q_dot_zone_lat,
            h_w=Q_(0, 'J / kg')
        )
        m_dot_supply = zone.m_da
        return_air = zone.air_out

        # The mass flow rate of supply air cannot become smaller than its
        # minimum limit.
        if m_dot_supply < self.m_dot_supply_min:
            m_dot_supply = self.m_dot_supply_min
            # To maintain the zone air setpoint temperature, the supply air
            # needs to be heated (see `_control_heating_coil`), unless the
            # mixed air has a higher temperature than the required supply air
            # temperature.
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

        # The mass flow rate of supply air cannot become greater than its
        # maximum limit:
        if m_dot_supply > self.m_dot_supply_max:
            m_dot_supply = self.m_dot_supply_max
            # To maintain the zone air setpoint temperature, the supply air
            # will need to be cooled extra (see `_control_cooling_coil`).
            zone = AirConditioningProcess(
                W_ai=supply_air.W,
                # Actually, `W_ai` will also change because the cooling
                # process line is sloped.
                T_ao=self.T_zone_sp,
                m_da=m_dot_supply,
                Q_sen=self.Q_dot_zone_sen,
                Q_lat=self.Q_dot_zone_lat,
                h_w=Q_(0, 'J / kg')
            )
            supply_air = zone.air_in  # new required state of supply air
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
            # Ventilation air flow rate is limited to its minimum required
            # rate (determined when designing the system).
            # --> OA damper in its minimum open position
            # --> RA damper in its maximum open position
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
            # Outdoor air is mixed with return air to reach the cooling coil
            # setpoint temperature.
            # --> the mixing controller controls the position of the OA and RA
            # damper.
            mixing_chamber = AdiabaticMixing(
                in1=AirStream(state=return_air),
                in2=AirStream(state=self.outdoor_air),
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
            # It will cost less energy to cool the outdoor air to the cooling
            # coil temperature. So, no mixing takes place.
            # --> OA damper in its maximum open position
            # --> RA damper in its minimum open position.
            m_dot_vent = m_dot_supply
            m_dot_recir = Q_(0.0, 'kg / s')
            mixed_air = self.outdoor_air
        return mixed_air, m_dot_vent, m_dot_recir

    def _control_cooling_coil(
        self,
        m_dot_supply: Quantity,
        mixed_air: HumidAir,
        req_supply_air: HumidAir
    ) -> tuple[HumidAir, AirConditioningProcess | None]:
        """Determines the operation of the cooling coil.

        Parameters
        ----------
        m_dot_supply:
            The mass flow rate of supply air to the zone that flows through the
            cooling coil. This was determined in method `_control_supply_fan`.
        mixed_air:
            The state of air leaving the mixing chamber and entering the cooling
            coil. This was determined in method `_control_mixing_air`.
        req_supply_air:
            The state of supply air required to maintain the zone air setpoint
            temperature. This was determined in method `_control_supply_fan`.

        Returns
        -------
        cooled_air:
            The state of air leaving the cooling coil and entering the heating
            coil if it is present in the air-cooling system.
        cooling_coil:
            If mechanical cooling is possible, the cooling coil (an
            `AirConditioningProcess` object) is returned. Otherwise, None is
            returned.
        """
        if self.has_var_cooling_setpoint:
            # The cooling coil controller controls the leaving air temperature
            # to the required supply air temperature.
            T_set_cool = req_supply_air.Tdb
        else:
            # The cooling coil controller controls the leaving air temperature
            # to a fixed value.
            T_set_cool = self.T_set_cool
        if mixed_air.Tdb > T_set_cool:
            # Only if the mixed air temperature is higher than the leaving air
            # temperature, the cooling coil can cool the air.
            if self.ADP_cc.W > mixed_air.W:
                # The ADP humidity ratio cannot be higher than the humidity
                # ratio of the entering air, as this would mean that the air
                # is humidified in the cooling coil.
                ADP_cc = HumidAir(Tdb=mixed_air.Tdp, RH=Q_(100, 'pct'))
            else:
                ADP_cc = self.ADP_cc
            cooling_coil = AirConditioningProcess(
                air_in=mixed_air,
                m_da=m_dot_supply,
                T_ao=T_set_cool,
                ADP=ADP_cc
            )
            cooled_air = cooling_coil.air_out
        else:
            # The cooling coil is turned off.
            cooling_coil = None
            cooled_air = mixed_air
        return cooled_air, cooling_coil

    @staticmethod
    def _control_heating_coil(
        m_dot_supply: Quantity,
        cooled_air: HumidAir,
        supply_air_req: HumidAir,
    ) -> tuple[HumidAir, AirConditioningProcess | None]:
        """Determines the operation of the heating coil.

        Parameters
        ----------
        m_dot_supply:
            The mass flow rate of supply air to the zone that flows through the
            heating coil. This was determined in method `_control_supply fan`.
        cooled_air:
            The state of air leaving the cooling coil.
        supply_air_req:
            The state of supply air required to maintain the zone air setpoint
            temperature. This was determined in method `_control_supply fan()`.

        Returns
        -------
        supply_air:
            The state of air leaving the heating coil.
        heating_coil:
            If the heating coil needs to be active to reach the required supply
            air temperature, the heating coil (`AirConditioningProcess` object)
            is returned. Otherwise, None is returned.
        """
        if cooled_air.Tdb < supply_air_req.Tdb:
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
            # If the cooled air leaving the cooling coil doesn't need heating,
            # return the state of air entering the heating coil:
            supply_air = cooled_air
            heating_coil = None
        return supply_air, heating_coil

    def _control_reversible_heat_pump(
        self,
        m_dot_supply: Quantity,
        mixed_air: HumidAir,
        supply_air_req: HumidAir
    ) -> tuple[HumidAir, AirConditioningProcess | None]:
        if supply_air_req.Tdb > mixed_air.Tdb:
            # requires heating
            return self._control_heating_coil(m_dot_supply, mixed_air, supply_air_req)
        elif supply_air_req.Tdb < mixed_air.Tdb:
            # requires cooling
            return self._control_cooling_coil(m_dot_supply, mixed_air, supply_air_req)
        else:
            return HumidAir(Tdb=supply_air_req.Tdb, W=mixed_air.W), None

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
            cooling coil. This was determined in method `_control_supply fan`.
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

        # Initially assume that air being supplied to the zone has its nominal
        # state:
        supply_air = HumidAir(Tdb=self.T_set_cool, W=self.W_set_cool)

        # Determine:
        # (1) the supply air mass flow rate to maintain the zone air temperature
        #     at its setpoint,
        # (2) the resulting state of the return air, and
        # (3) the required, new state of the supply air, should the supply air
        #     mass flow rate be below its minimum or above its maximum limit.
        m_dot_supply, return_air, req_supply_air = self._control_supply_fan(supply_air)

        # Mixing control:
        mixed_air, m_dot_vent, m_dot_recir = self._control_mixing_air(
            m_dot_supply,
            return_air
        )

        if not self.has_rev_heatpump:
            # Cooling coil control:
            cooled_air, cooling_coil = self._control_cooling_coil(
                m_dot_supply,
                mixed_air,
                req_supply_air
            )
            # Heating coil control:
            supply_air, heating_coil = self._control_heating_coil(
                m_dot_supply,
                cooled_air,
                req_supply_air
            )
        else:
            cooled_air = mixed_air
            supply_air, heat_pump = self._control_reversible_heat_pump(
                m_dot_supply,
                mixed_air,
                req_supply_air
            )
            cooling_coil = heating_coil = heat_pump

        # Zone:
        return_air = self._zone_air(
            m_dot_supply,
            supply_air
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
            self.Q_dot_cc = -cooling_coil.Q if cooling_coil.Q < 0 else Q_(0.0, 'W')
            self.SHR_cc = cooling_coil.SHR
        if heating_coil is None:
            self.Q_dot_hc = Q_(0.0, 'W')
        else:
            self.Q_dot_hc = heating_coil.Q_sen if heating_coil.Q > 0 else Q_(0.0, 'W')
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
    """Helper class to prepare the cooling load data of the single-zone building
    needed to analyze the `VAVSingleZoneAirCoolingSystem` model on an hourly
    basis on a given day.
    """
    def __init__(
        self,
        zone: FixedTemperatureZone,
        design_data: DesignData
    ) -> None:
        """Creates an instance of class `CoolingSimData`.

        Parameters
        ----------
        zone:
            `ConditionedZone` object that represents the thermal model of the
            single-zone building, which is already initialized with the
            weather data for the day the simulation is to be run.
        design_data:
            `DesignData` object that contains data from the design
            calculations of the single-zone VAV air-cooling system, including
            the sensible and latent design cooling load.

        On instantiation of the `CoolingSimData` object the cooling load of the
        zone is calculated for each hour of the selected day. Also the relative
        sensible and latent cooling loads are calculated with respect to the
        design cooling load given in `design_data`.
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
        # Calculate the zone loads (heat gains) for each hour of the day:
        self.df_Q_zone = zone.solve(unit='kW', num_cycles=10)
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
        self.design_data = design_data
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

        These tuples will be used as arguments when calling method `analyze` of
        class `VAVSingleZoneAirCoolingSystem`.
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
        """Returns for the given hour (index 0 to 23) a tuple of four elements:
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
