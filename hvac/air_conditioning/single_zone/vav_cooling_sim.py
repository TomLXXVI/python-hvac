"""
Steady-state part-load analysis of a single-zone VAV air-cooling system equipped
with a DX air-cooling coil and optionally an air-heating coil.
"""
from collections.abc import Iterator
import warnings
from dataclasses import dataclass
import pandas as pd
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid, CoolPropWarning, CP_HUMID_AIR
import hvac.heat_transfer.heat_exchanger.fin_tube.core as core
from hvac.air_conditioning import AirConditioningProcess, AirStream, AdiabaticMixing
from hvac.charts.psychrometric_chart import PsychrometricChart, StatePoint
from .cooling_design import Output
from hvac.logging import ModuleLogger


warnings.filterwarnings('ignore', category=CoolPropWarning)


Q_ = Quantity
HexCore = core.plain_fin_tube.PlainFinTubeHeatExchangerCore
logger = ModuleLogger.get_logger(__name__)


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
    """Model for analyzing the steady-state part-load operation of a single-zone
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
    outdoor air flow rate is restricted to what's required for zone ventilation,
    only and recirculation of return air is now at its maximum.
    """
    def __init__(
        self,
        design_data: DesignData,
        T_set_cool: Quantity | None = None,
        W_min_cool: Quantity | None = None,
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
            Setpoint of air temperature at the cooling coil outlet. If None, the
            supply air temperature in `design_data` will be used.
        W_min_cool: optional
            The minimum humidity ratio of the cooled air that can be attained
            at the set temperature of the cooled air. If None, the supply air
            humidity ratio in `design_data` will be used.
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

        self.T_set_cool = (
            T_set_cool if T_set_cool is not None
            else self.design_data.supply_air.Tdb
        )
        self.W_min_cool = (
            W_min_cool if W_min_cool is not None
            else self.design_data.supply_air.W
        )

        self.heating_coil_present = heating_coil_present

        if not self.heating_coil_present and variable_cooling_setpoint:
            self.variable_cooling_setpoint = True
        else:
            self.variable_cooling_setpoint = False

        self.rel_m_dot_supply_min = rel_m_dot_supply_min

        self.units = units or {}

        self.outdoor_air: HumidAir | None = None
        self.rel_Q_zone_sen: Quantity | None = None
        self.rel_Q_zone_lat: Quantity | None = None
        self.T_set_zone: Quantity | None = None
        self.release_cooling: bool = True

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

        # The minimum supply air mass flow rate must not be less than the
        # minimum ventilation air requirement or not be less than 60 % of the
        # design value to ensure proper mixing of supply air with zone air:
        self.m_dot_supply_min = max(
            self.design_data.m_dot_vent,
            rel_m_dot_supply_min * self.design_data.m_dot_supply
        )
        # The maximum supply air mass flow rate is limited to its design value
        # (considered to be the maximum flow rate that the supply fan can
        # displace).
        self.m_dot_supply_max = self.design_data.m_dot_supply

    def _control_supply_fan(
        self,
        supply_air: HumidAir,
    ) -> tuple[Quantity, HumidAir, HumidAir]:
        # Calculate the mass flow rate of supply air needed to maintain the
        # setpoint zone air temperature with the cooling coil being active.
        # The moisture content of the zone and return air will depend on the
        # moisture content of the supply air and on the latent load of the zone.
        zone = AirConditioningProcess(
            air_in=supply_air,
            T_ao=self.T_set_zone,
            Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
            Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
        )
        m_dot_supply = zone.m_da
        return_air = zone.air_out

        if m_dot_supply < self.m_dot_supply_min:
            # The mass flow rate of supply air cannot become less than the
            # minimum limit:
            m_dot_supply = self.m_dot_supply_min

            if self.heating_coil_present:
                # For the setpoint zone air temperature to be maintained, the
                # supply air needs to be heated (see `_control_heating_coil`).
                # Heating the supply-air does not change its moisture content,
                # only its temperature is raised.
                zone = AirConditioningProcess(
                    W_ai=supply_air.W,
                    T_ao=self.T_set_zone,
                    m_da=m_dot_supply,
                    Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
                    Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
                )
                supply_air = zone.air_in  # new state of supply air
                return_air = zone.air_out
            else:
                # If no heating coil is present:
                if self.variable_cooling_setpoint:
                    # The cooling coil controller is able to increase the cooling
                    # air setpoint temperature to maintain the zone air setpoint
                    # temperature.
                    # --> To what extent is this possible?
                    zone = AirConditioningProcess(
                        W_ai=supply_air.W,
                        # Actually, `W_ai` will also change because the cooling
                        # process line is sloped.
                        T_ao=self.T_set_zone,
                        m_da=m_dot_supply,
                        Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
                        Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
                    )
                    supply_air = zone.air_in  # new state of supply air
                    return_air = zone.air_out
                else:
                    # The setpoint zone air temperature cannot be maintained.
                    # Use outdoor air as supply air (and turn off the cooling
                    # coil --> see `_control_cooling_coil(...)`).
                    supply_air = self.outdoor_air  # new state of supply air
                    zone = AirConditioningProcess(
                        air_in=supply_air,
                        m_da=m_dot_supply,
                        Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
                        Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
                    )
                    return_air = zone.air_out
                    try:
                        _ = return_air.RH
                    except ValueError:
                        return_air = HumidAir(Tdb=return_air.Tdb, RH=Q_(100, 'pct'))

        if m_dot_supply > self.m_dot_supply_max:
            # The mass flow rate of supply air cannot become greater than the
            # maximum limit:
            m_dot_supply = self.m_dot_supply_max

            if self.variable_cooling_setpoint:
                # The cooling coil controller is able to decrease the cooling
                # air setpoint temperature to maintain the zone air setpoint
                # temperature.
                # --> To what extent is this possible?
                zone = AirConditioningProcess(
                    W_ai=supply_air.W,
                    # Actually, `W_ai` will also change because the cooling
                    # process line is sloped.
                    T_ao=self.T_set_zone,
                    m_da=m_dot_supply,
                    Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
                    Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
                )
                supply_air = zone.air_in  # new state of supply air
                return_air = zone.air_out
            else:
                # The zone air temperature will inevitably raise above its
                # setpoint.
                zone = AirConditioningProcess(
                    air_in=supply_air,
                    m_da=m_dot_supply,
                    Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
                    Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
                )
                return_air = zone.air_out

        return m_dot_supply, return_air, supply_air

    def _control_mixing_air(
        self,
        m_dot_supply: Quantity,
        return_air: HumidAir
    ) -> tuple[HumidAir, Quantity, Quantity]:
        # Mixing control. Determine the state of mixed air, the mass flow rate
        # of outdoor ventilation air and recirculated return air.
        if self.outdoor_air.Tdb >= self.T_set_zone:
            # If the outdoor air temperature is higher than the zone air
            # setpoint temperature: outdoor ventilation air flow rate is
            # restricted to the minimum ventilation rate required (determined
            # when designing the system) --> OA-damper in its minimum open
            # position (NC) and RA-damper in its maximum open position (NO).
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
        elif self.outdoor_air.Tdb < self.T_set_cool:
            # No mechanical cooling is needed.
            # If outdoor air temperature is less than the setpoint supply
            # air temperature: outdoor ventilation air is mixed with return
            # air to get air with a temperature equal to the setpoint
            # supply air temperature --> air mixing controller controls
            # the position of the OA- and RA-damper.
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
            # If outdoor air temperature is between the setpoint supply air
            # temperature and the zone air setpoint temperature: mixing cannot
            # establish the setpoint supply air temperature, because mixing
            # outdoor air with return air would inevitably raise the mixed air
            # temperature above the outdoor air temperature.
            # It will cost less cooling energy to cool the outdoor air to the
            # setpoint supply air temperature.
            # So, no mixing takes place --> OA-damper in its maximum open
            # position and RA-damper in its minimum open position.
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
        # Cooling coil control. Determine if the cooling coil can be active and
        # the state of air at the cooling coil outlet.
        if (m_dot_supply == self.m_dot_supply_min) and not self.heating_coil_present:
            # If the mass flow rate of supply air is at its minimum limit, and
            # no heating coil is present, the zone air temperature will drop
            # below its setpoint if the cooling coil is active, unless the
            # cooling coil controller can increase the cooling air setpoint
            # temperature to get at the required state of supply air needed to
            # maintain the zone air temperature setpoint.
            if self.variable_cooling_setpoint:
                if mixed_air.Tdb > supply_air_req.Tdb:
                    # Only if mixed air has a higher temperature than the
                    # required supply air temperature, the mixed air can be
                    # cooled to the required supply air temperature.
                    cooling_coil = AirConditioningProcess(
                        air_in=mixed_air,
                        m_da=m_dot_supply,
                        air_out=supply_air_req
                    )
                    cooled_air = supply_air_req
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
            # cooling coil controller can decrease the cooling air setpoint
            # temperature to get at the required state of supply air needed to
            # maintain the zone air temperature setpoint.
            if mixed_air.Tdb > supply_air_req.Tdb:
                cooling_coil = AirConditioningProcess(
                    air_in=mixed_air,
                    m_da=m_dot_supply,
                    air_out=supply_air_req
                )
                cooled_air = supply_air_req
            else:
                cooling_coil = None
                cooled_air = mixed_air
        elif mixed_air.Tdb > self.T_set_cool:
            # If the mass flow rate of supply air is between its minimum and
            # maximum limit, and the mixed air temperature is higher than the
            # cooling air setpoint temperature, the cooling coil controller
            # keeps the cooled air state fixed at its normal setpoint (both
            # temperature `T_set_cool` and humidity ratio `W_min_cool`).
            cooled_air = HumidAir(
                Tdb=self.T_set_cool,
                W=self.W_min_cool
            )
            cooling_coil = AirConditioningProcess(
                air_in=mixed_air,
                air_out=cooled_air,
                m_da=m_dot_supply
            )
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
        # Heating coil control.
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
            supply_air = cooled_air
            heating_coil = None
        return supply_air, heating_coil

    def _zone_air(
        self,
        m_dot_supply: Quantity,
        supply_air: HumidAir
    ) -> HumidAir:
        #
        zone = AirConditioningProcess(
            air_in=supply_air,
            m_da=m_dot_supply,
            Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
            Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
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
        T_set_zone: Quantity | None = None
    ) -> Output:
        """Analyzes the operating state of the system for the external
        conditions given.

        Parameters
        ----------
        outdoor_air:
            Outdoor air state.
        rel_Q_zone_sen:
            Fraction of the sensible zone load at design conditions (either
            a fraction between 0 and 1, or a percentage between 0 and 100 %).
        rel_Q_zone_lat:
            Fraction of the latent zone load at design conditions (either
            a fraction between 0 and 1, or a percentage between 0 and 100 %).
        T_set_zone: optional
            Setpoint of zone air temperature. If None, the zone air temperature
            in `design_data` will be used.

        Returns
        -------
        `Output` object (see class `Output` in module `cooling_design`).
        """
        self.outdoor_air = outdoor_air
        self.rel_Q_zone_sen = rel_Q_zone_sen
        self.rel_Q_zone_lat = rel_Q_zone_lat
        self.T_set_zone = (
            T_set_zone if T_set_zone is not None
            else self.design_data.zone_air.Tdb
        )

        # The cooling coil controller keeps the temperature of the air that
        # leaves the cooling coil at its setpoint if the cooling coil is active.
        # We further assume that the humidity ratio of the cooled air remains
        # fixed at `W_min_cool`.
        supply_air_ini = HumidAir(
            Tdb=self.T_set_cool,
            W=self.W_min_cool
        )

        # With the initial supply air state being selected, determine the needed
        # mass flow rate of supply air to maintain the zone air temperature at
        # its setpoint, the resulting state of return air, and the required
        # new state of supply air, should the needed mass flow rate of supply air
        # could not be attained with the initial state of supply air (i.e., if
        # the supply air mass flow rate is at its minimum or maximum limit).
        m_dot_supply, return_air, supply_air_req = self._control_supply_fan(
            supply_air_ini,
        )

        logger.debug(
            "return air 1: "
            f"{return_air.Tdb.to('degC'):~P.3f} DB, "
            f"{return_air.W.to('g / kg'):~P.3f} AH "
            f"({return_air.RH.to('pct'):~P.3f} RH)"
        )

        # Economizer operation: Determine the mass flow rate of outdoor
        # ventilation air and recirculated return air, and determine the state
        # of mixed air at the cooling coil inlet.
        mixed_air, m_dot_vent, m_dot_recir = self._control_mixing_air(
            m_dot_supply,
            return_air
        )

        per_m_dot_recir = m_dot_recir.to('kg / hr') / m_dot_supply.to('kg / hr')
        per_m_dot_vent = m_dot_vent.to('kg / hr') / m_dot_supply.to('kg / hr')

        logger.debug(
            "percentage mass flow rate recirculated return air: "
            f"{per_m_dot_recir.to('pct'):~P.3f}"
        )
        logger.debug(
            "percentage mass flow rate outdoor ventilation air: "
            f"{per_m_dot_vent.to('pct'):~P.3f}"
        )
        logger.debug(
            "mixed air: "
            f"{mixed_air.Tdb.to('degC'):~P.3f} DB, "
            f"{mixed_air.W.to('g / kg'):~P.3f} AH "
            f"({mixed_air.RH.to('pct'):~P.3f} RH)"
        )

        # Cooling coil.
        cooled_air, cooling_coil = self._control_cooling_coil(
            m_dot_supply,
            mixed_air,
            supply_air_req
        )

        logger.debug(
            "cooled air: "
            f"{cooled_air.Tdb.to('degC'):~P.3f} DB, "
            f"{cooled_air.W.to('g / kg'):~P.3f} AH "
            f"({cooled_air.RH.to('pct'):~P.3f} RH)"
        )

        # Heating coil.
        supply_air, heating_coil = self._control_heating_coil(
            m_dot_supply,
            cooled_air,
            supply_air_req
        )

        logger.debug(
            "supply air: "
            f"{supply_air.Tdb.to('degC'):~P.3f} DB, "
            f"{supply_air.W.to('g / kg'):~P.3f} AH "
            f"({supply_air.RH.to('pct'):~P.3f} RH)"
        )

        return_air = self._zone_air(
            m_dot_supply,
            supply_air
        )

        logger.debug(
            "return air 2: "
            f"{return_air.Tdb.to('degC'):~P.3f} DB, "
            f"{return_air.W.to('g / kg'):~P.3f} AH "
            f"({return_air.RH.to('pct'):~P.3f} RH)"
        )

        self.m_dot_supply = m_dot_supply
        self.m_dot_vent = m_dot_vent
        self.m_dot_recir = m_dot_recir
        self.cooled_air = cooled_air
        self.supply_air = supply_air
        self.return_air = return_air
        self.mixed_air = mixed_air
        if cooling_coil is None:
            self.Q_dot_cc = Q_(0.0, 'W')
            self.SHR_cc = Q_(float('nan'), 'frac')
        else:
            self.Q_dot_cc = -cooling_coil.Q
            Q_dot_cc_sen = (
                self.m_dot_supply * CP_HUMID_AIR
                * (self.mixed_air.Tdb - self.cooled_air.Tdb)
            )
            self.SHR_cc = Q_dot_cc_sen.to('W') / self.Q_dot_cc.to('W')
        if heating_coil is None:
            self.Q_dot_hc = Q_(0.0, 'W')
        else:
            self.Q_dot_hc = heating_coil.Q_sen

        return Output(
            self.m_dot_supply, self.m_dot_vent, self.m_dot_recir,
            self.mixed_air, self.cooled_air, self.supply_air, self.return_air,
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
        T_outdoor_db_rng: list[Quantity],
        T_outdoor_wb_rng: list[Quantity],
        df_Q_zone: pd.DataFrame,
        design_data: DesignData
    ) -> None:
        """Creates an instance of class `CoolingSimData`.

        Parameters
        ----------
        T_outdoor_db_rng:
            List of hourly dry-bulb outdoor air temperatures ordered from 0 to 23 h.
        T_outdoor_wb_rng:
            List of hourly wet-bulb outdoor air temperatures ordered from 0 to 23 h.
        df_Q_zone:
            Pandas-DataFrame object retrieved from a `Building`, a `BuildingEntity`,
            a `VentilationZone`, or a `Space` object in subpackage
            `hvac.cooling_load_calc` that contains the hourly sensible and latent
            cooling loads.
        design_data:
            `DesignData` object that contains data coming from the design calculations
            of the single-zone VAV air-cooling system, including the sensible and
            latent design cooling load.
        """
        self.T_outdoor_db_rng = T_outdoor_db_rng
        self.T_outdoor_wb_rng = T_outdoor_wb_rng
        self.df_dQ_zone = df_Q_zone
        self.design_data = design_data

        self.outdoor_air_rng = [
            HumidAir(Tdb=T_outdoor_db, Twb=T_outdoor_wb)
            for T_outdoor_db, T_outdoor_wb in zip(T_outdoor_db_rng, T_outdoor_wb_rng)
        ]

        self.Q_zone_sen_rng = [
            Q_(Q, 'kW')
            for Q in df_Q_zone['Q_sen_load'].values
        ]
        self.Q_zone_lat_rng = [
            Q_(Q, 'kW')
            for Q in df_Q_zone['Q_lat_load'].values
        ]

        self.rel_Q_zone_sen_rng = [
            Q_zone_sen.to('kW') / design_data.Q_zone_sen.to('kW')
            for Q_zone_sen in self.Q_zone_sen_rng
        ]
        self.rel_Q_zone_lat_rng = [
            Q_zone_lat.to('kW') / design_data.Q_zone_lat.to('kW')
            for Q_zone_lat in self.Q_zone_lat_rng
        ]

    def __call__(self) -> Iterator[list[tuple[HumidAir, Quantity, Quantity]]]:
        """Returns an iterator over tuples which have three elements:
        -   the state of outdoor air
        -   the corresponding relative sensible cooling load (ratio of the
            sensible part-load to the sensible design-load)
        -   the corresponding relative latent cooling load (ratio of the latent
            part-load to the latent design-load)

        These tuples can be used as arguments in the method `analyze` of class
        `VAVSingleZoneAirCoolingSystem`.
        """
        return zip(
            self.outdoor_air_rng,
            self.rel_Q_zone_sen_rng,
            self.rel_Q_zone_lat_rng
        )

    def get_sim_data(
        self,
        hour: int,
        is_relative: bool = True
    ) -> tuple[HumidAir, Quantity, Quantity]:
        """Returns for the given hour a tuple of three elements:
        -   the state of outdoor air
        -   the corresponding relative sensible cooling load (ratio of the
            sensible part-load to the sensible design-load)
        -   the corresponding relative latent cooling load (ratio of the latent
            part-load to the latent design-load)
        """
        outdoor_air = self.outdoor_air_rng[hour]
        if is_relative:
            rel_Q_zone_sen = self.rel_Q_zone_sen_rng[hour]
            rel_Q_zone_lat = self.rel_Q_zone_lat_rng[hour]
            return outdoor_air, rel_Q_zone_sen, rel_Q_zone_lat
        else:
            Q_zone_sen = self.Q_zone_sen_rng[hour]
            Q_zone_lat = self.Q_zone_lat_rng[hour]
            return outdoor_air, Q_zone_sen, Q_zone_lat
