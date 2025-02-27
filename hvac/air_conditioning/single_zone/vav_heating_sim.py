"""
Steady-state part-load analysis of a single-zone VAV air heating system
optionally equipped with a steam humidifier.
"""
import warnings
from dataclasses import dataclass
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid, CoolPropWarning
from hvac.air_conditioning import AirConditioningProcess, AirStream, AdiabaticMixing
from hvac.charts import PsychrometricChart, StatePoint
from .heating_design import Output
from hvac.logging import ModuleLogger

warnings.filterwarnings('ignore', category=CoolPropWarning)
logger = ModuleLogger.get_logger(__name__)
logger.setLevel(ModuleLogger.ERROR)

Q_ = Quantity


# Define standard (NTP) air:
Air = Fluid('Air')
standard_air = Air(
    T=Q_(20, 'degC'),
    P=Q_(101_325, 'Pa')
)
Water = Fluid('Water')


@dataclass
class DesignData:
    """Groups the design data of the single-zone VAV air-heating system.

    Attributes
    ----------
    Q_dot_zone_sen: Quantity
        Sensible zone load (heat loss) at design conditions.
    Q_dot_zone_lat: Quantity
        Latent zone load at design conditions.
    V_dot_vent_ntp:
        Minimum required volume flow rate of outdoor ventilation air referred
        to NTP air (dry air at 20 °C and 101.325 kPa).
    V_dot_supply_ntp:
        Required supply air volume flow rate at winter peak design conditions
        referred to NTP air.
    outdoor_air:
        State of outdoor air at design conditions.
    zone_air:
        Required state of zone air.
    supply_air:
        Required state of supply air at design conditions.
    T_steam:
        Temperature of injected steam for supply air humidification.
    Q_dot_zone_ihg:
        Total internal heat gain in the zone (assumed to be constant).
    ----------------------------------------------------------------------------
    m_dot_vent:
        Outdoor ventilation air mass flow rate. Internally calculated.
    m_dot_supply:
        Supply air mass flow rate. Internally calculated.
    m_dot_recir:
        Recirculation air mass flow rate. Internally calculated.
    V_dot_recir_ntp:
        Recirculation air volume flow rate. Internally calculated.
    K_zone:
        Global heat loss coefficient of the zone. Internally calculated
        and assumed to be a constant.
    T_outdoor_bal:
        Balance point temperature of the zone, i.e., the outdoor air temperature
        at which the sensible heat loss of the zone balances with the internal
        heat gain in the zone. Internally calculated and assumed to be constant.
    K_supply:
        Slope of the outdoor reset line to determine the supply air temperature
        as a function of the outdoor air temperature. Internally calculated.
    """
    Q_dot_zone_sen: Quantity
    Q_dot_zone_lat: Quantity
    V_dot_vent_ntp: Quantity
    V_dot_supply_ntp: Quantity
    outdoor_air: HumidAir
    zone_air: HumidAir
    supply_air: HumidAir
    T_steam: Quantity | None = None
    Q_dot_zone_ihg: Quantity = Q_(0, 'W')

    def __post_init__(self):
        self.m_dot_vent = self.V_dot_vent_ntp * standard_air.rho
        self.m_dot_supply = self.V_dot_supply_ntp * standard_air.rho
        self.m_dot_recir = self.m_dot_supply - self.m_dot_vent
        self.V_dot_recir_ntp = self.m_dot_recir / standard_air.rho
        self.K_zone = self.Q_dot_zone_sen / (self.outdoor_air.Tdb - self.zone_air.Tdb)
        self.T_outdoor_bal = self.zone_air.Tdb - self.Q_dot_zone_ihg / self.K_zone
        self.K_supply = (
            (self.supply_air.Tdb - self.zone_air.Tdb)
            / (self.outdoor_air.Tdb - self.T_outdoor_bal)
        )


class VAVSingleZoneAirHeatingSystem:
    """Model for analyzing steady-state part-load operation of a single-zone VAV
     air-heating system optionally equipped with a steam humidifier.
    """
    def __init__(
        self,
        design_data: DesignData,
        T_zone: Quantity | None = None,
        RH_zone: Quantity | None = None,
        humidifier_present: bool = True,
        T_steam: Quantity | None = None,
        units: dict[str, tuple[str, int]] | None = None
    ) -> None:
        """
        Creates a `VAVSingleZoneAirHeatingSystem` object.

        Parameters
        ----------
        design_data:
            `DesignData` instance that groups the design data of the AC system
            (see class `DesignData` in this module).
        T_zone: optional
            Setpoint of zone air temperature. If None, the zone air temperature
            in `design_data` will be used.
        RH_zone: optional
            Setpoint of zone air relative humidity. If None, the zone air
            relative humidity in `design_data` will be used.
        humidifier_present: default True
            Indicates of a steam humidifier is present downstream of the heating
            coil or not.
        T_steam: optional
            Temperature of saturated steam injected to humidify the air supplied
            to the zone. If None, the steam temperature in `design_data` will
            be used.
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
        self.T_zone = T_zone if T_zone is not None else self.design_data.zone_air.Tdb
        self.RH_zone = RH_zone if RH_zone is not None else self.design_data.zone_air.RH
        self.humidifier_present = humidifier_present
        self.T_steam = T_steam if T_steam is not None else self.design_data.T_steam
        self.units = units or {}

        self.zone_air = HumidAir(Tdb=self.T_zone, RH=self.RH_zone)
        self.W_zone = self.zone_air.W
        self.m_dot_vent = self.design_data.m_dot_vent

        if self.humidifier_present:
            if self.T_steam is not None:
                self.steam = Water(T=self.T_steam, x=Q_(1, 'frac'))
            else:
                raise ValueError('Humidifier present, but `T_steam` is `None`.')

        self.outdoor_air: HumidAir | None = None
        self.rel_Q_zone_sen: Quantity | None = None
        self.rel_Q_zone_lat: Quantity | None = None

        self.m_dot_supply: Quantity | None = None
        self.m_dot_recir: Quantity | None = None
        self.mixed_air: HumidAir | None = None
        self.heated_air: HumidAir | None = None
        self.supply_air: HumidAir | None = None
        self.return_air: HumidAir | None = None
        self.Q_dot_hc: Quantity | None = None
        self.m_dot_steam: Quantity | None = None

    def analyze(
        self,
        outdoor_air: HumidAir,
        rel_Q_zone_sen: Quantity | None,
        rel_Q_zone_lat: Quantity,
        i_max: int = 50,
        T_tol: Quantity = Q_(0.01, 'K'),
        W_tol: Quantity = Q_(0.01, 'g / kg')
    ) -> Output | None:
        """Analyzes VAV system heating operation at the external conditions
        given.

        The setpoint of supply air temperature is determined as a function of
        outdoor air temperature based on the outdoor reset line (outdoor reset
        control).

        Parameters
        ----------
        outdoor_air:
            State of outdoor air.
        rel_Q_zone_sen:
            Fraction of the sensible zone load with respect to the sensible
            zone load under winter peak design conditions. If set to `None`,
            the sensible zone load is internally calculated as a function of
            the given outdoor air temperature.
        rel_Q_zone_lat:
            Fraction of the latent zone load with respect to the latent zone
            load under winter peak design conditions.
        i_max: optional, default 50
            Maximum number of iterations to calculate the steady-state
            operation under the given conditions.
        T_tol: optional, default 0.01 K
            Tolerance for the supply air temperature.
        W_tol: optional, default 0.01 g/kg
            Tolerance for the supply air absolute humidity.

        Returns
        -------
        `Output` object (see class `Output` in module `heating_design`).
        """
        self.outdoor_air = outdoor_air
        if rel_Q_zone_sen is not None:
            self.rel_Q_zone_sen = rel_Q_zone_sen
        else:
            # If `rel_Q_zone_sen` was set to None, determine the sensible
            # heat loss of the zone as a function of outdoor air temperature:
            Q_zone_sen = (
                self.design_data.K_zone
                * (self.outdoor_air.Tdb - self.design_data.T_outdoor_bal)
            )
            self.rel_Q_zone_sen = Q_zone_sen / self.design_data.Q_dot_zone_sen
        self.rel_Q_zone_lat = rel_Q_zone_lat
        # The supply air temperature is determined by the outdoor reset line
        # of the zone:
        T_supply = (
            self.design_data.supply_air.Tdb
            - self.design_data.K_supply
            * (self.design_data.outdoor_air.Tdb - self.outdoor_air.Tdb)
        )
        # With the supply air temperature being determined, the mass flow rate
        # of supply air can be calculated from a sensible heat balance:
        zone = AirConditioningProcess(
            T_ai=T_supply,
            T_ao=self.T_zone,
            Q_sen=self.rel_Q_zone_sen * self.design_data.Q_dot_zone_sen
        )
        # However, the supply air mass flow rate cannot be greater than its
        # design value (the size of the supply fan is fixed), and may not become
        # less than 60 % of its design value (for proper air mixing in the zone)
        # or less than the minimum required ventilation air mass flow rate:
        m_dot_supply = max(
            min(zone.m_da, self.design_data.m_dot_supply),
            self.m_dot_vent,
            0.60 * self.design_data.m_dot_supply
        )
        # Recirculation mass flow rate:
        m_dot_recir = m_dot_supply - self.m_dot_vent
        # We initially guess that the desired humidity level in the zone is
        # maintained, and with this assumption we determine the initial humidity
        # of the supply air:
        zone = AirConditioningProcess(
            T_ai=T_supply,
            m_da=m_dot_supply,
            W_ao=self.zone_air.W,
            Q_sen=self.rel_Q_zone_sen * self.design_data.Q_dot_zone_sen,
            Q_lat=self.rel_Q_zone_lat * self.design_data.Q_dot_zone_lat
        )
        supply_air = zone.air_in
        return_air = zone.air_out
        for i in range(i_max):
            logger.debug(
                "return air: "
                f"{return_air.Tdb.to('degC'):~P.2f} DB, "
                f"{return_air.W.to('g / kg'):~P.2f} AH "
                f"({return_air.RH.to('pct'):~P.2f} RH)"
            )
            # Mixing return air with outdoor air:
            mixing_chamber = AdiabaticMixing(
                in1=AirStream(
                    state=return_air,
                    m_da=m_dot_recir
                ),
                in2=AirStream(
                    state=self.outdoor_air,
                    m_da=self.m_dot_vent
                ),
                out=AirStream(m_da=m_dot_supply)
            )
            mixed_air = mixing_chamber.stream_out.state
            if mixed_air.Tdb >= supply_air.Tdb:
                # No heating needed.
                heated_air = mixed_air
                preheat_coil = None
            else:
                heated_air = HumidAir(Tdb=supply_air.Tdb, W=mixed_air.W)
                preheat_coil = AirConditioningProcess(
                    air_in=mixed_air,
                    air_out=heated_air,
                    m_da=m_dot_supply
                )
            if mixed_air.W >= supply_air.W or not self.humidifier_present:
                # Humidification is not needed or possible.
                humidified_air = heated_air
                humidifier = None
            else:
                humidified_air = supply_air
                humidifier = AirConditioningProcess(
                    air_in=heated_air,
                    air_out=humidified_air,
                    m_da=m_dot_supply,
                    Q=Q_(0, 'W'),
                    h_w=self.steam.h
                )
            # The state of humidified air must be equal to the state of supply
            # air:
            T_dev = abs((humidified_air.Tdb - supply_air.Tdb).to('K'))
            W_dev = abs((humidified_air.W - supply_air.W).to('g / kg'))
            if T_dev < T_tol.to('K') and W_dev < W_tol.to('g / kg'):
                self.m_dot_supply = m_dot_supply
                self.m_dot_recir = m_dot_recir
                self.mixed_air = mixed_air
                self.heated_air = heated_air
                self.supply_air = humidified_air
                self.return_air = return_air
                if preheat_coil is None:
                    self.Q_dot_hc = Q_(0, 'W')
                else:
                    self.Q_dot_hc = preheat_coil.Q_sen
                if humidifier is None:
                    self.m_dot_steam = Q_(0.0, 'kg / s')
                else:
                    self.m_dot_steam = humidifier.m_w
                return Output(
                    self.m_dot_supply, self.m_dot_vent, self.m_dot_recir,
                    self.mixed_air, self.heated_air, self.supply_air, self.return_air,
                    self.Q_dot_hc, self.m_dot_steam,
                    units=self.units
                )
            supply_air = humidified_air
            zone = AirConditioningProcess(
                air_in=supply_air,
                m_da=m_dot_supply,
                Q_sen=self.rel_Q_zone_sen * self.design_data.Q_dot_zone_sen,
                Q_lat=self.rel_Q_zone_lat * self.design_data.Q_dot_zone_lat
            )
            return_air = zone.air_out

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
            name='air heating',
            start_point=StatePoint(
                self.mixed_air.Tdb,
                self.mixed_air.W
            ),
            end_point=StatePoint(
                self.heated_air.Tdb,
                self.heated_air.W
            )
        )
        psy_chart.plot_process(
            name='air humidification',
            start_point=StatePoint(
                self.heated_air.Tdb,
                self.heated_air.W
            ),
            end_point=StatePoint(
                self.supply_air.Tdb,
                self.supply_air.W
            )
        )
        psy_chart.plot_process(
            name='zone',
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
