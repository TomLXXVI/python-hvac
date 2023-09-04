from dataclasses import dataclass, field
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid
from hvac.air_conditioning import (
    AirConditioningProcess,
    AirStream,
    AdiabaticMixing
)
from hvac.charts import PsychrometricChart, StatePoint


Q_ = Quantity
Water = Fluid('Water')
Air = Fluid('Air')
standard_air = Air(T=Q_(20, 'degC'), P=Q_(101_325, 'Pa'))


@dataclass
class Output:
    m_dot_supply: Quantity
    m_dot_vent: Quantity
    m_dot_recir: Quantity
    mixed_air: HumidAir
    heated_air: HumidAir
    supply_air: HumidAir
    return_air: HumidAir
    Q_dot_hc: Quantity
    m_dot_steam: Quantity
    units: dict[str, tuple[str, int]] = field(default_factory=dict)

    def __post_init__(self):
        self.V_dot_vent_ntp = self.m_dot_vent / standard_air.rho
        self.V_dot_supply_ntp = self.m_dot_supply / standard_air.rho
        self.V_dot_recir_ntp = self.m_dot_recir / standard_air.rho
        self.default_units = {
            'm_dot': ('kg / hr', 3),
            'V_dot': ('m ** 3 / hr', 3),
            'T': ('degC', 3),
            'W': ('g / kg', 3),
            'RH': ('pct', 3),
            'Q_dot': ('kW', 3),
            'SHR': ('frac', 3)
        }
        self.default_units.update(self.units)
        self.units = self.default_units

    def __str__(self):
        return (
            "supply air mass flow rate = "
            f"{self.m_dot_supply.to(self.units['m_dot'][0]):~P.{self.units['m_dot'][1]}f}\n"
            "ventilation air mass flow rate = "
            f"{self.m_dot_vent.to(self.units['m_dot'][0]):~P.{self.units['m_dot'][1]}f}\n"
            "recirculation air mass flow rate = "
            f"{self.m_dot_recir.to(self.units['m_dot'][0]):~P.{self.units['m_dot'][1]}f}\n"
            "supply air volume flow rate (NTP) = "
            f"{self.V_dot_supply_ntp.to(self.units['V_dot'][0]):~P.{self.units['V_dot'][1]}f}\n"
            "ventilation air volume flow rate (NTP) = "
            f"{self.V_dot_vent_ntp.to(self.units['V_dot'][0]):~P.{self.units['V_dot'][1]}f}\n"
            "recirculation air volume flow rate (NTP) = "
            f"{self.V_dot_recir_ntp.to(self.units['V_dot'][0]):~P.{self.units['V_dot'][1]}f}\n"
            "mixed air = "
            f"{self.mixed_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.mixed_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.mixed_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "heated air = "
            f"{self.heated_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.heated_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.heated_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "supply air = "
            f"{self.supply_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.supply_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.supply_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "return air = "
            f"{self.return_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.return_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.return_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "heating coil load = "
            f"{self.Q_dot_hc.to(self.units['Q_dot'][0]):~P.{self.units['Q_dot'][1]}f}\n"
            "steam mass flow rate = "
            f"{self.m_dot_steam.to(self.units['m_dot'][0]):~P.{self.units['m_dot'][1]}f}"
        )


class AircoSystem:
    """Class to design an AC system (CAV or VAV) for heating and humidification
    of the supply air to a single-zone building in order to meet the sensible
    and latent zone load. Its purpose is to give the results that are needed for
    sizing and selecting the main components of the AC system.

    Attributes
    ----------
    m_dot_supply:
        The required mass flow rate of supply air to the zone.
    supply_air:
        The state of supply air under winter peak design conditions.
    m_dot_steam:
        The required saturated steam mass flow rate for humidification of the
        supply air to the zone.
    Q_dot_hc:
        The required heating capacity of the preheat-coil to heat the supply air
        to the zone. For winter peak design, it is assumed that air-heating only
        takes place in the preheat-coil. An additional heating coil is not being
        considered under winter peak design conditions.
    heated_air:
        The state of air at the preheat-coil outlet and entering the humidifier
        under peak design conditions.
    mixed_air:
        The state of air after mixing outdoor air with return air from the zone
        and entering the preheat-coil under winter peak design conditions.
    m_dot_recir:
        The mass flow rate of recirculated air under winter peak design
        conditions.
    m_dot_vent:
        The mass flow rate of outdoor ventilation air according to the minimum
        ventilation requirements for the zone.
    V_dot_supply_ntp:
        The volume flow rate of supply air referred to NTP conditions.
    V_dot_vent:
        The outdoor ventilation air volume flow rate referred to the actual
        state of outdoor air.
    """
    def __init__(
        self,
        zone_air: HumidAir,
        outdoor_air: HumidAir,
        Q_dot_zone: Quantity,
        SHR_zone: Quantity,
        V_dot_vent_ntp: Quantity,
        m_dot_supply_cooling: Quantity | None,
        T_supply_max: Quantity = Q_(40, 'degC'),
        T_steam: Quantity = Q_(100, 'degC'),
        units: dict[str, tuple[str, int]] | None = None
    ) -> None:
        """Creates an instance of `AircoSystem`, while performing the necessary
        design calculations based on the parameters given.

        Parameters
        ----------
        zone_air:
            The desired state of the air in the zone.
        outdoor_air:
            State of outdoor air for winter peak design.
        Q_dot_zone:
            The total peak heating load of the zone.
        SHR_zone:
            The sensible heat ratio of the total heating load under peak
            design conditions.
        V_dot_vent_ntp:
            The minimum required outdoor ventilation air volume flow rate
            referred to NTP conditions (dry air at 20 °C and 101.325 kPa).
        m_dot_supply_cooling: optional
            The design mass flow rate of supply air for summer peak design
            (cooling). If this mass flow rate should be greater than the
            required supply air mass flow rate calculated for winter peak
            design, this value will be used for the calculations (since the same
            supply fan will be used as well for summer as for winter operation).
            If None, the calculated value for winter peak design will be used.
        T_supply_max: default 40 °C
            The maximum allowable supply air temperature to avoid air
            stratification in the zone.
        T_steam: default 100 °C
            The temperature of saturated steam injected in the supply air stream
            to the zone.
        units: optional
            Units to be used for displaying the results. Keys are:
            'm_dot' for mass flow rate, 'V_dot' for volume flow rate, 'Q_dot'
            for heat transfer rate, 'T' for temperature, 'W' for absolute
            humidity, and 'RH' for relative humidity. Values are 2-tuples:
            the first element is the unit (a string), the second element is
            the number of decimals to be displayed (an int).
        """
        self.zone_air = self.return_air = zone_air
        self.outdoor_air = outdoor_air
        self.Q_dot_zone = Q_dot_zone
        self.SHR_zone = SHR_zone
        self.V_dot_vent_ntp = V_dot_vent_ntp
        self.m_dot_supply_cooling = m_dot_supply_cooling
        self.T_supply_max = T_supply_max
        self.T_steam = T_steam
        self.units = units or {}

        self.m_dot_vent: Quantity | None = None
        self.supply_air: HumidAir | None = None
        self.m_dot_supply: Quantity | None = None
        self.m_dot_recir: Quantity | None = None
        self.mixed_air: HumidAir | None = None
        self.heated_air: HumidAir | None = None
        self.m_dot_steam: Quantity | None = None
        self.Q_dot_hc: Quantity | None = None
        self.V_dot_supply_ntp: Quantity | None = None
        self.V_dot_vent: Quantity | None = None

    def design(self) -> Output:
        """Runs the design calculations and returns the results in an `Output`
        instance that holds:
        m_dot_supply:
            The required mass flow rate of the supply air to the zone.
        m_dot_vent:
            The mass flow rate of outdoor ventilation air.
        m_dot_recir:
            The mass flow rate of recirculated air.
        mixed_air:
            The state of air downstream of the mixing plenum and at the cooling
            coil inlet.
        heated_air:
            The required state of air at the heating coil outlet.
        supply_air:
            The required state of the supply air to compensate for the zone
            cooling loads.
        return_air:
            The state of air returned to the air handler. This is also the state
            of the zone air.
        Q_dot_hc:
            The required heating coil capacity at design conditions.
        m_dot_steam:
            The mass flow rate of saturated steam injected in the supply air.
        """
        self.m_dot_vent = self.V_dot_vent_ntp * standard_air.rho
        self.supply_air, self.m_dot_supply = self._determine_supply_air()
        self.m_dot_recir = self.m_dot_supply - self.m_dot_vent
        self.mixed_air = self._determine_mixed_air()
        self.heated_air, self.m_dot_steam = self._determine_heated_air()
        self.Q_dot_hc = self._determine_heating_coil_load()
        self.V_dot_supply_ntp = self.m_dot_supply / standard_air.rho
        self.V_dot_vent = self.m_dot_vent / self.outdoor_air.rho
        return Output(
            self.m_dot_supply, self.m_dot_vent, self.m_dot_recir,
            self.mixed_air, self.heated_air, self.supply_air, self.return_air,
            self.Q_dot_hc, self.m_dot_steam, self.units
        )

    def _determine_supply_air(self) -> tuple[HumidAir, Quantity]:
        zone = AirConditioningProcess(
            T_ai=self.T_supply_max,
            air_out=self.zone_air,
            Q=self.Q_dot_zone,
            SHR=self.SHR_zone
        )
        if zone.m_da > self.m_dot_supply_cooling:
            m_dot_supply = zone.m_da
            supply_air = zone.air_in
        else:
            m_dot_supply = self.m_dot_supply_cooling
            zone = AirConditioningProcess(
                air_out=self.zone_air,
                Q=self.Q_dot_zone,
                SHR=self.SHR_zone,
                m_da=m_dot_supply
            )
            supply_air = zone.air_in
        return supply_air, m_dot_supply

    def _determine_mixed_air(self) -> HumidAir:
        mixing_chamber = AdiabaticMixing(
            in1=AirStream(
                state=self.return_air,
                m_da=self.m_dot_recir
            ),
            in2=AirStream(
                state=self.outdoor_air,
                m_da=self.m_dot_vent
            ),
            out=AirStream(m_da=self.m_dot_supply)
        )
        return mixing_chamber.stream_out.state

    def _determine_heated_air(self) -> tuple[HumidAir, Quantity]:
        steam_sat = Water(T=self.T_steam, x=Q_(1, 'frac'))
        humidification = AirConditioningProcess(
            W_ai=self.mixed_air.W,
            air_out=self.supply_air,
            m_da=self.m_dot_supply,
            Q=Q_(0, 'W'),
            h_w=steam_sat.h
        )
        preheated_air = humidification.air_in
        m_dot_steam = humidification.m_w
        return preheated_air, m_dot_steam

    def _determine_heating_coil_load(self) -> Quantity:
        preheat_coil = AirConditioningProcess(
            T_ai=self.mixed_air.Tdb,
            T_ao=self.heated_air.Tdb,
            m_da=self.m_dot_supply
        )
        return preheat_coil.Q_sen

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
