from hvac import Quantity
from hvac.fluids import HumidAir, Fluid
from hvac.air_conditioning import (
    AirConditioningProcess,
    AirStream,
    AdiabaticMixing
)


Q_ = Quantity
Water = Fluid('Water')
Air = Fluid('Air')
air_ntp = Air(T=Q_(20, 'degC'), P=Q_(101_325, 'Pa'))


class AircoSystem:
    """Class to design an AC system (CAV or VAV) for the heating and the
    humidification of supply air to a single-zone building in order to meet the
    sensible and latent zone load. Its purpose is for sizing and selecting the
    main components of the AC system.

    After instantiation of the class, the following calculated attributes are
    available:

    Attributes
    ----------
    m_dot_supply:
        The required mass flow rate of supply air to the zone.
    supply_air:
        The state of supply air under peak design conditions.
    m_dot_steam:
        The required saturated steam mass flow rate for humidification of the
        supply air to the zone.
    Q_dot_preheat:
        The required heating capacity of the preheat-coil to heat the supply air
        to the zone. For winter peak design, it is assumed that air-heating only
        takes place in the preheat-coil. An additional heating coil is not being
        considered under peak design conditions.
    preheated_air:
        The state of air at the preheat-coil outlet and entering the humidifier
        under peak design conditions.
    mixed_air:
        The state of air after mixing outdoor air with return air from the zone
        and entering the preheat-coil under peak design conditions.
    m_dot_recir:
        The mass flow rate of recirculated air under peak design conditions.
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
        """
        self.zone_air = self.return_air = zone_air
        self.outdoor_air = outdoor_air
        self.Q_dot_zone = Q_dot_zone
        self.SHR_zone = SHR_zone
        self.V_dot_vent_ntp = V_dot_vent_ntp
        self.m_dot_supply_cooling = m_dot_supply_cooling
        self.T_supply_max = T_supply_max
        self.T_steam = T_steam

        self.m_dot_vent = self.V_dot_vent_ntp * air_ntp.rho
        self.supply_air, self.m_dot_supply = self._determine_supply_air()
        self.m_dot_recir = self.m_dot_supply - self.m_dot_vent
        self.mixed_air = self._determine_mixed_air()
        self.preheated_air, self.m_dot_steam = self._determine_preheated_air()
        self.Q_dot_preheat = self._determine_preheat_coil_load()

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

    def _determine_preheated_air(self) -> tuple[HumidAir, Quantity]:
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

    def _determine_preheat_coil_load(self) -> Quantity:
        preheat_coil = AirConditioningProcess(
            T_ai=self.mixed_air.Tdb,
            T_ao=self.preheated_air.Tdb,
            m_da=self.m_dot_supply
        )
        return preheat_coil.Q_sen

    @property
    def V_dot_supply_ntp(self) -> Quantity:
        """Get the volume flow rate of supply air referred to NTP."""
        return self.m_dot_supply / air_ntp.rho

    @property
    def V_dot_vent(self) -> Quantity:
        """Get the volume flow rate of outdoor ventilation air referred to
        the actual state of outdoor air.
        """
        return self.m_dot_vent / self.outdoor_air.rho
