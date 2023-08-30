from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning import (
    AirConditioningProcess,
    AirStream,
    AdiabaticMixing
)

Q_ = Quantity


class CAVSystem:
    """Class for analyzing a CAV system at part-load conditions.

    In a CAV system:
    1.  The mass flow rate of supply air remains constant (constant-speed fan).
    2.  The air condition at the cooling coil outlet is controlled to
        keep it constant. The cooling coil's main purpose is to dehumidify the
        air.
    3.  As the state of the air leaving the cooling coil is kept constant, a
        reheat coil is necessary to control the supply air temperature in order
        to maintain the desired zone air temperature at zone loads which are
        smaller than the design load.

    Attributes
    ----------
    supply_air:
        State of supply air at part-load.
    zone_air:
        State of zone air at part-load.
    mixed_air:
        State of mixed air at part-load.
    Q_cc:
        Cooling coil load at part-load.
    SHR_cc:
        Cooling coil sensible heat ratio at part-load.
    Q_rh:
        Reheat-coil load at part-load.
    """
    def __init__(
        self,
        T_zone: Quantity,
        outdoor_air: HumidAir,
        Q_zone: Quantity,
        SHR_zone: Quantity,
        V_vent: Quantity,
        m_supply: Quantity,
        cooled_air: HumidAir
    ) -> None:
        """
        Parameters
        ----------
        T_zone:
            Zone air temperature (fixed by zone thermostat)
        outdoor_air:
            State of outdoor air.
        Q_zone:
            Total cooling load of the zone.
        SHR_zone:
            Sensible heat ratio of the zone load.
        V_vent:
            Volume flow rate of ventilation air at outdoor air state.
        m_supply:
            Mass flow rate of supply air (fixed)
        cooled_air:
            State of air at the cooling coil outlet (fixed by controller).
        """
        self.T_zone = T_zone
        self.outdoor_air = outdoor_air
        self.Q_zone = Q_zone
        self.SHR_zone = SHR_zone
        self.V_vent = V_vent
        self.m_supply = m_supply
        self.cooled_air = cooled_air

        self.supply_air, self.zone_air = self._determine_supply_air()
        self.mixed_air = self._determine_mixed_air()
        self.Q_cc, self.SHR_cc = self._determine_cooling_coil_load()
        self.Q_rh = self._determine_reheat_coil_load()

    def _determine_supply_air(self) -> tuple[HumidAir, HumidAir]:
        zone = AirConditioningProcess(
            W_ai=self.cooled_air.W,
            T_ao=self.T_zone,
            m_da=self.m_supply,
            Q=self.Q_zone,
            SHR=self.SHR_zone,
            h_w=Q_(0, 'J / kg')
        )
        return zone.air_in, zone.air_out

    def _determine_mixed_air(self) -> HumidAir:
        m_vent = self.outdoor_air.rho * self.V_vent
        m_recir = self.m_supply - m_vent
        mixing_chamber = AdiabaticMixing(
            in1=AirStream(self.zone_air, m_recir),
            in2=AirStream(self.outdoor_air, m_vent),
            out=AirStream(m_da=self.m_supply)
        )
        return mixing_chamber.stream_out.state

    def _determine_cooling_coil_load(self) -> tuple[Quantity, Quantity]:
        cooling_coil = AirConditioningProcess(
            air_in=self.mixed_air,
            air_out=self.cooled_air,
            m_da=self.m_supply
        )
        return cooling_coil.Q, cooling_coil.SHR

    def _determine_reheat_coil_load(self) -> Quantity:
        reheat_coil = AirConditioningProcess(
            air_in=self.cooled_air,
            air_out=self.supply_air,
            m_da=self.m_supply
        )
        return reheat_coil.Q_sen


class VAVSystem:
    """Class for analyzing a VAV system at part-load conditions.

    In a VAV system:
    1.  The air condition at the cooling coil outlet is kept constant to
        dehumidify the air. If sensible heating due to inefficiencies in the
        supply fan is ignored, the supply air condition corresponds with the
        air condition at the cooling coil outlet.
    2.  The space thermostat modulates the mass flow rate of supply air to
        maintain the desired zone air temperature (variable-speed fan).

    Attributes
    ----------
    m_supply:
        Required mass flow rate of supply air to offset the sensible load of
        the zone.
    zone_air:
        State of zone air at part-load.
    mixed_air:
        State of mixed air at part-load.
    Q_cc:
        Cooling coil load at part-load.
    SHR_cc:
        Cooling coil sensible heat ratio at part-load.
    """

    def __init__(
        self,
        T_zone: Quantity,
        outdoor_air: HumidAir,
        Q_zone: Quantity,
        SHR_zone: Quantity,
        V_vent: Quantity,
        cooled_air: HumidAir
    ) -> None:
        """
        Parameters
        ----------
        T_zone:
            Zone air temperature (fixed by zone thermostat)
        outdoor_air:
            State of outdoor air.
        Q_zone:
            Total cooling load of the zone.
        SHR_zone:
            Sensible heat ratio of the zone load.
        V_vent:
            Volume flow rate of ventilation air at outdoor air state.
        cooled_air:
            State of air at the cooling coil outlet (fixed by controller).
        """
        self.T_zone = T_zone
        self.outdoor_air = outdoor_air
        self.Q_zone = Q_zone
        self.SHR_zone = SHR_zone
        self.V_vent = V_vent
        self.cooled_air = cooled_air

        self.supply_air = cooled_air
        self.m_supply, self.zone_air = self._determine_supply_air_flow_rate()
        self.mixed_air = self._determine_mixed_air()
        self.Q_cc, self.SHR_cc = self._determine_cooling_coil_load()

    def _determine_supply_air_flow_rate(self) -> tuple[Quantity, HumidAir]:
        zone = AirConditioningProcess(
            air_in=self.supply_air,
            T_ao=self.T_zone,
            Q=self.Q_zone,
            SHR=self.SHR_zone
        )
        return zone.m_da, zone.air_out

    def _determine_mixed_air(self) -> HumidAir:
        m_vent = self.outdoor_air.rho * self.V_vent
        m_recir = self.m_supply - m_vent
        mixing_chamber = AdiabaticMixing(
            in1=AirStream(self.zone_air, m_recir),
            in2=AirStream(self.outdoor_air, m_vent),
            out=AirStream(m_da=self.m_supply)
        )
        return mixing_chamber.stream_out.state

    def _determine_cooling_coil_load(self) -> tuple[Quantity, Quantity]:
        cooling_coil = AirConditioningProcess(
            air_in=self.mixed_air,
            air_out=self.cooled_air,
            m_da=self.m_supply
        )
        return cooling_coil.Q, cooling_coil.SHR
