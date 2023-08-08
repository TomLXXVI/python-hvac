from hvac import Quantity
from hvac.fluids import HumidAir, Fluid
from hvac.air_conditioning import (
    AirConditioningProcess,
    AirStream,
    AdiabaticMixing,
    Fan
)

Q_ = Quantity
Air = Fluid('Air')
air_ntp = Air(T=Q_(20, 'degC'), P=Q_(101_325, 'Pa'))


class AircoSystem:
    """Class for determining the cooling coil load of a single-zone CAV or
    VAV system at design conditions needed for sizing the cooling coil.

    The design calculations are executed on instantiation of the class.
    After instantiation, the following attributes are available:

    Attributes
    ----------
    outdoor_air: HumidAir
        If heat recovery is taken into account, the state of the outdoor air
        downstream of the heat recovery device and upstream of the mixing
        plenum.
    supply_air: HumidAir
        Required state of the supply air to compensate for the zone cooling
        load.
    m_supply: PlainQuantity
        Required mass flow rate of supply air to the zone.
    mixed_air: HumidAir
        State of air downstream of the mixing plenum, at the cooling coil inlet.
    m_recir: PlainQuantity
        Mass flow rate of return air at design conditions.
    cooled_air: HumidAir
        Required state of air at the cooling coil outlet.
    Q_cc: PlainQuantity
        Required cooling coil capacity at design conditions.
    SHR_cc: PlainQuantity
        Sensible heat ratio of the cooling coil.

    Notes
    -----
    It is implicitly assumed that the supply fan is downstream from the cooling
    coil, i.e., a draw-through arrangement of the supply fan is assumed.
    """

    def __init__(
        self,
        zone_air: HumidAir,
        outdoor_air: HumidAir,
        Q_zone: Quantity,
        SHR_zone: Quantity,
        V_vent: Quantity,
        T_supply: Quantity,
        dP_fan: Quantity | None = None,
        eta_fan: Quantity | None = None,
        eta_motor: Quantity | None = None,
        eps_hr_h: Quantity | None = None,
        eps_hr_W: Quantity | None = None
    ) -> None:
        """Creates an `AircoSystem` instance.

        Parameters
        ----------
        zone_air:
            Desired state of zone air. Should be the same state that was used
            for the cooling load calculation of the building.
        outdoor_air:
            State of outdoor air on which the design calculations will be
            based. Should be the same state that was used for the cooling load
            calculation of the building.
        Q_zone:
            Total cooling load of the zone according to the cooling load
            calculation of the building.
        SHR_zone:
            Sensible heat ratio of the zone according to the cooling load
            calculation of the building.
        V_vent:
            Volume flow rate of outdoor air needed to ventilate the zone based
            on the minimum ventilation requirement of the zone.
        T_supply:
            Dry-bulb temperature of supply air to the zone selected for the
            design.
        dP_fan: optional
            Supply fan pressure. The fan pressure can be determined from a
            pressure drop calculation when sizing the duct system.
        eta_fan: optional
            Supply fan efficiency. The fan efficiency can be determined after
            the supply fan has been selected, i.e., after sizing the duct system
            first.
        eta_motor: optional
            Supply fan motor efficiency. The efficiency of the fan motor needs
            to be taken into account if the fan motor is also located in the
            air stream.
        eps_hr_h: optional
            Enthalpy effectiveness of the heat recovery system.
        eps_hr_W: optional
            Humidity ratio effectiveness of the heat recovery system.
        """
        self.zone_air = zone_air
        self.outdoor_air = outdoor_air
        self.Q_zone = Q_zone
        self.SHR_zone = SHR_zone
        self.V_vent = V_vent
        self.T_supply = T_supply

        self.outdoor_air = self._determine_outdoor_air(eps_hr_h, eps_hr_W)
        self.supply_air, self.m_supply = self._determine_supply_air()
        self.mixed_air, self.m_recir = self._determine_mixed_air()
        self.cooled_air = self._determine_cooled_air(eta_fan, eta_motor, dP_fan)
        self.Q_cc, self.SHR_cc = self.determine_cooling_coil_load()

    def _determine_outdoor_air(
        self,
        eps_hr_h: Quantity,
        eps_hr_W: Quantity
    ) -> HumidAir:
        # get the state of outdoor air after heat recovery if it is present
        if eps_hr_h is not None and eps_hr_W is not None:
            # heat recovery present
            dh_max = self.outdoor_air.h - self.zone_air.h
            dh = eps_hr_h * dh_max
            h_ao = self.outdoor_air.h - dh
            dW_max = self.outdoor_air.W - self.zone_air.W
            dW = eps_hr_W * dW_max
            W_ao = self.outdoor_air.W - dW
            outdoor_air = HumidAir(h=h_ao, W=W_ao)
            return outdoor_air
        else:
            # no heat recovery taken into account
            return self.outdoor_air

    def _determine_supply_air(self) -> tuple[HumidAir, Quantity]:
        # determine the required state of the supply air and the required
        # mass flow rate of supply air based on the selection of the dry-bulb
        # supply temperature
        zone = AirConditioningProcess(
            T_ai=self.T_supply,
            air_out=self.zone_air,
            Q=self.Q_zone,
            SHR=self.SHR_zone
        )
        return zone.air_in, zone.m_da

    def _determine_mixed_air(self) -> tuple[HumidAir, Quantity]:
        # determine the state of air at the inlet of the cooling coil after
        # ventilation outdoor air and return air from the zone have been mixed
        # in the mixing plenum of the AHU. Also determine the resulting mass
        # flow rate of return air, that follows from the difference between the
        # required mass flow rate of supply air and the required ventilation
        # air mass flow rate.
        m_vent = self.outdoor_air.rho * self.V_vent
        m_recir = self.m_supply - m_vent
        mixing_chamber = AdiabaticMixing(
            in1=AirStream(self.zone_air, m_recir),
            in2=AirStream(self.outdoor_air, m_vent),
            out=AirStream(m_da=self.m_supply)
        )
        return mixing_chamber.stream_out.state, m_recir

    def _determine_cooled_air(
        self,
        eta_fan: Quantity,
        eta_motor: Quantity,
        dP_fan: Quantity
    ) -> HumidAir:
        # determine the state of air at the cooling coil outlet: if supply
        # fan air heating is taken into account and assuming a draw-through
        # arrangement, the state of the cooled air is determined in order to
        # compensate for air heating due to supply fan inefficiencies.
        if all((eta_fan, eta_motor, dP_fan)):
            supply_fan = Fan(
                air_out=self.supply_air,
                eta_fan=eta_fan,
                eta_motor=eta_motor,
                dP_fan=dP_fan
            )
            return supply_fan.air_in
        else:
            return self.supply_air

    def determine_cooling_coil_load(self) -> tuple[Quantity, Quantity]:
        # determine the required cooling coil capacity and cooling coil sensible
        # heat ratio
        cooling_coil = AirConditioningProcess(
            air_in=self.mixed_air,
            air_out=self.cooled_air,
            m_da=self.m_supply,
            h_w=Q_(0.0, 'J / kg')
        )
        return cooling_coil.Q, cooling_coil.SHR

    @property
    def V_supply_ntp(self) -> Quantity:
        """Get the volume flow rate of supply air referred to NTP."""
        return self.m_supply / air_ntp.rho
