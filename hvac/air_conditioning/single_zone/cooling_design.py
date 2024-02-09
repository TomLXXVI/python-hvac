from dataclasses import dataclass, field
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid
from hvac.air_conditioning import (
    AirConditioningProcess,
    AirStream,
    AdiabaticMixing,
    Fan
)
from hvac.charts import PsychrometricChart, StatePoint


Q_ = Quantity
Air = Fluid('Air')
standard_air = Air(T=Q_(20, 'degC'), P=Q_(101_325, 'Pa'))


@dataclass
class Output:
    m_dot_supply: Quantity
    m_dot_vent: Quantity
    m_dot_recir: Quantity
    outdoor_air: HumidAir
    mixed_air: HumidAir
    cooled_air: HumidAir
    supply_air: HumidAir
    return_air: HumidAir
    zone_air_sp: HumidAir
    Q_dot_sen_zone: Quantity
    Q_dot_lat_zone: Quantity
    SHR_zone: Quantity
    Q_dot_cc: Quantity
    SHR_cc: Quantity
    Q_dot_hc: Quantity | None = None
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
        output = (
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
            "outdoor air = "
            f"{self.outdoor_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.outdoor_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.outdoor_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "mixed air = "
            f"{self.mixed_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.mixed_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.mixed_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "cooled air = "
            f"{self.cooled_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.cooled_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.cooled_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "supply air = "
            f"{self.supply_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.supply_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.supply_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "return air = "
            f"{self.return_air.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.return_air.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.return_air.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "zone air (setpoint) = "
            f"{self.zone_air_sp.Tdb.to(self.units['T'][0]):~P.{self.units['T'][1]}f} DB, "
            f"{self.zone_air_sp.W.to(self.units['W'][0]):~P.{self.units['W'][1]}f} AH "
            f"({self.zone_air_sp.RH.to(self.units['RH'][0]):~P.{self.units['RH'][1]}f} RH)\n"
            "zone cooling load = "
            f"{self.Q_dot_sen_zone.to(self.units['Q_dot'][0]):~P.{self.units['Q_dot'][1]}f} (S), "
            f"{self.Q_dot_lat_zone.to(self.units['Q_dot'][0]):~P.{self.units['Q_dot'][1]}f} (L), "
            f"{self.SHR_zone.to(self.units['SHR'][0]):~P.{self.units['SHR'][1]}f}\n"
            "cooling coil load = "
            f"{self.Q_dot_cc.to(self.units['Q_dot'][0]):~P.{self.units['Q_dot'][1]}f}, "
            f"{self.SHR_cc.to(self.units['SHR'][0]):~P.{self.units['SHR'][1]}f}"
        )
        if self.Q_dot_hc is not None:
            output += '\n'
            output += (
                "heating coil load = "
                f"{self.Q_dot_hc.to(self.units['Q_dot'][0]):~P.{self.units['Q_dot'][1]}f}"
            )
        return output


class AircoSystem:
    """Class for determining the cooling coil load of a single-zone CAV or
    VAV system at summer peak design conditions needed to size the cooling
    coil.

    Attributes
    ----------
    outdoor_air: HumidAir
        If heat recovery is taken into account, the state of the outdoor air
        downstream of the heat recovery device and upstream of the mixing
        plenum.
    supply_air: HumidAir
        Required state of the supply air to compensate for the zone cooling
        load.
    m_dot_supply: PlainQuantity
        Required mass flow rate of supply air to the zone.
    mixed_air: HumidAir
        State of air downstream of the mixing plenum, at the cooling coil inlet.
    m_dot_recir: PlainQuantity
        Mass flow rate of return air at design conditions.
    cooled_air: HumidAir
        Required state of air at the cooling coil outlet.
    Q_dot_cc: PlainQuantity
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
        Q_dot_zone: Quantity,
        SHR_zone: Quantity,
        V_dot_vent_ntp: Quantity,
        T_supply: Quantity,
        dP_fan: Quantity | None = None,
        eta_fan: Quantity | None = None,
        eta_motor: Quantity | None = None,
        eps_hr_h: Quantity | None = None,
        eps_hr_W: Quantity | None = None,
        units: dict[str, tuple[str, int]] | None = None
    ) -> None:
        """Creates an `AircoSystem` instance.

        Parameters
        ----------
        zone_air:
            Desired state of zone air. Should be the same state used
            for the cooling load calculation of the building.
        outdoor_air:
            State of outdoor air on which the design calculations will be
            based. Should be the same state used for the cooling load
            calculation of the building.
        Q_dot_zone:
            Total cooling load of the zone according to the cooling load
            calculation of the building.
        SHR_zone:
            Sensible heat ratio of the zone according to the cooling load
            calculation of the building.
        V_dot_vent_ntp:
            Volume flow rate of outdoor air referred to the NTP air state that is
            needed to fulfill the minimum ventilation requirement of the zone.
        T_supply:
            Design value selected for the supply air dry-bulb temperature.
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
        units: optional
            Units to be used for displaying the results. Keys are:
            'm_dot' for mass flow rate, 'V_dot' for volume flow rate, 'Q_dot'
            for heat transfer rate, 'T' for temperature, 'W' for absolute
            humidity, 'RH' for relative humidity, and 'SHR' for the sensible heat
            ratio. Values are 2-tuples: the first element is the unit
            (a string), the second element is the number of decimals to be
            displayed (an int).
        """
        self.zone_air = self.return_air = zone_air
        self.outdoor_air = outdoor_air
        self.Q_dot_zone = Q_dot_zone
        self.SHR_zone = SHR_zone
        self.V_dot_vent_ntp = V_dot_vent_ntp
        self.T_supply = T_supply
        self.dP_fan = dP_fan
        self.eta_fan = eta_fan
        self.eta_motor = eta_motor
        self.eps_hr_h = eps_hr_h
        self.eps_hr_W = eps_hr_W
        self.units = units or {}

        self.Q_dot_zone_sen = self.SHR_zone * self.Q_dot_zone
        self.Q_dot_zone_lat = self.Q_dot_zone - self.Q_dot_zone_sen

        self.m_dot_vent: Quantity | None = None
        self.supply_air: HumidAir | None = None
        self.m_dot_supply: Quantity | None = None
        self.mixed_air: Quantity | None = None
        self.m_dot_recir: Quantity | None = None
        self.cooled_air: HumidAir | None = None
        self.Q_dot_cc: Quantity | None = None
        self.SHR_cc: Quantity | None = None
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
        cooled_air:
            The required state of air at the cooling coil outlet.
        supply_air:
            The required state of the supply air to compensate for the zone
            cooling loads.
        return_air:
            The state of air returned to the air handler. This is also the state
            of the zone air.
        Q_dot_cc:
            The required cooling coil capacity at design conditions.
        SHR_cc:
            The sensible heat ratio of the cooling coil.
        """
        self.m_dot_vent = self.V_dot_vent_ntp * standard_air.rho
        self.outdoor_air = self._determine_outdoor_air(self.eps_hr_h, self.eps_hr_W)
        self.supply_air, self.m_dot_supply = self._determine_supply_air()
        self.mixed_air, self.m_dot_recir = self._determine_mixed_air()
        self.cooled_air = self._determine_cooled_air(self.eta_fan, self.eta_motor, self.dP_fan)
        self.Q_dot_cc, self.SHR_cc = self.determine_cooling_coil_load()
        self.V_dot_supply_ntp = self.m_dot_supply / standard_air.rho
        self.V_dot_vent = self.m_dot_vent / self.outdoor_air.rho
        return Output(
            self.m_dot_supply, self.m_dot_vent, self.m_dot_recir,
            self.outdoor_air, self.mixed_air, self.cooled_air,
            self.supply_air, self.return_air, self.zone_air,
            self.Q_dot_zone_sen, self.Q_dot_zone_lat, self.SHR_zone,
            self.Q_dot_cc, self.SHR_cc, None, self.units
        )

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
            Q=self.Q_dot_zone,
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
        m_dot_recir = self.m_dot_supply - self.m_dot_vent
        mixing_chamber = AdiabaticMixing(
            in1=AirStream(self.zone_air, m_dot_recir),
            in2=AirStream(self.outdoor_air, self.m_dot_vent),
            out=AirStream(m_da=self.m_dot_supply)
        )
        return mixing_chamber.stream_out.state, m_dot_recir

    def _determine_cooled_air(
        self,
        eta_fan: Quantity,
        eta_motor: Quantity,
        dP_fan: Quantity
    ) -> HumidAir:
        # determine the state of air at the cooling coil outlet: if supply
        # fan air heating is taken into account and assuming a draw-through
        # arrangement, the state of the cooled air is determined to
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
            m_da=self.m_dot_supply,
            h_w=Q_(0.0, 'J / kg')
        )
        return cooling_coil.Q, cooling_coil.SHR

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
            name='fan heating',
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
