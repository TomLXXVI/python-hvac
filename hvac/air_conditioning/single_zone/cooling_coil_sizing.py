import warnings
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid
from hvac.air_conditioning import (
    AirStream,
    AdiabaticMixing,
    AirConditioningProcess
)
Q_ = Quantity
Air = Fluid('Air')
air_ntp = Air(T=Q_(20, 'degC'), P=Q_(101_325, 'Pa'))


class DxAirCoolingCoilSizer:
    """Class for sizing a DX air-cooling coil used for cooling and dehumidifying
     air in a single-zone air-cooling system.

    The DX air cooler is selected from a catalog/datasheet having a known
    cooling capacity at the specified operating conditions.
    The design procedure determines the required face area of this air cooler,
    while keeping the face velocity and the refrigerant temperature (evaporation
    temperature at the specified values given in the datasheet. As such, the
    heat and the mass transfer effectiveness of the air cooler can be assumed
    to remain the same as the values that are valid under the given set of
    operating conditions in the datasheet.
    """
    class WetDxAirCooler:

        def __init__(
            self,
            air_in_ref: HumidAir,
            air_out_ref: HumidAir,
            T_rfg: Quantity,
            v_fa_ntp: Quantity,
            dP_ntp: Quantity
        ) -> None:
            """
            Class that represents a DX air-cooler selected from a catalog with
            a known set of operating conditions. A fully wetted air-side heat
            transfer surface is assumed.

            Parameters
            ----------
            air_in_ref:
                The state of the air at the air cooler inlet taken from the
                datasheet.
            air_out_ref:
                The state of the air at the air cooler outlet taken from the
                datasheet.
            T_rfg:
                The refrigerant's evaporation temperature taken from the
                datasheet.
            v_fa_ntp:
                The face velocity of air referred to NTP conditions taken from
                the datasheet.
            dP_ntp:
                Pressure loss over the air cooler which corresponds with the
                NTP face velocity, taken from the datasheet.
            """
            self.air_out_sat = HumidAir(Tdb=T_rfg, RH=Q_(100, 'pct'))
            self.eps_h = self._enthalpy_effectiveness(air_in_ref, air_out_ref)
            self.eps_W = self._humidity_effectiveness(air_in_ref, air_out_ref)
            self.v_fa_ntp = v_fa_ntp
            self.dP_ntp = dP_ntp
            self.air_in: HumidAir | None = None

        def _enthalpy_effectiveness(
            self,
            air_in: HumidAir,
            air_out: HumidAir
        ) -> float:
            # Determine the total energy transfer effectiveness from the known
            # inlet and outlet air conditions.
            dh = (air_out.h - air_in.h).to('J / kg')
            dh_max = (self.air_out_sat.h - air_in.h).to('J / kg')
            eps_h = dh / dh_max
            return eps_h.m

        def _humidity_effectiveness(
            self,
            air_in: HumidAir,
            air_out: HumidAir
        ) -> float:
            # Determine the mass transfer effectiveness from the known
            # inlet and outlet air conditions.
            dW = (air_out.W - air_in.W).to('kg / kg')
            dW_max = (self.air_out_sat.W - air_in.W).to('kg / kg')
            eps_W = dW / dW_max
            return eps_W.m

        @property
        def air_out(self) -> HumidAir:
            """Get the actual air state at the cooler outlet based on the
            assumption that heat and mass transfer effectiveness are constants.
            """
            if isinstance(self.air_in, HumidAir):
                h_ao = (
                    self.air_in.h + self.eps_h
                    * (self.air_out_sat.h - self.air_in.h)
                )
                W_ao = (
                    self.air_in.W + self.eps_W
                    * (self.air_out_sat.W - self.air_in.W)
                )
                air_out = HumidAir(h=h_ao, W=W_ao)
                return air_out

        @property
        def q(self) -> Quantity:
            """Get the specific cooling capacity of the air cooler at the actual
            condition of the inlet air, i.e., the cooling capacity per unit face
            area.
            """
            if isinstance(self.air_in, HumidAir):
                q = (
                    self.eps_h * air_ntp.rho * self.v_fa_ntp
                    * (self.air_out_sat.h - self.air_in.h)
                )
                return q

    def __init__(
        self,
        zone_air: HumidAir,
        outdoor_air: HumidAir,
        Q_zone: Quantity,
        SHR_zone: Quantity,
        f_vent: Quantity,
        air_cooler: WetDxAirCooler,
        i_max: int = 20,
        W_tol: Quantity = Q_(5.e-4, 'g / kg')
    ) -> None:
        """Creates an instance of `DxAirCoolingCoilSizer` initialized with the
        design information for the single-zone air-cooling system. All the
        calculations are done on instantiation of the object. Method `info()`
        returns a string with the calculation results.

        Parameters
        ----------
        zone_air:
            Desired zone air state at design load conditions of the zone.
        outdoor_air:
            State of outdoor air at design load conditions of the zone.
        Q_zone:
            Total cooling load of the zone at design load conditions.
        SHR_zone:
            Sensible heat ratio of the zone's design cooling load.
        f_vent:
            Ventilation air flow rate expressed as a fraction of the supply
            air mass flow rate to the zone.
        air_cooler:
            Instance of class `WetDxAirCooler`. The air cooler selected from
            the catalog that will be used to cool and dehumidify the supply
            air to the zone.
        i_max:
            Maximum number of iterations to find the final zone air state.
        W_tol:
            Allowable tolerance or deviation for the final zone air humidity
            ratio. If after the maximum number of iterations, the humidity ratio
            falls outside the tolerance margin, a warning will be given.

        Attributes
        ----------
        self.A_fa:
            The required face area of the air-cooling coil.
        self.Q_cc:
            The resulting cooling coil load.
        self.W_fan_min:
            Aeraulic power to be delivered to the supply air needed to establish
            the required mass flow rate of supply air to the zone.
        self.m_supply:
            The required mass flow rate of supply air to the zone.
        self.V_supply_ntp:
            The required volume flow rate of supply air referred to NTP.
        self.m_vent:
            The resulting mass flow rate of outdoor ventilation air.
        self.V_vent_ntp:
            The resulting volume flow rate of outdoor ventilation air referred
            to NTP.
        self.m_recir:
            The resulting mass flow rate of recirculation air.
        self.V_recir_ntp:
            The resulting volume flow rate of recirculation air referred to NTP.
        self.supply_air:
            The required state of the supply air to the zone.
        self.zone_air:
            The resulting state of the air in the zone.
        self.mixed_air:
            The resulting state of air entering the cooling coil after adiabatic
            mixing of recirculated zone air and outdoor air.
        """
        self.zone_air = zone_air
        self.outdoor_air = outdoor_air
        self.Q_zone = Q_zone
        self.SHR_zone = SHR_zone.to('frac')
        self.f_vent = f_vent.to('frac')
        self.air_cooler = air_cooler
        self.W_tol = W_tol.to('kg / kg')

        # Determine the required mass flow rate of supply air to the zone.
        # Iteration is used to find the final state of the zone air. The
        # iteration starts with the desired zone air state set by the user as an
        # initial guess and stops when the zone air humidity ratio doesn't change
        # anymore between two successive loops, i.e. when the difference between
        # the present and previous value is smaller than the allowable tolerance,
        # or when the maximum number of iterations has been reached.
        for i in range(i_max):
            self.mixed_air = self._determine_mixed_air()
            self.supply_air = self._determine_supply_air()
            self.m_supply, self.V_supply_ntp = self._determine_supply_air_flow_rate()
            zone_air_new = self._determine_zone_air()
            if self._check_zone_humidity_ratio(zone_air_new):
                break
            self.zone_air = zone_air_new
        else:
            warnings.warn(
                'the maximum number of iterations to find the final zone air '
                'state is exceeded', category=RuntimeWarning
            )

        self.A_fa = self._determine_face_area()
        self.Q_cc = self._determine_cooling_capacity()
        self.W_fan_min = self._determine_aeraulic_fan_power()
        self.m_vent, self.V_vent_ntp = self._determine_ventilation_flow_rate()
        self.m_recir, self.V_recir_ntp = self._determine_recirculation_flow_rate()

    def _determine_mixed_air(self) -> HumidAir:
        # To determine the state of the mixed air, we impose that the
        # ventilation air mass flow rate should be a given fraction of the total
        # supply air mass flow rate. The fraction of recirculation air mass
        # flow rate is then also determined. Based on these fractions, the state
        # of mixed air can be determined, as the zone air state and outdoor air
        # state are known.
        mixing_chamber = AdiabaticMixing(
            in1=AirStream(self.zone_air, m_da_fractional=(1 - self.f_vent)),
            in2=AirStream(self.outdoor_air, m_da_fractional=self.f_vent)
        )
        return mixing_chamber.stream_out.state

    def _determine_supply_air(self) -> HumidAir:
        # The mixed air goes to the inlet of the air cooler. The state of the
        # supply air at the outlet of the air cooler is determined by the
        # air cooler property `air_out` (see class `WetDxAirCooler`).
        self.air_cooler.air_in = self.mixed_air
        return self.air_cooler.air_out

    def _determine_supply_air_flow_rate(self) -> tuple[Quantity, Quantity]:
        # Determine the mass flow rate of supply air from a sensible heat
        # balance of the zone. This is the mass flow rate of supply air that is
        # needed to maintain the desired zone air temperature and compensates
        # for the sensible cooling load of the zone.
        zone = AirConditioningProcess(
            air_in=self.supply_air,
            T_ao=self.zone_air.Tdb,
            Q_sen=self.SHR_zone * self.Q_zone
        )
        m_supply = zone.m_da
        V_supply_ntp = m_supply / air_ntp.rho
        return m_supply, V_supply_ntp

    def _determine_zone_air(self) -> HumidAir:
        # Determine the state of the zone air.
        zone = AirConditioningProcess(
            air_in=self.supply_air,
            m_da=self.m_supply,
            Q=self.Q_zone,
            SHR=self.SHR_zone
        )
        return zone.air_out

    def _check_zone_humidity_ratio(self, zone_air: HumidAir) -> bool:
        # Check if the difference between the current and previous humidity
        # ratio of the zone air is within the allowable tolerance margin.
        W_new = zone_air.W.to('kg / kg')
        W = self.zone_air.W.to('kg / kg')
        if abs(W_new - W) < self.W_tol:
            return True
        return False

    def _determine_face_area(self) -> Quantity:
        # Determine the face area of the air cooler such that for the fixed
        # face velocity taken from the datasheet the required mass flow rate of
        # supply air is delivered to the zone.
        A_fa = self.V_supply_ntp / self.air_cooler.v_fa_ntp
        return A_fa

    def _determine_cooling_capacity(self) -> Quantity:
        Q_cc = self.A_fa * self.air_cooler.q
        return Q_cc

    def _determine_aeraulic_fan_power(self) -> Quantity:
        W_fan_min = self.V_supply_ntp * self.air_cooler.dP_ntp
        return W_fan_min

    def _determine_ventilation_flow_rate(self) -> tuple[Quantity, Quantity]:
        m_vent = self.f_vent * self.m_supply
        V_vent_ntp = m_vent / air_ntp.rho
        return m_vent, V_vent_ntp

    def _determine_recirculation_flow_rate(self) -> tuple[Quantity, Quantity]:
        m_recir = (1 - self.f_vent) * self.m_supply
        V_recir_ntp = m_recir / air_ntp.rho
        return m_recir, V_recir_ntp

    def info(self) -> str:
        """Returns an overview of the calculation results.

        It should be checked that the zone air humidity is within acceptable
        limits and that the resulting ventilation air flow rate fulfills the
        ventilation requirement of the zone.
        """
        info = [
            (
                "- mixed air: "
                f"{self.mixed_air.Tdb.to('degC'):~P.1f} DB "
                f"{self.mixed_air.RH.to('pct'):~P.1f} RH"
            ),
            (
                "- supply air: "
                f"{self.supply_air.Tdb.to('degC'):~P.1f} DB "
                f"{self.supply_air.RH.to('pct'):~P.0f} RH"
            ),
            (
                "- zone air: "
                f"{self.zone_air.Tdb.to('degC'):~P.1f} DB "
                f"{self.zone_air.RH.to('pct'):~P.0f} RH"
            ),
            (
                "- mass flow rate supply air: "
                f"{self.m_supply.to('kg / s'):~P.1f}"
            ),
            (
                "- volume flow rate supply air @ NTP: "
                f"{self.V_supply_ntp.to('m ** 3 / hr'):~P.1f}"
            ),
            (
                "- mass flow rate ventilation air: "
                f"{self.m_vent.to('kg / s'):~P.1f}"
            ),
            (
                "- volume flow rate ventilation air @ NTP: "
                f"{self.V_vent_ntp.to('m ** 3 / hr'):~P.1f}"
            ),
            (
                "- mass flow rate recirculation air: "
                f"{self.m_recir.to('kg / s'):~P.1f}"
            ),
            (
                "- volume flow rate recirculation air @ NTP: "
                f"{self.V_recir_ntp.to('m ** 3 / hr'):~P.1f}"
            ),
            (
                "- required face area of air cooler: "
                f"{self.A_fa.to('m ** 2'):~P.3f}"
            ),
            (
                "- cooling capacity of air cooler: "
                f"{self.Q_cc.to('kW'):~P.3f}"
            ),
            (
                "- aeraulic power needed for air flow through cooler: "
                f"{self.W_fan_min.to('W'):~P.0f}"
            )
        ]
        return "\n".join(info)
