import warnings
import numpy as np
from scipy import optimize
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState, CP_HUMID_AIR
import hvac.heat_transfer.heat_exchanger.fin_tube.core as fin_tube_core
from hvac.air_conditioning import AirConditioningProcess


Q_ = Quantity
HexCore = fin_tube_core.plain_fin_tube.PlainFinTubeHeatExchangerCore


class DXAirCoolingCoilWarning(Warning):
    pass


class DXAirCoolingCoil:
    """Simplified model of a controlled DX coil for cooling air in an air
    conditioning system. The DX coil is a plain fin-tube, counter-flow heat
    exchanger.

    Known inputs on DX coil model:
    - mass flow rate of air through DX coil (determined by the sensible cooling
    load of the space)
    - state of entering air (i.e. the state of air after adiabatic mixing)
    - setpoint temperature of leaving air (controlled to remain fixed)
    - state of refrigerant entering the DX coil on the refrigerant side.

    Only boiling refrigerant is being considered (superheating of refrigerant
    vapor is ignored to keep the model simple). It is assumed that refrigerant
    always leaves the DX coil as saturated vapor (as controlled by the expansion
    device on the cooling machine). To keep the temperature of the air leaving
    the cooling coil constant, the mass flow rate of refrigerant will be
    modulated (compressor with variable speed drive). However, for a given
    setpoint of the leaving air temperature, there will be a lower limit on the
    absolute humidity ratio of the leaving air (the straight line on the
    psychrometric chart that connects the entering and leaving air state, must
    intersect or at least touch the saturation line, i.e, the cooling coil will
    always have an Apparatus Dew Point (ADP)).

    The DX air cooling coil model calculates:
    - the heat transfer effectiveness of the DX coil.
    - the heat transfer rate of the DX coil.
    - the required mass flow rate of refrigerant so that refrigerant leaves the
    cooling coil as saturated vapor.
    - the full state of the air leaving the DX coil, i.e., also including its
    humidity level

    Attributes
    ----------
    eps: Quantity
        Heat transfer effectiveness of the DX coil.
    Q: Quantity
        Heat transfer rate of the DX coil.
    m_dot_rfg: Quantity
        Mass flow rate of refrigerant.
    air_out: HumidAir
        Full state of leaving air including its humidity
    BF: Quantity
        Bypass factor of the DX coil.
    CF: Quantity
        Contact factor of the DX coil.
    SHR: Quantity
        Sensible heat ratio of the DX coil.
    ADP: Quantity
        Apparatus Dew Point of the DX coil.
    """
    def __init__(
        self,
        m_dot_air: Quantity,
        air_in: HumidAir,
        T_air_out: Quantity,
        rfg_in: FluidState,
        hex_core: HexCore
    ) -> None:
        """
        Parameters
        ----------
        m_dot_air:
            Mass flow rate of air through heat exchanger.
        air_in:
            State of air entering the heat exchanger.
        T_air_out:
            Temperature of air leaving the heat exchanger.
        rfg_in:
            State of refrigerant entering the heat exchanger.
        hex_core:
            Heat exchanger core.

        """
        self.m_dot_air = m_dot_air.to('kg / s')
        self.air_in = air_in
        self.T_air_out = T_air_out.to('K')
        self.rfg_in = rfg_in
        self.hex_core = hex_core
        self.Rfg = rfg_in.fluid
        self.rfg_out = self.Rfg(P=rfg_in.P, x=Q_(1, 'frac'))  # saturated vapor
        self.Q_max = self.__Q_max__()

        self.m_dot_rfg: Quantity | None = None
        self.eps: Quantity | None = None
        self.Q: Quantity | None = None
        self.air_out: HumidAir | None = None
        self.BF: Quantity | None = None
        self.CF: Quantity | None = None
        self.SHR: Quantity | None = None
        self.ADP: HumidAir | None = None

        self.__calculate_operating_state__()

    @staticmethod
    def __mean_air__(
        air_in: HumidAir,
        air_out: HumidAir,
        rfg_in: FluidState,
        rfg_out: FluidState
    ) -> HumidAir:
        """Calculates the mean state of air along the heat transfer surface."""
        # Saturated air with refrigerant outlet temperature:
        air_sat_rfg_out = HumidAir(
            Tdb=rfg_out.T,
            RH=Q_(100, 'pct')
        )
        # Saturated air with refrigerant inlet temperature:
        air_sat_rfg_in = HumidAir(
            Tdb=rfg_in.T,
            RH=Q_(100, 'pct')
        )
        # Average enthalpy of saturated air along heat transfer surface:
        h_avg_air_sat = (air_sat_rfg_out.h + air_sat_rfg_in.h) / 2
        # Enthalpy potential at air inlet (= refrigerant outlet):
        dh_air_in = air_in.h - air_sat_rfg_out.h
        # Enthalpy potential at air outlet (= refrigerant inlet):
        dh_air_out = air_out.h - air_sat_rfg_in.h
        # Calculate LMED:
        dh_max = max(dh_air_in, dh_air_out)
        dh_min = min(dh_air_in, dh_air_out)
        lmed = (dh_max - dh_min) / np.log(dh_max / dh_min)
        # Mean enthalpy of air along heat transfer surface:
        h_avg_air = h_avg_air_sat + lmed
        # Average temperature of refrigerant along heat transfer surface:
        T_rfg_avg = (rfg_in.T + rfg_out.T) / 2
        # Temperature potential at air inlet (= refrigerant outlet):
        dT_air_in = air_in.Tdb - rfg_out.T
        # Temperature potential at air outlet (= refrigerant inlet):
        dT_air_out = air_out.Tdb - rfg_in.T
        # Calculate LMTD:
        dT_max = max(dT_air_in, dT_air_out)
        dT_min = min(dT_air_in, dT_air_out)
        lmtd = (dT_max - dT_min) / np.log(dT_max / dT_min)
        # Mean temperature of air along the heat transfer surface:
        T_avg_air = T_rfg_avg + lmtd
        # Mean air state along the heat transfer surface:
        air_avg = HumidAir(Tdb=T_avg_air, h=h_avg_air)
        return air_avg

    @staticmethod
    def __mean_refrigerant__(
        rfg_in: FluidState,
        rfg_out: FluidState
    ) -> FluidState:
        """Calculates the mean state of refrigerant along the heat transfer
        surface.
        """
        # Average enthalpy of refrigerant along heat transfer surface:
        h_avg_rfg = (rfg_in.h + rfg_out.h) / 2
        # Average pressure of refrigerant along heat transfer surface:
        P_avg_rfg = (rfg_in.P + rfg_out.P) / 2
        # Mean refrigerant state along the heat transfer surface:
        rfg_avg = rfg_in.fluid(h=h_avg_rfg, P=P_avg_rfg)
        return rfg_avg

    def __Q_max__(self) -> Quantity:
        """Calculates the theoretic maximum heat transfer rate that would
        be possible given the entering state of air, the entering state of
        refrigerant, and the mass flow rate of air through the cooling coil.
        """
        air_out_sat = HumidAir(Tdb=self.rfg_in.T, RH=Q_(100, 'pct'))
        Q_max = self.m_dot_air * (self.air_in.h - air_out_sat.h)
        return Q_max.to('W')

    def __h_ext__(
        self,
        air_mean: HumidAir,
        rfg_mean: FluidState,
        m_dot_rfg: Quantity
    ) -> Quantity:
        """Calculates the air-side convection heat transfer coefficient based on
        the mean states of air and refrigerant, and the mass flow rates of air
        and refrigerant.
        """
        self.hex_core.m_dot_ext = self.m_dot_air
        self.hex_core.ext.fluid_mean = air_mean
        self.hex_core.m_dot_int = m_dot_rfg
        self.hex_core.int.fluid_mean = rfg_mean
        return self.hex_core.ext.h * self.hex_core.ext.eta

    def __UA_ext_wet__(self, h_ext: Quantity) -> Quantity:
        """Calculates the enthalpy-based overall conductance of the cooling coil
        referred to the external, air-side heat transfer surface.
        """
        h_ext_wet = h_ext / CP_HUMID_AIR   # kg /(m**2 * s)
        R_ext_wet = 1 / h_ext_wet  # m**2 * s / kg
        R_int_wet = Q_(0, 'm ** 2 * s / kg')  # xi = 0 (boiling refrigerant)
        R_wet = R_ext_wet + R_int_wet
        U_wet = 1 / R_wet  # kg /(m**2 * s)
        A_ext = self.hex_core.ext.geo.A
        UA_ext_wet = U_wet * A_ext  # kg / s
        return UA_ext_wet.to('kg / s')

    def __NTU__(self, UA_ext_wet: Quantity) -> float:
        """Calculates the Number of Transfer Units (NTU)."""
        NTU = UA_ext_wet / self.m_dot_air
        return NTU.m

    def __eps__(self, h_ext: Quantity) -> float:
        """Calculates the enthalpy-based heat transfer effectiveness of the
        air-cooling coil.
        """
        UA_ext_wet = self.__UA_ext_wet__(h_ext)
        NTU = self.__NTU__(UA_ext_wet)
        eps = 1 - np.exp(-NTU)
        return eps

    def __ADP__(self, air_out: HumidAir) -> tuple[HumidAir, Quantity, Quantity]:
        """Calculates the resulting Apparatus Dew Point (ADP) of the cooling coil
        that follows from the given state of air entering at the cooling coil
        inlet and the calculated state of air leaving at the cooling coil outlet.
        Also returns the contact factor and bypass factor of the cooling coil
        that go with these air states.
        """
        p = AirConditioningProcess(
            air_in=self.air_in,
            air_out=air_out
        )
        try:
            ADP = p.ADP
        except ValueError:
            warnings.warn(
                message=(
                    "No intersection with saturation line for entering "
                    f"refrigerant temperature {self.rfg_in.T.to('degC')}. "
                    "ADP temperature is set to dew point temperature of "
                    "leaving air."
                ),
                category=DXAirCoolingCoilWarning
            )
            ADP = HumidAir(Tdb=air_out.Tdp, RH=Q_(100, 'pct'))
        CF = p.beta
        BF = 1 - CF
        return ADP, CF, BF

    def __W_air_out_min__(self) -> tuple[Quantity, Quantity]:
        """Returns the lowest absolute humidity ratio that is possible at
        the currently set dry-bulb temperature of the air leaving the DX-coil
        (`self.T_air_out`). Also returns the ADP temperature that corresponds
        with this minimum humidity ratio.
        """

        def _W_air_out(T_adp: float) -> float:
            """Given `T_air_in` and `W_air_in` from `self.air_in`, and also
            `self.T_air_out` and `T_adp`, returns the corresponding `W_air_out`
            (following the definition of the contact factor of an air-cooling
            coil).
            """
            T_adp = Q_(T_adp, 'degC').to('K')
            W_adp = HumidAir(Tdb=T_adp, RH=Q_(100, 'pct')).W
            k1 = (self.air_in.W - W_adp) / (self.air_in.Tdb - T_adp)
            k2 = self.T_air_out - self.air_in.Tdb
            W_air_out = self.air_in.W + k1 * k2
            return W_air_out.to('g / kg').m

        # Find the ADP temperature for which the humidity ratio at the
        # leaving air temperature `T_air_out` is minimal. The possible maximum
        # ADP temperature can certainly never be higher than the cooling coil
        # leaving air temperature.
        T_adp_min = -30.0
        T_adp_max = self.T_air_out.to('degC').m
        r = optimize.minimize_scalar(_W_air_out, bracket=(T_adp_min, T_adp_max))
        T_adp_r = r.x
        W_air_out_min = _W_air_out(T_adp_r)
        return Q_(W_air_out_min, 'g / kg'), Q_(T_adp_r, 'degC')

    def __calculate_operating_state__(
        self,
        i_max: int = 50,
        eps_tol: Quantity = Q_(0.005, 'frac')
    ) -> None:
        """Calculates the resulting operating state of the DX air cooling coil
        under the operating conditions set in `__init__`.

        We know:
        - the state of air entering the cooling coil
        - the mass flow rate of air through the cooling coil
        - the setpoint temperature of air leaving the cooling coil.
        - the state of refrigerant entering the cooling coil
        - the state the refrigerant must have when leaving the cooling coil
        (saturated vapor)

        To determine the resulting operating state, we use an iterative solving
        technique:
        1. We initially guess the state (humidity) of the air leaving the cooling
        coil. This allows us to determine the heat transfer rate to the air
        stream in the cooling coil. The theoretic maximum heat transfer rate is
        known from the inputs to the cooling coil model. Now, an initial value
        for the effectiveness of the cooling coil can be determined.

        2. The effectiveness of the cooling coil is also determined by the heat
        transfer characteristics of the cooling coil's heat exchanger core.
        For this, we first need to calculate the mean states of air and
        refrigerant and the refrigerant mass flow rate. The needed mass flow
        rate of refrigerant can be determined since we have imposed that the
        leaving refrigerant must be saturated vapor. Now, we can calculate the
        air-side convection heat transfer coefficient, and with it, a new value
        for the effectiveness can be determined.

        3. Should we have guessed the leaving air state directly right, both
        effectiveness values would be the same. We check the deviation between
        the two values. If the deviation is still too big, we re-calculate the
        heat transfer rate and leaving air state, and then we go back to step 2.
        This procedure is repeated until the deviation has become small enough
        (i.e., smaller than the set tolerance).

        """
        # Before solving for the actual operating state of the cooling coil, get
        # the minimum absolute humidity ratio that the leaving air can attain
        # when its temperature is controlled to stay at `self.T_air_out`.
        W_air_out_min, _ = self.__W_air_out_min__()
        # Initially assume that the air leaving the cooling coil is at its
        # setpoint temperature and is fully saturated:
        air_out = HumidAir(Tdb=self.T_air_out, RH=Q_(100, 'pct'))
        # From this assumption, we can get an initial guess of the heat
        # transfer rate from the air to the refrigerant and of the effectiveness
        # `eps` of the cooling coil:
        Q = self.m_dot_air * (self.air_in.h - air_out.h)
        eps = Q.to('W').m / self.Q_max.to('W').m
        for i in range(i_max):
            # Calculate the mean state of air and refrigerant along the heat
            # transfer surface:
            air_mean = self.__mean_air__(self.air_in, air_out, self.rfg_in, self.rfg_out)
            rfg_mean = self.__mean_refrigerant__(self.rfg_in, self.rfg_out)
            # Based on the heat transfer rate from air to refrigerant, calculate
            # the required mass flow rate of refrigerant so that refrigerant
            # leaves the cooling coil as saturated vapor:
            m_dot_rfg = Q / (self.rfg_out.h - self.rfg_in.h)
            # Calculate the air-side convection heat transfer coefficient:
            h_ext = self.__h_ext__(air_mean, rfg_mean, m_dot_rfg)
            # Get a new value for `eps`:
            eps_new = self.__eps__(h_ext)
            # Check convergence:
            dev = eps_new - eps
            if abs(dev) < eps_tol.to('frac').m:
                # Absolute humidity of `air_out` cannot be lower than `W_air_out_min`.
                # If that's the case, the controller will maintain the leaving air
                # state at dry-bulb temperature `self.T_air_out` and absolute
                # humidity ratio `W_air_out_min`.
                if air_out.W.to('g / kg') < W_air_out_min:
                    self.air_out = HumidAir(Tdb=self.T_air_out, W=W_air_out_min)
                    self.Q = self.m_dot_air * (self.air_in.h - self.air_out.h)
                    self.m_dot_rfg = self.Q / (self.rfg_out.h - self.rfg_in.h)
                    self.eps = self.Q / self.Q_max
                    Q_sen = self.m_dot_air * CP_HUMID_AIR * (self.air_in.Tdb - self.air_out.Tdb)
                    self.SHR = Q_sen.to('W') / self.Q.to('W')
                    self.ADP, self.CF, self.BF = self.__ADP__(self.air_out)
                else:
                    self.eps = Q_(eps, 'frac')
                    self.Q = Q
                    self.m_dot_rfg = m_dot_rfg
                    self.air_out = air_out
                    Q_sen = self.m_dot_air * CP_HUMID_AIR * (self.air_in.Tdb - self.air_out.Tdb)
                    self.SHR = Q_sen.to('W') / self.Q.to('W')
                    self.ADP, self.CF, self.BF = self.__ADP__(self.air_out)
                break
            # If `dev` not within tolerance band: repeat loop with the new value
            # of `eps`:
            eps = eps_new
            # Re-calculate `Q` for next loop:
            Q = eps * self.Q_max
            # Re-calculate `air_out` for next loop:
            h_air_out = self.air_in.h - Q / self.m_dot_air
            air_out = HumidAir(Tdb=self.T_air_out, h=h_air_out)
        else:
            raise ValueError(
                "Full operating state could not be determined "
                f"after {i_max} iterations."
            )
