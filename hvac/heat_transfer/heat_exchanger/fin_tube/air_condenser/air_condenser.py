"""
AIR CONDENSER MODEL
Plain fin-tube counter-flow heat exchanger.

Known inputs:
- State of air entering the condenser.
- Mass flow rate of air entering the condenser.
- State of refrigerant entering the condenser.
- Mass flow rate of refrigerant entering the condenser.

Purpose:
Determine the state of the air and the state of the refrigerant leaving the
condenser. Once these states have been determined, other operational
characteristics are also determined.
"""
import numpy as np
from scipy import optimize
from hvac import Quantity
from hvac.fluids import FluidState, HumidAir, CP_HUMID_AIR, CoolPropError
from hvac.heat_transfer.heat_exchanger.fin_tube import PlainFinTubeHeatExchangerCore
from hvac.heat_transfer.heat_exchanger.eps_ntu import CounterFlowHeatExchanger
from hvac.logging import ModuleLogger

Q_ = Quantity
logger = ModuleLogger.get_logger(__name__)
logger.setLevel(ModuleLogger.ERROR)


class CondenserError(Exception):
    pass


class DesuperheatingError(CondenserError):
    pass


class CondensingError(CondenserError):
    pass


class SubcoolingError(CondenserError):
    pass


def _get_lmtd(
    air_in: HumidAir,
    air_out: HumidAir,
    rfg_in: FluidState,
    rfg_out: FluidState
) -> Quantity:
    """
    Calculates the LMTD of a heat exchanger.
    """
    dT_in = rfg_out.T - air_in.Tdb
    dT_out = rfg_in.T - air_out.Tdb
    dT_max = max(dT_in, dT_out)
    dT_min = min(dT_in, dT_out)
    if dT_min == dT_max:
        return dT_min
    if dT_min.m <= 0.0:
        logger.debug(
            f"Calculated LMTD. `dT_min` was {dT_min.to('K'):~P.3g}, but cannot "
            f"be zero or negative. It has been changed to a positive value near "
            f"zero."
        )
        dT_min = Q_(1.e-12, 'K')
    lmtd = (dT_max - dT_min) / np.log(dT_max / dT_min)
    return lmtd


class DesuperheatingRegion:

    def __init__(
        self,
        W_fro: Quantity,
        H_fro: Quantity,
        S_trv: Quantity,
        S_lon: Quantity,
        D_int: Quantity,
        D_ext: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity
    ) -> None:
        self.core = PlainFinTubeHeatExchangerCore(
            L1=W_fro,
            L3=H_fro,
            S_t=S_trv,
            S_l=S_lon,
            D_i=D_int,
            D_o=D_ext,
            t_f=t_fin,
            N_f=N_fin,
            k_fin=k_fin
        )
        # Known parameters:
        self.rfg_in: FluidState | None = None   # discharge gas from compressor
        self.rfg_out: FluidState | None = None  # saturated vapor
        self.rfg_m_dot: Quantity | None = None
        self.air_m_dot: Quantity | None = None

        # Unknown parameters to be solved for:
        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.L_flow: Quantity | None = None

    def __fun__(
        self,
        L_flow: float,
        air_mean: HumidAir,
        rfg_mean: FluidState,
        counter: list[int]
    ) -> float:
        """
        Calculates a new value for the flow length of the desuperheating region
        and returns the deviation between this new value and the current value
        `L_flow`.
        """
        # Set the parameters on the heat exchanger core needed to determine
        # the overall heat transfer conductance of the desuperheating region:
        L_flow = Q_(L_flow, 'mm')
        self.core.L2 = L_flow
        self.core.m_dot_ext = self.air_m_dot
        self.core.m_dot_int = self.rfg_m_dot
        self.core.ext.fluid_mean = air_mean
        self.core.int.fluid_mean = rfg_mean

        # Determine a new value the flow length:
        # First, the convection heat transfer coefficient `h_int` on the
        # refrigerant side of the heat exchanger is determined.
        # Then, the total thermal resistance `R_int` between the refrigerant and
        # the heat exchanger wall is determined based on the known heat rate
        # `self.Q_dot` rejected by the refrigerant stream in the desuperheating
        # region.
        # From these, the internal surface area `A_int` of the heat exchanger
        # can be determined.
        # Finally, in the function `_get_flow_length` the flow length of the
        # desuperheating region is calculated that corresponds with `A_int`.
        h_int = self.core.int.h
        R_int = (rfg_mean.T - self.core.T_wall) / self.Q_dot
        A_int = 1 / (R_int * h_int)
        L_flow_new = self._get_flow_length(A_int.to('m ** 2'))

        dev = (L_flow_new - L_flow).to('mm')
        # Should we have passed the correct flow length `L_flow` immediately,
        # the deviation with the new value would be zero.

        i = counter[0]
        logger.debug(
            f"Desuperheating region/Iteration {i + 1}: "
            f"Flow length deviation = {dev:~P.3f}."
        )
        counter[0] += 1
        return dev.magnitude

    def solve(
        self,
        L_flow_max: Quantity,
        x_tol: float = 0.1,
        r_tol: float = 0.01,
        i_max: int = 20
    ) -> Quantity:
        """
        Returns the flow length of the desuperheating region such that the heat
        transfer rate through the heat exchanger balances with the heat rate
        being rejected by the refrigerant stream (and, at the same time, being
        absorbed by the air stream in the desuperheating region).
        """
        # Determine the state of air leaving the desuperheating region:
        h_air_out = self.air_in.h + self.Q_dot / self.air_m_dot
        self.air_out = HumidAir(h=h_air_out, W=self.air_in.W)

        # Determine the mean states of air and refrigerant in the desuperheating
        # region:
        air_mean, rfg_mean = self._get_fluid_mean_states()

        # Initialize iteration counter and set maximum number of iterations:
        counter = [0]

        # Within the given bracket, find the flow length of the desuperheating
        # region for which function `__fun__` returns a deviation near zero.
        try:
            sol = optimize.root_scalar(
                self.__fun__,
                args=(air_mean, rfg_mean, counter),
                method='brentq',
                bracket=(1.0, L_flow_max.to('mm').m),
                xtol=x_tol,
                rtol=r_tol,
                maxiter=i_max
            )
        except ValueError:
            raise DesuperheatingError(
                "The refrigerant cannot be desuperheated into saturated vapor"
                "in the condenser under the current operating conditions."
            ) from None
        self.L_flow = Q_(sol.root, 'mm')
        return self.L_flow

    def _get_fluid_mean_states(self) -> tuple[HumidAir, FluidState]:
        """
        Determines the mean state of the air stream and the refrigerant stream
        in the desuperheating region of the condenser.
        """
        T_air_avg = (self.air_in.Tdb + self.air_out.Tdb) / 2
        air_avg = HumidAir(Tdb=T_air_avg, W=self.air_in.W)

        Rfg = self.rfg_in.fluid
        P_cnd = self.rfg_in.P
        T_rfg_avg = (self.rfg_in.T + self.rfg_out.T) / 2
        rfg_avg = Rfg(P=P_cnd, T=T_rfg_avg)

        C_air = self.air_m_dot * air_avg.cp
        C_rfg = self.rfg_m_dot * rfg_avg.cp
        C_max = max(C_air, C_rfg)
        C_min = min(C_air, C_rfg)
        C_rat = C_min / C_max

        if C_rat >= 0.5:
            return air_avg, rfg_avg
        else:
            lmtd = _get_lmtd(
                self.air_in, self.air_out,
                self.rfg_in, self.rfg_out
            )
            if C_max == C_rfg:
                T_air_avg = rfg_avg.T - lmtd
                air_avg = HumidAir(Tdb=T_air_avg, W=self.air_in.W)
                return air_avg, rfg_avg
            else:
                T_rfg_avg = air_avg.Tdb + lmtd
                rfg_avg = Rfg(T=T_rfg_avg, P=P_cnd)
                return air_avg, rfg_avg

    def _get_flow_length(self, A_int: Quantity) -> Quantity:
        """
        Calculates the flow length that corresponds with the interior heat
        transfer surface area `A_int`.
        """
        D_i = self.core.int.geometry.D_i
        L1 = self.core.int.geometry.L1
        L3 = self.core.int.geometry.L3
        S_trv = self.core.int.geometry.S_t
        S_lon = self.core.int.geometry.S_l
        num = A_int / (np.pi * D_i * L1) - 0.5
        den = L3 / (S_trv * S_lon) - 1 / (2 * S_lon)
        L2 = num / den
        return L2.to('mm')

    @property
    def Q_dot(self) -> Quantity:
        """
        Returns the heat transfer rate in the desuperheating region of the
        condenser.
        """
        Q_dot = self.rfg_m_dot * (self.rfg_in.h - self.rfg_out.h)
        return Q_dot

    @property
    def Q_dot_max(self) -> Quantity:
        """
        Returns the theoretically maximum heat transfer rate between the
        air and refrigerant stream in the desuperheating region of the
        condenser.
        """
        air_mean, rfg_mean = self._get_fluid_mean_states()
        C_air = self.air_m_dot * air_mean.cp
        C_rfg = self.rfg_m_dot * rfg_mean.cp
        C_min = min(C_air, C_rfg)
        Q_dot_max = C_min * (self.rfg_in.T - self.air_in.Tdb)
        return Q_dot_max

    @property
    def dP_air(self) -> Quantity:
        """
        Returns the air-side pressure drop along the desuperheating region of
        the condenser.
        """
        dP_air = self.core.ext.get_pressure_drop(self.air_in, self.air_out)
        return dP_air


class CondensingRegion:

    def __init__(
        self,
        W_fro: Quantity,
        H_fro: Quantity,
        S_trv: Quantity,
        S_lon: Quantity,
        D_int: Quantity,
        D_ext: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity
    ) -> None:
        self.core = PlainFinTubeHeatExchangerCore(
            L1=W_fro,
            L3=H_fro,
            S_t=S_trv,
            S_l=S_lon,
            D_i=D_int,
            D_o=D_ext,
            t_f=t_fin,
            N_f=N_fin,
            k_fin=k_fin,
            condensing=True
        )
        # Known parameters:
        self.rfg_in: FluidState | None = None   # saturated vapor
        self.rfg_out: FluidState | None = None  # saturated liquid
        self.rfg_m_dot: Quantity | None = None
        self.air_m_dot: Quantity | None = None

        # Unknown parameters to be solved for:
        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.L_flow: Quantity | None = None

    def __fun__(
        self,
        L_flow: float,
        air_mean: HumidAir,
        rfg_mean: FluidState,
        counter: list[int]
    ) -> float:
        """
        Calculates a new value for the flow length of the condensing region
        and returns the deviation between this new value and the current value
        `L_flow`.
        """
        # Set the parameters on the heat exchanger core needed to determine
        # the overall heat transfer conductance of the condensing region:
        L_flow = Q_(L_flow, 'mm')
        self.core.L2 = L_flow
        self.core.m_dot_ext = self.air_m_dot
        self.core.m_dot_int = self.rfg_m_dot
        self.core.ext.fluid_mean = air_mean
        self.core.int.fluid_mean = rfg_mean

        # Determine a new value the flow length:
        # First, the convection heat transfer coefficient `h_int` on the
        # refrigerant side of the heat exchanger is determined.
        # Then, the total thermal resistance `R_int` between the refrigerant and
        # the heat exchanger wall is determined based on the known heat rate
        # `self.Q_dot` rejected by the refrigerant stream in the condensing
        # region.
        # From these, the internal surface area `A_int` of the heat exchanger
        # can be determined.
        # Finally, in the function `_get_flow_length` the flow length of the
        # condensing region is calculated that corresponds with `A_int`.
        h_int = self.core.int.h
        R_int = (rfg_mean.T - self.core.T_wall) / self.Q_dot
        A_int = 1 / (R_int * h_int)
        L_flow_new = self._get_flow_length(A_int.to('m ** 2'))

        dev = (L_flow_new - L_flow).to('mm')
        # Should we have passed the correct flow length `L_flow` immediately,
        # the deviation with the new value would be zero.

        i = counter[0]
        logger.debug(
            f"Condensing region/Iteration {i + 1}: "
            f"Flow length deviation = {dev:~P.3f}."
        )
        counter[0] += 1
        return dev.magnitude

    def solve(
        self,
        L_flow_max: Quantity,
        x_tol: float = 0.1,
        r_tol: float = 0.01,
        i_max: int = 20
    ) -> Quantity:
        """
        Returns the flow length of the condensing region such that the heat
        transfer rate through the heat exchanger balances with the heat rate
        being rejected by the refrigerant stream (and, at the same time, being
        absorbed by the air stream in the condensing region.)
        """
        # Determine the state of air leaving the condensing region:
        h_air_out = self.air_in.h + self.Q_dot / self.air_m_dot
        self.air_out = HumidAir(h=h_air_out, W=self.air_in.W)

        # Determine the mean states of air and refrigerant in the condensing
        # region:
        air_mean, rfg_mean = self._get_fluid_mean_states()

        # Initialize iteration counter:
        counter = [0]

        # Within the given bracket, find the flow length of the condensing
        # region for which function `__fun__` returns a deviation near zero.
        try:
            sol = optimize.root_scalar(
                self.__fun__,
                args=(air_mean, rfg_mean, counter),
                method='brentq',
                bracket=(1.0, L_flow_max.to('mm').m),
                xtol=x_tol,
                rtol=r_tol,
                maxiter=i_max
            )
        except ValueError:
            raise CondensingError(
                "The refrigerant cannot be fully condensed in the "
                "condenser under the current operating conditions."
            ) from None
        self.L_flow = Q_(sol.root, 'mm')
        return self.L_flow

    def _get_fluid_mean_states(self) -> tuple[HumidAir, FluidState]:
        """
        Determines the mean state of the air stream and the refrigerant stream
        in the condensing region of the condenser.
        """
        Rfg = self.rfg_in.fluid
        P_cnd = self.rfg_in.P
        rfg_avg = Rfg(P=P_cnd, x=Q_(0.5, 'frac'))
        lmtd = _get_lmtd(
            self.air_in, self.air_out,
            rfg_avg, rfg_avg
        )
        T_air_avg = rfg_avg.T - lmtd
        air_avg = HumidAir(Tdb=T_air_avg, W=self.air_in.W)
        return air_avg, rfg_avg

    @property
    def Q_dot(self) -> Quantity:
        """
        Returns the heat transfer rate in the desuperheating region of the
        condenser.
        """
        Q_dot = self.rfg_m_dot * (self.rfg_in.h - self.rfg_out.h)
        return Q_dot

    def _get_flow_length(self, A_int: Quantity) -> Quantity:
        """
        Calculates the flow length that corresponds with the interior heat
        transfer surface area `A_int`.
        """
        D_i = self.core.int.geometry.D_i
        L1 = self.core.int.geometry.L1
        L3 = self.core.int.geometry.L3
        S_trv = self.core.int.geometry.S_t
        S_lon = self.core.int.geometry.S_l
        num = A_int / (np.pi * D_i * L1) - 0.5
        den = L3 / (S_trv * S_lon) - 1 / (2 * S_lon)
        L2 = num / den
        return L2.to('mm')

    @property
    def Q_dot_max(self) -> Quantity:
        """
        Returns the theoretically maximum heat transfer rate between the
        air and refrigerant stream in the condensing region of the
        condenser.
        """
        air_mean, _ = self._get_fluid_mean_states()
        C_min = self.air_m_dot * air_mean.cp
        Q_dot_max = C_min * (self.rfg_in.T - self.air_in.Tdb)
        return Q_dot_max

    @property
    def dP_air(self) -> Quantity:
        """
        Returns the air-side pressure drop along the condensing region of
        the condenser.
        """
        dP_air = self.core.ext.get_pressure_drop(self.air_in, self.air_out)
        return dP_air


class SubcoolingRegion:

    def __init__(
        self,
        W_fro: Quantity,
        H_fro: Quantity,
        S_trv: Quantity,
        S_lon: Quantity,
        D_int: Quantity,
        D_ext: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity
    ) -> None:
        self.core = PlainFinTubeHeatExchangerCore(
            L1=W_fro,
            L3=H_fro,
            S_t=S_trv,
            S_l=S_lon,
            D_i=D_int,
            D_o=D_ext,
            t_f=t_fin,
            N_f=N_fin,
            k_fin=k_fin
        )
        # Known parameters:
        self.rfg_in: FluidState | None = None  # saturated liquid
        self.rfg_m_dot: Quantity | None = None
        self.air_in: HumidAir | None = None
        self.air_m_dot: Quantity | None = None

        # Unknown parameters to be solved:
        self.rfg_out: FluidState | None = None
        self.air_out: HumidAir | None = None
        self.Q_dot: Quantity | None = None
        self.L_flow: Quantity | None = None

    def _get_rfg_out(self, T_rfg_out: Quantity) -> FluidState:
        Rfg = self.rfg_in.fluid
        P_cnd = self.rfg_in.P
        try:
            rfg_out = Rfg(P=P_cnd, T=T_rfg_out)
        except CoolPropError:
            # This exception will be raised when the refrigerant is not
            # fully condensed and still a two-phase mixture.
            rfg_out = Rfg(P=P_cnd, x=Q_(0, 'frac'))
        return rfg_out

    def __fun__(
        self,
        T_rfg_out: Quantity,
        L_flow: Quantity,
        counter: list[int]
    ) -> tuple[Quantity, Quantity]:
        """
        Calculates for the given flow length `L_flow` of the subcooling region,
        a new value for the temperature of the refrigerant leaving the
        condenser, so that the heat transfer rate through the heat exchanger
        would balance with the heat rate rejected by the refrigerant in the
        subcooling region.
        """
        # Current state of the refrigerant leaving the condenser:
        rfg_out = self._get_rfg_out(T_rfg_out)

        # Heat rate rejected by the refrigerant in the subcooling region for
        # the current state of the leaving refrigerant:
        Q_dot = self.rfg_m_dot * (self.rfg_in.h - rfg_out.h)

        # Determine the mean states of air and refrigerant in the subcooling
        # region:
        T_air_out = self.air_in.Tdb + Q_dot / (CP_HUMID_AIR * self.air_m_dot)
        air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
        air_mean, rfg_mean = self._get_fluid_mean_states(air_out, rfg_out)

        # Set heat exchanger core parameters to calculate UA:
        self.core.L2 = L_flow
        self.core.m_dot_ext = self.air_m_dot
        self.core.m_dot_int = self.rfg_m_dot
        self.core.ext.fluid_mean = air_mean
        self.core.int.fluid_mean = rfg_mean

        # Determine a new value for the heat transfer rate through the heat
        # exchanger core of the subcooling region:
        cnt_flw_hex = CounterFlowHeatExchanger(
            C_cold=self.air_m_dot * air_mean.cp,
            C_hot=self.rfg_m_dot * rfg_mean.cp,
            T_cold_in=self.air_in.Tdb,
            T_hot_in=self.rfg_in.T,
            UA=self.core.UA
        )
        Q_dot_new = cnt_flw_hex.Q

        # Determine the new state of the refrigerant leaving the condenser:
        h_rfg_out_new = self.rfg_in.h - Q_dot_new / self.rfg_m_dot
        Rfg = self.rfg_in.fluid
        P_cnd = self.rfg_in.P
        rfg_out_new = Rfg(P=P_cnd, h=h_rfg_out_new)

        # Temperature deviation between the new state and the current state of
        # the refrigerant leaving the condenser:
        dev = (rfg_out_new.T - rfg_out.T).to('K')

        i = counter[0]
        logger.debug(
            f"Subcooling region/Iteration {i + 1}: "
            f"Flow length = {L_flow.to('mm'):~P.3f}. "
            f"Temperature deviation of leaving refrigerant = {dev:~P.3f}."
        )
        counter[0] += 1
        return dev, rfg_out_new.T

    def solve(
        self,
        L_flow: Quantity,
        tol: Quantity = Q_(0.1, 'K'),
        i_max: int = 20
    ) -> HumidAir:
        """
        For the given subcooling flow length `L_flow`, solve for the output
        parameters of the subcooling region (attributes `air_out`, `rfg_out`,
        and `Q_dot`).

        Returns
        -------
        The state of air leaving the subcooling region and entering the
        condensing region of the condenser.
        """
        T_rfg_out_min = self.air_in.Tdb.to('degC')
        T_rfg_out_max = self.rfg_in.T.to('degC')

        counter = [0]

        # Initial guess of the temperature of the refrigerant leaving the
        # condenser:
        T_rfg_out = (T_rfg_out_min.to('K') + T_rfg_out_max.to('K')) / 2

        for i in range(i_max):
            dev, T_rfg_out_new = self.__fun__(T_rfg_out, L_flow, counter)
            if abs(dev) < tol:
                Rfg = self.rfg_in.fluid
                P_cnd = self.rfg_in.P
                self.rfg_out = Rfg(P=P_cnd, T=T_rfg_out_new)
                self.Q_dot = self.rfg_m_dot * (self.rfg_in.h - self.rfg_out.h)
                T_air_out = self.air_in.Tdb + self.Q_dot / (CP_HUMID_AIR * self.air_m_dot)
                self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
                return self.air_out
            T_rfg_out = T_rfg_out_new
        else:
            raise SubcoolingError(
                f"No acceptable solution found after {i_max} iterations."
            ) from None

    def _get_fluid_mean_states(
        self,
        air_out: HumidAir,
        rfg_out: FluidState
    ) -> tuple[HumidAir, FluidState]:
        """
        Determines the mean state of the air stream and the refrigerant stream
        in the subcooling region of the condenser.
        """
        T_air_avg = (self.air_in.Tdb + air_out.Tdb) / 2
        air_avg = HumidAir(Tdb=T_air_avg, W=self.air_in.W)

        T_rfg_avg = (self.rfg_in.T + rfg_out.T) / 2
        rfg_avg = self._get_rfg_out(T_rfg_avg)

        C_air = self.air_m_dot * air_avg.cp
        C_rfg = self.rfg_m_dot * rfg_avg.cp
        C_max = max(C_air, C_rfg)
        C_min = min(C_air, C_rfg)
        C_rat = C_min / C_max

        if C_rat >= 0.5:
            return air_avg, rfg_avg
        else:
            lmtd = _get_lmtd(
                self.air_in, air_out,
                self.rfg_in, rfg_out
            )
            if C_max == C_rfg:
                T_air_avg = rfg_avg.T - lmtd
                air_avg = HumidAir(Tdb=T_air_avg, W=self.air_in.W)
                return air_avg, rfg_avg
            else:
                T_rfg_avg = air_avg.Tdb + lmtd
                rfg_avg = self._get_rfg_out(T_rfg_avg)
                return air_avg, rfg_avg

    @property
    def Q_dot_max(self) -> Quantity:
        """
        Returns the theoretically maximum heat transfer rate between
        air and refrigerant in the subcooling region of the condenser.
        """
        air_mean, rfg_mean = self._get_fluid_mean_states(self.air_out, self.rfg_out)
        C_air = self.air_m_dot * air_mean.cp
        C_rfg = self.rfg_m_dot * rfg_mean.cp
        C_min = min(C_air, C_rfg)
        Q_dot_max = C_min * (self.rfg_in.T - self.air_in.Tdb)
        return Q_dot_max

    @property
    def dP_air(self) -> Quantity:
        """
        Returns the air-side pressure drop along the subcooling region of
        the condenser.
        """
        dP_air = self.core.ext.get_pressure_drop(self.air_in, self.air_out)
        return dP_air


class PlainFinTubeCounterFlowAirCondenser:
    """
    Model class for a plain fin-tube counter-flow air condenser.

    Attributes
    ----------
    air_in: HumidAir (known)
    air_m_dot: Quantity (known)
    rfg_in: FluidState (known)
    rfg_m_dot: Quantity (known)
    P_cnd: Quantity
    T_cnd: Quantity
    air_out: HumidAir
    rfg_out: FluidState
    Q_dot: Quantity
    Q_dot_max: Quantity
    eps: Quantity
    air_dP: Quantity
    """
    def __init__(
        self,
        W_fro: Quantity,
        H_fro: Quantity,
        N_rows: int,
        S_trv: Quantity,
        S_lon: Quantity,
        D_int: Quantity,
        D_ext: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)')
    ) -> None:
        """
        Creates the plain fin-tube counter-flow air condenser.

        Parameters
        ----------
        W_fro:
            Width of the frontal area of the evaporator.
        H_fro:
            Height of the frontal area.
        N_rows:
            The number of rows in the evaporator.
        S_trv:
            Transversal pitch, i.e., the spacing between tubes in a single row.
        S_lon:
            Longitudinal pitch, i.e., the spacing between tubes of adjacent rows.
        D_int:
            Inside diameter of the tubes.
        D_ext:
            Outside diameter of the tubes.
        t_fin:
            Thickness of the plain fins.
        N_fin:
            Number of fins per unit length of tube (i.e., the inverse of fin
            density).
        k_fin:
            Thermal conductivity of the fin material. The default value applies
            to aluminum.
        """
        self.desuperheating_region = DesuperheatingRegion(
            W_fro, H_fro, S_trv, S_lon, D_int, D_ext,
            t_fin, N_fin, k_fin
        )
        self.condensing_region = CondensingRegion(
            W_fro, H_fro, S_trv, S_lon, D_int, D_ext,
            t_fin, N_fin, k_fin
        )
        self.subcooling_region = SubcoolingRegion(
            W_fro, H_fro, S_trv, S_lon, D_int, D_ext,
            t_fin, N_fin, k_fin
        )
        self.L_flow = N_rows * S_lon

        # Known parameters:
        self.rfg_in: FluidState | None = None   # discharge gas from compressor
        self.rfg_m_dot: Quantity | None = None
        self.air_in: HumidAir | None = None
        self.air_m_dot: Quantity | None = None
        self.P_cnd: Quantity | None = None
        self.T_cnd: Quantity | None = None

        # Unknown parameters to be solved for:
        self.rfg_out: FluidState | None = None
        self.air_out: HumidAir | None = None
        self.Q_dot: Quantity | None = None
        self.Q_dot_max: Quantity | None = None
        self.eps: Quantity | None = None
        self.dT_sc: Quantity | None = None
        self.air_dP: Quantity | None = None

    def __fun__(self, L_flow_sub: float, counter: list[int]) -> float:
        """
        Calculates the deviation between the calculated total flow length of the
        condenser and its actual total flow length, starting from a given
        flow length for the subcooling region only.
        """
        L_flow_sub = Q_(L_flow_sub, 'mm')
        i = counter[0]

        logger.debug(
            f"Subcooling region/Iteration {i + 1}: "
            f"Try with subcooling flow length {L_flow_sub:~P.3f}."
        )

        # Solve the subcooling region with the given subcooling flow length:
        air_out = self.subcooling_region.solve(L_flow_sub)

        logger.debug(
            f"Subcooling region/Iteration {i + 1}: "
            f"Leaving air temperature = {air_out.Tdb.to('degC'):~P.3f}. "
            "Determine condensing flow length..."
        )

        # Solve the condensing region for its condensing flow length:
        self.condensing_region.air_in = air_out
        L_flow_cnd = self.condensing_region.solve(self.L_flow)

        logger.debug(
            f"Condensing region/Iteration {i + 1}: "
            f"Flow length = {L_flow_cnd.to('mm'):~P.3f}. "
            "Determine desuperheating flow length..."
        )

        # Solve the desuperheating region for its desuperheating flow length:
        self.desuperheating_region.air_in = self.condensing_region.air_out
        L_flow_dsh = self.desuperheating_region.solve(self.L_flow)

        logger.debug(
            f"Desuperheating region/Iteration {i + 1}: "
            f"Flow length = {L_flow_dsh.to('mm'):~P.3f}."
        )

        # Determine the deviation between the currently calculated total flow
        # length of the condenser and its actual total flow length:
        L_flow_new = L_flow_sub + L_flow_cnd + L_flow_dsh
        dev = (L_flow_new - self.L_flow).to('mm')

        logger.debug(
            f"Condenser/Iteration {i + 1}: "
            f"Deviation with actual condenser flow length = {dev:~P.3f}."
        )

        counter[0] += 1
        return dev.m

    def solve(
        self,
        air_in: HumidAir,
        air_m_dot: Quantity,
        rfg_in: FluidState,
        rfg_m_dot: Quantity,
        L_flow_sub_min: Quantity = Q_(0.1, 'mm'),
        x_tol: float = 0.1,
        r_tol: float = 0.01,
        i_max: int = 20
    ) -> tuple[HumidAir, FluidState]:
        """
        Solves for the operating state of the condenser under the given
        operating conditions.

        The operating state is fully determined once the states of the air
        and the refrigerant leaving the condenser are determined.

        Parameters
        ----------
        air_in:
            State of the air entering the condenser.
        air_m_dot:
            Mass flow rate of air through the condenser.
        rfg_in:
            State of the refrigerant entering the condenser.
        rfg_m_dot:
            Mass flow rate of refrigerant through the condenser.
        L_flow_sub_min: optional
            The smallest subcooling flow length to try with. By default, 1 mm.
        x_tol: optional
            Absolute tolerance for root-finding algorithm
            (`scipy.optimize.root_scalar` with 'brentq' method).
        r_tol: optional
            Relative tolerance for root-finding algorithm.
        i_max: optional
            Maximum number of iterations for root-finding algorithm

        Returns
        -------
        The states of the air and refrigerant leaving the condenser.

        Notes
        -----
        The condenser has three regions. First, the refrigerant is desuperheated.
        Next, the refrigerant condenses. Finally, the refrigerant is subcooled.
        We know the mass flow rate of refrigerant, and we know the state of air
        and refrigerant entering the subcooling region, respectively the
        desuperheating region.
        We can choose a subcooling flow length and solve for the state of air
        leaving the subcooling region and entering the condensing region.
        The states of refrigerant at the entrance of the desuperheating region,
        the condensing region, and the subcooling region are all known. Since we
        also know the refrigerant mass flow rate, we can determine the heat
        transfer in the condensing and the desuperheating region. With these, we
        can determine the flow length of the condensing region and of the
        desuperheating region, and also the state of air leaving the condensing
        region and the desuperheating region.
        The sum of the chosen subcooling flow length and the flow lengths of
        the condensing and the desuperheating region calculated from it, must
        ultimately be equal to the actual total flow length of the condenser.
        To find the subcooling flow length, for which the sum of the flow lengths
        of the three regions equals the actual total condenser flow length, we
        use a root-finding algorithm.
        """
        self._set_air_in(air_in)
        self._set_air_m_dot(air_m_dot)
        self._set_rfg_in(rfg_in)
        self._set_rfg_m_dot(rfg_m_dot)

        # Find the flow length of the subcooling region so that the sum of the
        # flow lengths of the subcooling, condensing, and desuperheating region
        # equals the actual, total flow length of the condenser:
        logger.debug(
            'Solving for the operating state of the condenser...'
        )
        L_flow_sub_min = L_flow_sub_min.to('mm').m
        L_flow_sub_max = self.L_flow.to('mm').m
        counter = [0]
        try:
            sol = optimize.root_scalar(
                self.__fun__,
                args=(counter,),
                method='brentq',
                bracket=(L_flow_sub_min, L_flow_sub_max),
                xtol=x_tol,  # mm
                rtol=r_tol,  # 1 %
                maxiter=i_max
            )
        except ValueError:
            message = (
                "The refrigerant cannot be subcooled in the condenser "
                "under the current operating conditions."
            )
            logger.error(message)
            raise CondenserError(message) from None
        except (SubcoolingError, CondensingError, DesuperheatingError) as err:
            logger.error(f"{type(err).__name__}: {err}")
            raise err
        else:
            logger.debug(
                f"Calculation finished: {sol.flag}"
            )
            self.subcooling_region.L_flow = Q_(sol.root, 'mm')
            self.rfg_out = self.subcooling_region.rfg_out
            self.air_out = self.desuperheating_region.air_out
            self.Q_dot = (
                self.desuperheating_region.Q_dot
                + self.condensing_region.Q_dot
                + self.subcooling_region.Q_dot
            )
            self.Q_dot_max = (
                self.desuperheating_region.Q_dot_max
                + self.condensing_region.Q_dot_max
                + self.subcooling_region.Q_dot_max
            )
            self.eps = self.Q_dot / self.Q_dot_max
            self.dT_sc = self.T_cnd - self.rfg_out.T
            self.air_dP = (
                self.desuperheating_region.dP_air
                + self.condensing_region.dP_air
                + self.subcooling_region.dP_air
            )
            return self.air_out, self.rfg_out

    def _set_air_in(self, air_in: HumidAir) -> None:
        self.air_in = air_in
        self.subcooling_region.air_in = air_in

    def _set_air_m_dot(self, air_m_dot: Quantity) -> None:
        self.air_m_dot = air_m_dot
        self.subcooling_region.air_m_dot = air_m_dot
        self.condensing_region.air_m_dot = air_m_dot
        self.desuperheating_region.air_m_dot = air_m_dot

    def _set_rfg_in(self, rfg_in: FluidState) -> None:
        self.rfg_in = rfg_in
        self.desuperheating_region.rfg_in = rfg_in
        self.P_cnd = self.rfg_in.P
        Rfg = rfg_in.fluid
        self.desuperheating_region.rfg_out = Rfg(P=self.P_cnd, x=Q_(1, 'frac'))
        self.condensing_region.rfg_in = self.desuperheating_region.rfg_out
        self.condensing_region.rfg_out = Rfg(P=self.P_cnd, x=Q_(0, 'frac'))
        self.T_cnd = self.condensing_region.rfg_out.T
        self.subcooling_region.rfg_in = self.condensing_region.rfg_out

    def _set_rfg_m_dot(self, rfg_m_dot: Quantity) -> None:
        self.rfg_m_dot = rfg_m_dot
        self.desuperheating_region.rfg_m_dot = rfg_m_dot
        self.condensing_region.rfg_m_dot = rfg_m_dot
        self.subcooling_region.rfg_m_dot = rfg_m_dot
