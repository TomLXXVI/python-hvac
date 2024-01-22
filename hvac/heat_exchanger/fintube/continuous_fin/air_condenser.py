import math
from abc import ABC, abstractmethod
from scipy import optimize
from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import HumidAir, FluidState, CoolPropError, CP_HUMID_AIR
import hvac.heat_transfer.forced_convection.internal_flow as single_phase_flow
import hvac.heat_transfer.condensation.flow_condensation as condensing_flow
from hvac.heat_exchanger.general.eps_ntu import CounterFlowHeatExchanger
from .geometry import ContinuousFinStaggeredTubeBank
from .air_to_water import ExternalSurface


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


class HeatExchangerCore:
    # Continuous plain fin, staggered tube bank.

    def __init__(
        self,
        width: Quantity,
        height: Quantity,
        num_rows: int,
        pitch_trv: Quantity,
        pitch_lon: Quantity,
        d_o: Quantity,
        d_i: Quantity,
        t_fin: Quantity,
        fin_density: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
        d_r: Quantity | None = None,
        condensing: bool = False,
        A_min_tot: Quantity | None = None,
        N_rows_tot: int | None = None
    ) -> None:
        self.A_min_tot = A_min_tot
        self.N_rows_tot = N_rows_tot
        self.geometry = ContinuousFinStaggeredTubeBank(
            width, height, num_rows, pitch_trv,
            pitch_lon, d_o, d_i, t_fin, fin_density,
            k_fin, d_r
        )
        if condensing:
            self.internal = CondensingPhaseInternalSurface(self)
        else:
            self.internal = SinglePhaseInternalSurface(self)
        self.external = ExternalSurface(self)

        self.T_wall: Quantity | None = None
        self.h_int: Quantity | None = None
        self.R_int: Quantity | None = None
        self.h_ext: Quantity | None = None
        self.R_ext: Quantity | None = None
        self.R_tot: Quantity | None = None
        self.UA: Quantity | None = None
        self.eta_surf: float | None = None
        self.dP_ext: Quantity | None = None

    def hex_properties(
        self,
        air_mean: HumidAir,
        air_in: HumidAir,
        air_out: HumidAir,
        rfg_mean: FluidState,
        air_m_dot: Quantity,
        rfg_m_dot: Quantity,
    ) -> dict[str, Quantity | float]:
        self.internal.m_dot = rfg_m_dot
        self.internal.rfg = rfg_mean
        self.external.m_dot = air_m_dot
        self.external.air = air_mean
        self.T_wall = self.__find_wall_temperature()
        self.h_int = self.internal.heat_transfer_coeff(self.T_wall)
        self.R_int = self.internal.thermal_resistance(self.h_int)
        self.h_ext = self.external.heat_transfer_coeff(self.T_wall)
        self.R_ext = self.external.thermal_resistance(self.h_ext)
        self.R_tot = self.R_int + self.R_ext
        # Note: The thermal conduction resistance of the heat exchanger body
        # is ignored, and fouling is not being taken into account.
        self.UA = 1 / self.R_tot
        self.dP_ext = self.external.pressure_drop(air_in.rho, air_out.rho, self.T_wall)
        self.eta_surf = self.external.eta_surf(self.h_ext)
        return {
            'h_int': self.h_int,
            'h_ext': self.h_ext,
            'R_int': self.R_int,
            'R_ext': self.R_ext,
            'R_tot': self.R_tot,
            'UA': self.UA,
            'eta_surf': self.eta_surf,
            'dP_ext': self.dP_ext
        }

    def __find_wall_temperature(self) -> Quantity:

        def __eq__(T_wall_: float) -> float:
            T_wall_ = Q_(T_wall_, 'K')
            h_int = self.internal.heat_transfer_coeff(T_wall_)
            R_int = self.internal.thermal_resistance(h_int)
            h_ext = self.external.heat_transfer_coeff(T_wall_)
            R_ext = self.external.thermal_resistance(h_ext)
            T_int_ = self.internal.rfg.T.to('K')
            T_ext_ = self.external.air.Tdb.to('K')
            n = T_int_ / R_int + T_ext_ / R_ext
            d = 1 / R_int + 1 / R_ext
            T_wall_new = n / d
            dev = (T_wall_new - T_wall_).to('K').m
            return dev

        T_int = self.internal.rfg.T.to('K').m   # hot fluid
        T_ext = self.external.air.Tdb.to('K').m  # cold fluid
        if T_int > T_ext:
            sol = optimize.root_scalar(
                __eq__,
                bracket=(T_ext + 1.e-12, T_int - 1.e-12)
            )
            T_wall = Q_(sol.root, 'K')
        else:
            T_wall = Q_(T_int, 'K')
        return T_wall


class InternalSurface(ABC):

    def __init__(self, parent: HeatExchangerCore):
        self.parent = parent
        self.rfg: FluidState | None = None
        self.m_dot: Quantity | None = None

    def _m_dot_tube(self) -> Quantity:
        # Here it is assumed that every tube in the first row is connected to
        # the supply header.
        A_min = self.parent.A_min_tot
        n_r = self.parent.N_rows_tot
        d_i = self.parent.geometry.d_i
        G = self.m_dot / (A_min / n_r)
        A_tube = math.pi * d_i ** 2 / 4
        m_dot_tube = G * A_tube
        return m_dot_tube

    def thermal_resistance(self, h: Quantity) -> Quantity:
        A_tot = self.parent.geometry.internal.A_tot
        R = 1 / (h * A_tot)
        return R

    @abstractmethod
    def heat_transfer_coeff(self, *args, **kwargs) -> Quantity:
        ...


class SinglePhaseInternalSurface(InternalSurface):

    def heat_transfer_coeff(self, *args, **kwargs) -> Quantity:
        tube = single_phase_flow.CircularTube(
            Di=self.parent.geometry.d_i,
            L=self.parent.geometry.length,
            fluid=self.rfg
        )
        tube.m_dot = self._m_dot_tube()
        h = tube.avg_heat_transfer_coefficient()
        if isinstance(h, tuple):
            h = h[0]  # if laminar flow, assume constant heat flux
        return h.to('W / (m ** 2 * K)')


class CondensingPhaseInternalSurface(InternalSurface):

    def heat_transfer_coeff(self, T_wall: Quantity) -> Quantity:
        tube = condensing_flow.HorizontalTube(
            D=self.parent.geometry.internal.d_h,
            fluid=self.rfg.fluid
        )
        tube.m_dot = self._m_dot_tube()
        h = tube.heat_trf_coeff(
            x=self.rfg.x,
            T_sat=self.rfg.T,
            T_surf=T_wall
        )
        return h.to('W / (m ** 2 * K)')


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
        logger.warning(
            f"Calculated LMTD. `dT_min` was {dT_min.to('K'):~P.3g}, but cannot "
            f"be zero or negative. It has been changed to a positive value near "
            f"zero."
        )
        dT_min = Q_(1.e-12, 'K')
    lmtd = (dT_max - dT_min) / math.log(dT_max / dT_min)
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
        k_fin: Quantity,
        A_min_tot: Quantity,
        N_rows_tot: int
    ) -> None:
        self.core = HeatExchangerCore(
            width=W_fro,
            height=H_fro,
            num_rows=0,
            pitch_trv=S_trv,
            pitch_lon=S_lon,
            d_i=D_int,
            d_o=D_ext,
            t_fin=t_fin,
            fin_density=N_fin,
            k_fin=k_fin,
            A_min_tot=A_min_tot,
            N_rows_tot=N_rows_tot
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
        self.core.geometry.length = L_flow
        self.core.hex_properties(
            air_mean=air_mean,
            air_in=self.air_in,
            air_out=self.air_out,
            rfg_mean=rfg_mean,
            air_m_dot=self.air_m_dot,
            rfg_m_dot=self.rfg_m_dot
        )

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
        R_int = (rfg_mean.T - self.core.T_wall) / self.Q_dot
        A_int = 1 / (R_int * self.core.h_int)
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

    def __fun__new(
        self,
        L_flow: float,
        air_mean: HumidAir,
        rfg_mean: FluidState,
        counter: list[int]
    ) -> float:
        """
        Calculates the heat transfer rate through the heat exchanger core of the
        desuperheating region and returns the deviation between this value and the
        heat rate that the refrigerant must reject to turn from superheated vapor
        into saturated vapor.
        """
        L_flow = Q_(L_flow, 'mm')
        self.core.geometry.length = L_flow
        self.core.hex_properties(
            air_mean=air_mean,
            air_in=self.air_in,
            air_out=self.air_out,
            rfg_mean=rfg_mean,
            air_m_dot=self.air_m_dot,
            rfg_m_dot=self.rfg_m_dot
        )

        cnt_flw_hex = CounterFlowHeatExchanger(
            C_cold=self.air_m_dot * air_mean.cp,
            C_hot=self.rfg_m_dot * rfg_mean.cp,
            T_cold_in=self.air_in.Tdb,
            T_hot_in=self.rfg_in.T,
            UA=self.core.UA
        )

        dev = (cnt_flw_hex.Q - self.Q_dot).to('kW')
        i = counter[0]
        logger.debug(
            f"Desuperheating region/Iteration {i + 1}: "
            f"Heat transfer deviation = {dev:~P.3f}."
        )
        counter[0] += 1
        return dev.magnitude

    def solve(
        self,
        L_flow_max: Quantity,
        x_tol: float = 0.001,
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
                self.__fun__new,
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
        D_i = self.core.geometry.d_i
        L1 = self.core.geometry.width
        L3 = self.core.geometry.height
        S_trv = self.core.geometry.pitch_trv
        S_lon = self.core.geometry.pitch_lon
        num = A_int / (math.pi * D_i * L1) - 0.5
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
        return self.core.dP_ext


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
        k_fin: Quantity,
        A_min_tot: Quantity,
        N_rows_tot: int
    ) -> None:
        self.core = HeatExchangerCore(
            width=W_fro,
            height=H_fro,
            num_rows=0,
            pitch_trv=S_trv,
            pitch_lon=S_lon,
            d_i=D_int,
            d_o=D_ext,
            t_fin=t_fin,
            fin_density=N_fin,
            k_fin=k_fin,
            condensing=True,
            A_min_tot=A_min_tot,
            N_rows_tot=N_rows_tot
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
        self.core.geometry.length = L_flow
        self.core.hex_properties(
            air_mean=air_mean,
            air_in=self.air_in,
            air_out=self.air_out,
            rfg_mean=rfg_mean,
            air_m_dot=self.air_m_dot,
            rfg_m_dot=self.rfg_m_dot
        )

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
        R_int = (rfg_mean.T - self.core.T_wall) / self.Q_dot
        A_int = 1 / (R_int * self.core.h_int)
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

    def __fun__new(
        self,
        L_flow: float,
        air_mean: HumidAir,
        rfg_mean: FluidState,
        counter: list[int]
    ) -> float:
        """
        Calculates the heat transfer rate through the heat exchanger core of the
        condensing region and returns the deviation between this value and the
        heat rate that the refrigerant must reject to turn from saturated vapor
        into saturated liquid.
        """
        L_flow = Q_(L_flow, 'mm')
        self.core.geometry.length = L_flow
        self.core.hex_properties(
            air_mean=air_mean,
            air_in=self.air_in,
            air_out=self.air_out,
            rfg_mean=rfg_mean,
            air_m_dot=self.air_m_dot,
            rfg_m_dot=self.rfg_m_dot
        )

        cnt_flw_hex = CounterFlowHeatExchanger(
            C_cold=self.air_m_dot * air_mean.cp,
            C_hot=Q_(float('inf'), 'W / K'),  # self.rfg_m_dot * rfg_mean.cp,
            T_cold_in=self.air_in.Tdb,
            T_hot_in=self.rfg_in.T,
            UA=self.core.UA
        )

        dev = (cnt_flw_hex.Q - self.Q_dot).to('kW')
        i = counter[0]
        logger.debug(
            f"Condensing region/Iteration {i + 1}: "
            f"Flow length = {L_flow:~P.3f}. "
            f"Heat transfer deviation = {dev:~P.3f}."
        )
        counter[0] += 1
        return dev.magnitude

    def solve(
        self,
        L_flow_max: Quantity,
        x_tol: float = 0.001,
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
                self.__fun__new,
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
        Returns the heat transfer rate in the condensing region of the
        condenser.
        """
        Q_dot = self.rfg_m_dot * (self.rfg_in.h - self.rfg_out.h)
        return Q_dot

    def _get_flow_length(self, A_int: Quantity) -> Quantity:
        """
        Calculates the flow length that corresponds with the interior heat
        transfer surface area `A_int`.
        """
        D_i = self.core.geometry.d_i
        L1 = self.core.geometry.width
        L3 = self.core.geometry.height
        S_trv = self.core.geometry.pitch_trv
        S_lon = self.core.geometry.pitch_lon
        num = A_int / (math.pi * D_i * L1) - 0.5
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
        return self.core.dP_ext


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
        k_fin: Quantity,
        A_min_tot: Quantity,
        N_rows_tot: int
    ) -> None:
        self.core = HeatExchangerCore(
            width=W_fro,
            height=H_fro,
            num_rows=0,
            pitch_trv=S_trv,
            pitch_lon=S_lon,
            d_i=D_int,
            d_o=D_ext,
            t_fin=t_fin,
            fin_density=N_fin,
            k_fin=k_fin,
            A_min_tot=A_min_tot,
            N_rows_tot=N_rows_tot
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
        self.core.geometry.length = L_flow
        self.core.hex_properties(
            air_mean=air_mean,
            air_in=self.air_in,
            air_out=air_out,
            rfg_mean=rfg_mean,
            air_m_dot=self.air_m_dot,
            rfg_m_dot=self.rfg_m_dot
        )

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
        return self.core.dP_ext


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
        # Create the geometry of the whole condenser to determine `A_min` of
        # the whole condenser; we need this, together with the number of rows
        # of the whole condenser, to determine the mass flow rate of
        # refrigerant in 1 tube (see class `InternalSurface`, method
        # `_m_dot_tube()`).
        geometry = ContinuousFinStaggeredTubeBank(
            W_fro, H_fro, N_rows, S_trv, S_lon,
            D_ext, D_int, t_fin, N_fin, k_fin
        )
        A_min_tot = geometry.internal.A_min
        # Create the geometry of the desuperheating region, the condensing
        # region, and the subcooling region:
        self.desuperheating_region = DesuperheatingRegion(
            W_fro, H_fro, S_trv, S_lon, D_int, D_ext,
            t_fin, N_fin, k_fin,
            A_min_tot, N_rows
        )
        self.condensing_region = CondensingRegion(
            W_fro, H_fro, S_trv, S_lon, D_int, D_ext,
            t_fin, N_fin, k_fin,
            A_min_tot, N_rows
        )
        self.subcooling_region = SubcoolingRegion(
            W_fro, H_fro, S_trv, S_lon, D_int, D_ext,
            t_fin, N_fin, k_fin,
            A_min_tot, N_rows
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
            f"Condenser/Iteration {i + 1}: "
            f"Try with subcooling flow length {L_flow_sub:~P.3f}."
        )

        # Solve the subcooling region with the given subcooling flow length for
        # the state of air leaving the subcooling region:
        air_out = self.subcooling_region.solve(L_flow_sub)

        logger.debug(
            f"Condenser(SCR)/Iteration {i + 1}: "
            f"Leaving air temperature = {air_out.Tdb.to('degC'):~P.3f}. "
            "Determine condensing flow length to reject "
            f"{self.condensing_region.Q_dot.to('kW'):~P.3f}..."
        )

        # Set the state of air entering the condensing region to the state of
        # air leaving the subcooling region, and solve the condensing region
        # for its condensing flow length:
        self.condensing_region.air_in = air_out
        L_flow_cnd = self.condensing_region.solve(self.L_flow)

        logger.debug(
            f"Condenser(CDR)/Iteration {i + 1}: "
            f"Flow length = {L_flow_cnd.to('mm'):~P.3f}. "
            "Determine desuperheating flow length to reject "
            f"{self.desuperheating_region.Q_dot.to('kW'):~P.3f}"
        )

        # Set the state of air entering the desuperheating region to the state
        # of air leaving the condensing region, and solve the desuperheating
        # region for its desuperheating flow length:
        self.desuperheating_region.air_in = self.condensing_region.air_out
        L_flow_dsh = self.desuperheating_region.solve(self.L_flow)

        logger.debug(
            f"Condenser(DSR)/Iteration {i + 1}: "
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
            The smallest subcooling flow length to try with. By default, 0.1 mm.
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
                rtol=r_tol,  # %
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
