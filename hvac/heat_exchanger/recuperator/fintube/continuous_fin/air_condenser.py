import math
from abc import ABC, abstractmethod

import numpy as np
from scipy import optimize
from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import HumidAir, FluidState, CoolPropError, CP_HUMID_AIR
import hvac.heat_transfer.forced_convection.internal_flow as single_phase_flow
import hvac.heat_transfer.condensation.flow_condensation as condensing_flow
from hvac.heat_exchanger.recuperator.general.eps_ntu import CounterFlowHeatExchanger
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
        num_circuits: int | None = None
    ) -> None:
        self.geometry = ContinuousFinStaggeredTubeBank(
            width, height, num_rows, pitch_trv,
            pitch_lon, d_o, d_i, t_fin, fin_density,
            k_fin, d_r
        )
        self.num_circuits = num_circuits
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
        self.tube = None

    def _get_m_dot_tube(self) -> Quantity:
        if self.parent.num_circuits is None:
            n = self.parent.geometry.num_tubes_1st_row
        else:
            n = self.parent.num_circuits
        m_dot_tube = self.m_dot / n
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
        self.tube = single_phase_flow.CircularTube(
            Di=self.parent.geometry.d_i,
            L=self.parent.geometry.width,
            fluid=self.rfg
        )
        self.tube.m_dot = self._get_m_dot_tube()
        h = self.tube.avg_heat_transfer_coefficient()
        if isinstance(h, tuple):
            h = h[0]  # if laminar flow, assume constant heat flux
        return h.to('W / (m ** 2 * K)')


class CondensingPhaseInternalSurface(InternalSurface):

    def heat_transfer_coeff(self, T_wall: Quantity) -> Quantity:
        self.tube = condensing_flow.HorizontalTube(
            D=self.parent.geometry.internal.d_h,
            fluid=self.rfg.fluid
        )
        self.tube.m_dot = self._get_m_dot_tube()
        h = self.tube.heat_trf_coeff(
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
        geometry: ContinuousFinStaggeredTubeBank,
        num_circuits: int | None = None
    ) -> None:
        self.core = HeatExchangerCore(
            width=geometry.width,
            height=geometry.height,
            num_rows=0,
            pitch_trv=geometry.pitch_trv,
            pitch_lon=geometry.pitch_lon,
            d_i=geometry.d_i,
            d_o=geometry.d_o,
            t_fin=geometry.t_fin,
            fin_density=geometry.fin_density,
            k_fin=geometry.k_fin,
            num_circuits=num_circuits
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
        Calculates the heat transfer rate through the heat exchanger core of the
        desuperheating region and returns the deviation between this value and 
        the heat rate that the refrigerant must reject to turn from superheated 
        vapor into saturated vapor.
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
        Q_dot = cnt_flw_hex.Q
        dev = (Q_dot - self.Q_dot).to('W')
        i = counter[0]
        logger.debug(
            f"{i + 1}. Desuperheating region. "
            f"Flow length = {L_flow:~P.3f} -> "
            f"heat transfer rate = {Q_dot.to('W'):~P.3f} "
            f"(deviation = {dev:~P.3f})."
        )
        counter[0] += 1
        return dev.m

    def solve(
        self,
        L_flow_max: Quantity,
        x_tol: float = 0.001,
        r_tol: float = 0.01,
        i_max: int = 20
    ) -> Quantity:
        """
        Returns the required flow length of the desuperheating region so that 
        the heat transfer rate through the heat exchanger core would balance 
        with the heat rate the refrigerant must reject to turn into saturated
        vapor.
        """
        # Determine the state of air leaving the desuperheating region.
        h_air_out = self.air_in.h + self.Q_dot / self.air_m_dot
        self.air_out = HumidAir(h=h_air_out, W=self.air_in.W)
        # Determine the mean states of air and refrigerant in the desuperheating
        # region.
        air_mean, rfg_mean = self._get_fluid_mean_states()
        # Initialize iteration counter and set maximum number of iterations.
        counter = [0]
        # Within the given bracket, find the flow length of the desuperheating
        # region for which function `__fun__` returns a deviation near zero.
        try:
            sol = optimize.root_scalar(
                self.__fun__,
                args=(air_mean, rfg_mean, counter),
                method='brentq',
                bracket=(1e-1, L_flow_max.to('mm').m),
                xtol=x_tol,
                rtol=r_tol,
                maxiter=i_max
            )
        except ValueError:
            raise DesuperheatingError(
                "Refrigerant cannot be desuperheated to saturated vapor "
                "under current operating conditions."
            ) from None
        self.L_flow = Q_(sol.root, 'mm')
        self.core.geometry.length = self.L_flow
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
        geometry: ContinuousFinStaggeredTubeBank,
        num_circuits: int | None = None
    ) -> None:
        self.core = HeatExchangerCore(
            width=geometry.width,
            height=geometry.height,
            num_rows=0,
            pitch_trv=geometry.pitch_trv,
            pitch_lon=geometry.pitch_lon,
            d_i=geometry.d_i,
            d_o=geometry.d_o,
            t_fin=geometry.t_fin,
            fin_density=geometry.fin_density,
            k_fin=geometry.k_fin,
            condensing=True,
            num_circuits=num_circuits
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
        Q_dot = cnt_flw_hex.Q
        dev = (Q_dot - self.Q_dot).to('W')
        i = counter[0]
        logger.debug(
            f"{i + 1}. Condensing region. "
            f"Flow length = {L_flow:~P.3f} -> "
            f"heat transfer rate = {Q_dot.to('W'):~P.3f} "
            f"(deviation = {dev:~P.3f})."
        )
        counter[0] += 1
        return dev.m

    def solve(
        self,
        L_flow_max: Quantity,
        x_tol: float = 0.001,
        r_tol: float = 0.01,
        i_max: int = 20
    ) -> Quantity:
        """
        Returns the required flow length of the condensing region so that the 
        heat transfer rate through the heat exchanger core would balance with 
        the heat rate that the refrigerant must reject in the condensing region 
        to turn from saturated vapor into saturated liquid.
        """
        # Determine the state of air leaving the condensing region.
        h_air_out = self.air_in.h + self.Q_dot / self.air_m_dot
        self.air_out = HumidAir(h=h_air_out, W=self.air_in.W)
        # Determine the mean states of air and refrigerant in the condensing
        # region.
        air_mean, rfg_mean = self._get_fluid_mean_states()
        # Initialize iteration counter.
        counter = [0]
        # Within the given bracket, find the flow length of the condensing
        # region for which function `__fun__` returns a deviation near zero.
        try:
            sol = optimize.root_scalar(
                self.__fun__,
                args=(air_mean, rfg_mean, counter),
                method='brentq',
                bracket=(1e-1, L_flow_max.to('mm').m),
                xtol=x_tol,
                rtol=r_tol,
                maxiter=i_max
            )
        except ValueError:
            raise CondensingError(
                "Refrigerant cannot be condensed to saturated liquid "
                "under current operating conditions."
            ) from None
        self.L_flow = Q_(sol.root, 'mm')
        self.core.geometry.length = self.L_flow
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
        geometry: ContinuousFinStaggeredTubeBank,
        num_circuits: int | None = None
    ) -> None:
        self.core = HeatExchangerCore(
            width=geometry.width,
            height=geometry.height,
            num_rows=0,
            pitch_trv=geometry.pitch_trv,
            pitch_lon=geometry.pitch_lon,
            d_i=geometry.d_i,
            d_o=geometry.d_o,
            t_fin=geometry.t_fin,
            fin_density=geometry.fin_density,
            k_fin=geometry.k_fin,
            num_circuits=num_circuits
        )
        # Known parameters:
        self.rfg_in: FluidState | None = None  # saturated liquid
        self.rfg_m_dot: Quantity | None = None
        self.air_in: HumidAir | None = None
        self.air_m_dot: Quantity | None = None

        self.dT_sc: Quantity | None = None

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

    def __fun_analysis__(
        self,
        T_rfg_out: Quantity,
        L_flow: Quantity,
        counter: list[int]
    ) -> tuple[Quantity, Quantity]:
        """
        With the given flow length `L_flow` of the subcooling region, calculates
        a new value for the temperature of the refrigerant leaving the 
        condenser.
        """
        # Current state of the refrigerant leaving the condenser.
        rfg_out = self._get_rfg_out(T_rfg_out)
        # Heat rate rejected by the refrigerant in the subcooling region for
        # the current state of the leaving refrigerant.
        Q_dot = self.rfg_m_dot * (self.rfg_in.h - rfg_out.h)
        # Determine the mean states of air and refrigerant in the subcooling
        # region.
        T_air_out = self.air_in.Tdb + Q_dot / (CP_HUMID_AIR * self.air_m_dot)
        air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
        air_mean, rfg_mean = self._get_fluid_mean_states(air_out, rfg_out)
        # Set heat exchanger core parameters to calculate UA.
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
        # exchanger core of the subcooling region.
        cnt_flw_hex = CounterFlowHeatExchanger(
            C_cold=self.air_m_dot * air_mean.cp,
            C_hot=self.rfg_m_dot * rfg_mean.cp,
            T_cold_in=self.air_in.Tdb,
            T_hot_in=self.rfg_in.T,
            UA=self.core.UA
        )
        Q_dot_new = cnt_flw_hex.Q
        # Determine the new state of refrigerant leaving the condenser.
        h_rfg_out_new = self.rfg_in.h - Q_dot_new / self.rfg_m_dot
        Rfg = self.rfg_in.fluid
        P_cnd = self.rfg_in.P
        rfg_out_new = Rfg(P=P_cnd, h=h_rfg_out_new)
        # Temperature deviation between the new state and the current state of
        # the refrigerant leaving the condenser.
        dev = (rfg_out_new.T - rfg_out.T).to('K')
        i = counter[0]
        logger.debug(
            f"{i + 1}. Subcooling region. "
            f"Temperature of leaving refrigerant = {rfg_out_new.T.to('degC'):~P.3f} "
            f"(deviation = {dev:~P.3f})."
        )
        counter[0] += 1
        return dev, rfg_out_new.T

    def solve_analysis(
        self,
        L_flow: Quantity,
        tol: Quantity = Q_(0.1, 'K'),
        i_max: int = 20
    ) -> HumidAir:
        """
        With the given subcooling flow length `L_flow`, solve for the output
        parameters of the subcooling region `air_out`, `rfg_out`, and `Q_dot`.
        """
        T_rfg_out_min = self.air_in.Tdb.to('degC')
        T_rfg_out_max = self.rfg_in.T.to('degC')
        counter = [0]
        # Initial guess of the temperature of the refrigerant leaving the
        # condenser.
        T_rfg_out = (T_rfg_out_min.to('K') + T_rfg_out_max.to('K')) / 2
        for i in range(i_max):
            # Determine the heat transfer rate in the subcooling region and from
            # this determine the corresponding new state of refrigerant leaving 
            # the condenser.
            # Repeat these calculations until the temperature of the leaving
            # refrigerant becomes nearly constant (or until the maximum number 
            # of loop iterations has been reached).
            dev, T_rfg_out_new = self.__fun_analysis__(T_rfg_out, L_flow, counter)
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
                f"State of refrigerant leaving the condenser "
                f"could not be determined after {i_max} iterations."
            ) from None

    def __fun_design__(
        self,
        L_flow: float,
        air_mean: HumidAir,
        rfg_mean: FluidState,
        counter: list[int]
    ) -> float:
        """
        Calculates the heat transfer rate through the heat exchanger core of the
        subcooling region and returns the deviation between this value and the
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
            C_hot=self.rfg_m_dot * rfg_mean.cp,
            T_cold_in=self.air_in.Tdb,
            T_hot_in=self.rfg_in.T,
            UA=self.core.UA
        )
        Q_dot = cnt_flw_hex.Q
        dev = (Q_dot - self.Q_dot).to('W')
        i = counter[0]
        logger.debug(
            f"{i + 1}. Subcooling region. "
            f"Flow length = {L_flow:~P.3f} -> "
            f"heat transfer rate = {Q_dot.to('W'):~P.3f} "
            f"(deviation = {dev:~P.3f})."
        )
        counter[0] += 1
        return dev.m

    def solve_design(
        self,
        L_flow_max: Quantity,
        x_tol: float = 0.001,
        r_tol: float = 0.01,
        i_max: int = 20
    ) -> Quantity:
        # Determine the heat rejection rate of the refrigerant in the subcooling
        # region.
        self.Q_dot = self.rfg_m_dot * (self.rfg_in.h - self.rfg_out.h)
        # Determine the state of air leaving the subcooling region.
        h_air_out = self.air_in.h + self.Q_dot / self.air_m_dot
        self.air_out = HumidAir(h=h_air_out, W=self.air_in.W)
        # Determine the mean states of air and refrigerant in the subcooling
        # region.
        air_mean, rfg_mean = self._get_fluid_mean_states(self.air_out, self.rfg_out)
        # Initialize iteration counter.
        counter = [0]
        # Within the given bracket, find the flow length of the subcooling
        # region for which function `__fun_design__` returns a deviation near
        # zero.
        try:
            sol = optimize.root_scalar(
                self.__fun_design__,
                args=(air_mean, rfg_mean, counter),
                method='brentq',
                bracket=(1e-1, L_flow_max.to('mm').m),
                xtol=x_tol,
                rtol=r_tol,
                maxiter=i_max
            )
        except ValueError:
            raise SubcoolingError(
                f"Refrigerant cannot be subcooled with {self.dT_sc.to('K'):~P.0f} "
                "under current operating conditions."
            ) from None
        self.L_flow = Q_(sol.root, 'mm')
        self.core.geometry.length = self.L_flow
        return self.L_flow

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
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
        num_circuits: int | None = None
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
        num_circuits: int | None, default None
            Number of refrigerant circuits. If None, the number of circuits
            is set equal to the number of tubes in the first row.
        """
        # Create the heat exchanger geometry of the condenser.
        self.geometry = ContinuousFinStaggeredTubeBank(
            W_fro, H_fro, N_rows, S_trv, S_lon,
            D_ext, D_int, t_fin, N_fin, k_fin
        )
        # Create the heat exchanger of the desuperheating region, the condensing
        # region, and the subcooling region:
        self.desuperheating_region = DesuperheatingRegion(
            self.geometry, num_circuits
        )
        self.condensing_region = CondensingRegion(
            self.geometry, num_circuits
        )
        self.subcooling_region = SubcoolingRegion(
            self.geometry, num_circuits
        )
        # Determine total flow length of the condenser.
        self.L_flow = N_rows * self.geometry.pitch_lon

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

    def solve(
        self,
        air_in: HumidAir,
        air_m_dot: Quantity,
        rfg_in: FluidState,
        rfg_m_dot: Quantity,
        dT_sc: Quantity | None = None
    ) -> tuple[HumidAir, FluidState] | tuple[Quantity, int]:
        """
        Function `solve(...)` can have two different applications:

        1.  Analysis:
            Under the given operating conditions, determines the operating state
            of the condenser of which the flow length/number of rows was
            specified at instantiation of the
            `PlainFinTubeCounterFlowAirCondenser`.

        2.  Design (sizing):
            Under the given operating conditions, determines the total flow
            length of the condenser needed so that refrigerant leaves the
            condenser with the required degree of subcooling.

        Parameters
        ----------
        air_in: HumidAir
            State of the air entering the condenser.
        air_m_dot: Quantity
            Mass flow rate of air through the condenser.
        rfg_in: FluidState
            State of the refrigerant entering the condenser.
        rfg_m_dot: Quantity
            Mass flow rate of refrigerant through the condenser.
        dT_sc: Quantity, optional
            Required degree of refrigerant subcooling. If specified, the design
            problem is solved. If `None`, the analysis problem is solved.

        Returns
        -------
        tuple[HumidAir, FluidState] | tuple[Quantity, int]
            If parameter `dT_sc` is None, the analysis problem is solved.
                air_out: HumidAir
                    State of air leaving the condenser in the desuperheating
                    region.
                rfg_out: FluidState
                    State of refrigerant leaving the condenser in the subcooled
                    region.
            If parameter `dT_sc` is specified, the design problem is solved.
                L_flow: Quantity
                    Total flow length of the condenser.
                N_rows: int
                    Number of rows of the condenser.
        """
        if dT_sc is None:
            return self.solve_analysis(
                air_in, air_m_dot,
                rfg_in, rfg_m_dot,
            )
        else:
            return self.solve_design(
                air_in, air_m_dot,
                rfg_in, rfg_m_dot,
                dT_sc
            )

    def __fun_analysis__(
        self,
        L_flow_sub: float,
        counter: list[int]
    ) -> float:
        """
        Calculates the deviation between the calculated total flow length of the
        condenser and its actual total flow length, starting with a given
        flow length for the subcooling region only.
        """
        L_flow_sub = Q_(L_flow_sub, 'mm')
        i = counter[0]
        logger.debug(
            f"Iteration {i + 1}"
        )
        logger.debug(        
            f"Try with subcooling flow length {L_flow_sub:~P.3f}."
        )
        # Solve the subcooling region with the given subcooling flow length for
        # the state of air leaving the subcooling region:
        air_out = self.subcooling_region.solve_analysis(L_flow_sub)
        logger.debug(
            f"Subcooling region. "
            f"Leaving air temperature = {air_out.Tdb.to('degC'):~P.3f}."
        )
        # Set the state of air entering the condensing region to the state of
        # air leaving the subcooling region, and solve the condensing region
        # for its condensing flow length:
        self.condensing_region.air_in = air_out
        L_flow_cnd = self.condensing_region.solve(self.L_flow)
        logger.debug(
            f"Condenser region. "
            f"Required flow length = {L_flow_cnd.to('mm'):~P.3f}."
        )
        # Set the state of air entering the desuperheating region to the state
        # of air leaving the condensing region, and solve the desuperheating
        # region for its desuperheating flow length:
        self.desuperheating_region.air_in = self.condensing_region.air_out
        L_flow_dsh = self.desuperheating_region.solve(self.L_flow)
        logger.debug(
            f"Desuperheating region. "
            f"Required flow length = {L_flow_dsh.to('mm'):~P.3f}."
        )
        # Determine the deviation between the currently calculated total flow
        # length of the condenser and its actual total flow length:
        L_flow_new = L_flow_sub + L_flow_cnd + L_flow_dsh
        dev = (L_flow_new - self.L_flow).to('mm')
        logger.debug(
            f"Condenser. "
            f"Total flow length = {L_flow_new.to('mm'):~P.3f} "
            f"(deviation = {dev:~P.3f})."
        )
        counter[0] += 1
        return dev.m

    def solve_analysis(
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
        Under the given operating conditions, determines the operating state of
        the condenser of which the flow length/number of rows was specified at
        instantiation of the `PlainFinTubeCounterFlowAirCondenser`.

        Parameters
        ----------
        air_in: HumidAir
            State of the air entering the condenser.
        air_m_dot: Quantity
            Mass flow rate of air through the condenser.
        rfg_in: FluidState
            State of the refrigerant entering the condenser.
        rfg_m_dot: Quantity
            Mass flow rate of refrigerant through the condenser.
        L_flow_sub_min: Quantity, optional
            The smallest subcooling flow length to try with. By default, 0.1 mm.
        x_tol: float, optional
            Absolute tolerance for root-finding algorithm
            (`scipy.optimize.root_scalar` with 'brentq' method).
        r_tol: float, optional
            Relative tolerance for root-finding algorithm.
        i_max: int, optional
            Maximum number of iterations for root-finding algorithm
        
        Returns
        -------
        tuple[HumidAir, FluidState]
            The state of air and the state of refrigerant leaving the condenser.

        Raises
        ------
        SubcoolingError:
            If the state of refrigerant leaving the condensor cannot be
            determined.
        CondensingError:
            If the required flow length for condensing the refrigerant cannot
            be determined (i.e. if the heat transfer rate through the heat 
            exchanger core of the condensing region cannot be balanced with the 
            heat rate the refrigerant must reject to turn from saturated vapor 
            into saturated liquid).
        DesuperheatingError:
            If the required flow length for desuperheating the refrigerant 
            cannot be determined (i.e. if the heat transfer rate through the 
            heat exchanger core of the desuperheating region cannot be balanced 
            with the heat rate the refrigerant must reject to turn from the 
            inlet state into saturated vapor into saturated).
        CondenserError:
            If the required flow length of the condenser cannot be made equal
            to the actual flow length of the condenser.
        
        Notes
        -----
        The condenser has three regions. First, the refrigerant is desuperheated
        and becomes saturated vapor. Next, the refrigerant condenses into
        saturated liquid. Finally, the refrigerant is subcooled before it leaves
        the condenser.

        We know the mass flow rate of refrigerant and the mass flow rate
        of air. We also know the state of air entering the subcooled region and
        the state of refrigerant entering the desuperheating region.

        We can choose a subcooling flow length and solve for the state of
        air leaving the subcooling region and entering the condensing region.

        The refrigerant state at the entrance of the desuperheating region, the
        entrance of the condensing region (saturated vapor), and the entrance of
        the subcooling region (saturated liquid) are known from the start.
        Since we also know the refrigerant mass flow rate, we can determine the
        heat rate the refrigerant must reject in the condensing and in the
        desuperheating region. From this, we can determine the required flow
        length of the condensing region and of the desuperheating region.

        A solution is found when the sum of the chosen subcooling flow length
        and the calculated flow lengths of the condensing and the desuperheating
        region equals the actual total flow length of the condenser.
        """
        # Assign/calculate what is known or can be calculated directly.
        self._set_air_in(air_in)
        self._set_air_m_dot(air_m_dot)
        self._set_rfg_in(rfg_in)
        self._set_rfg_m_dot(rfg_m_dot)
        # Find the flow length of the subcooling region so that the sum of the
        # flow lengths of the subcooling, condensing, and desuperheating region
        # equals the actual, total flow length of the condenser.
        L_flow_sub_min = L_flow_sub_min.to('mm').m
        L_flow_sub_max = self.L_flow.to('mm').m
        counter = [0]
        try:
            sol = optimize.root_scalar(
                self.__fun_analysis__,
                args=(counter,),
                method='brentq',
                bracket=(L_flow_sub_min, L_flow_sub_max),
                xtol=x_tol,  # mm
                rtol=r_tol,  # %
                maxiter=i_max
            )
        except ValueError:
            message = (
                "The refrigerant cannot be subcooled under "
                "current operating conditions."
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

    def solve_design(
        self,
        air_in: HumidAir,
        air_m_dot: Quantity,
        rfg_in: FluidState,
        rfg_m_dot: Quantity,
        dT_sc: Quantity
    ) -> tuple[Quantity, int]:
        """
        Under the given operating conditions, determines the total flow length
        of the condenser so that refrigerant leaves the condenser with the
        desired degree of subcooling.

        Parameters
        ----------
        air_in: HumidAir
            State of the air entering the condenser.
        air_m_dot: Quantity
            Mass flow rate of air through the condenser.
        rfg_in: FluidState
            State of the refrigerant entering the condenser.
        rfg_m_dot: Quantity
            Mass flow rate of refrigerant through the condenser.
        dT_sc: Quantity
            Desired degree of refrigerant subcooling.

        Returns
        -------
        L_flow: Quantity
            Total flow length of the condenser.
        N_rows: int
            Number of rows of the condenser.
        """
        # Assign/calculate what is known or can be calculated directly:
        self._set_air_in(air_in)
        self._set_air_m_dot(air_m_dot)
        self._set_rfg_in(rfg_in, dT_sc)
        self._set_rfg_m_dot(rfg_m_dot)

        # Subcooling region: determine the required flow length to subcool
        # the refrigerant from saturated liquid to subcooled liquid having the
        # required degree of subcooling.
        L_flow_subcool = None
        for k in range(4):
            L_flow_max = (5 ** k) * self.L_flow
            try:
                L_flow_subcool = self.subcooling_region.solve_design(
                    L_flow_max=L_flow_max,
                )
            except SubcoolingError:
                logger.debug('Try again...')
                continue
            else:
                break
        if L_flow_subcool is None:
            raise SubcoolingError(
                "Required flow length of subcooling region "
                "could not be determined."
            ) from None
        logger.debug(
            f"Required subcooling flow length = "
            f"{L_flow_subcool.to('mm'):~P.3f}"
        )
        # Condensing region: determine the required flow length to completely
        # condense the refrigerant from saturated vapor to saturated liquid.
        L_flow_cond = None
        self.condensing_region.air_in = self.subcooling_region.air_out
        for k in range(4):
            L_flow_max = (5 ** k) * self.L_flow
            try:
                L_flow_cond = self.condensing_region.solve(
                    L_flow_max=L_flow_max
                )
            except CondensingError:
                logger.debug('Try again...')
                continue
            else:
                break
        if L_flow_cond is None:
            raise CondensingError(
                "Required flow length of condensing region "
                "could not be determined"
            ) from None
        logger.debug(
            f"Required condensing flow length = "
            f"{L_flow_cond.to('mm'):~P.3f}"
        )
        # Desuperheating region: determine the required flow length to
        # desuperheat the refrigerant from superheated vapor to saturated
        # liquid.
        L_flow_desuper = None
        self.desuperheating_region.air_in = self.condensing_region.air_out
        for k in range(4):
            L_flow_max = (5 ** k) * self.L_flow
            try:
                L_flow_desuper = self.desuperheating_region.solve(
                    L_flow_max=L_flow_max
                )
            except DesuperheatingError:
                logger.debug('Try again...')
                continue
            else:
                break
        if L_flow_desuper is None:
            raise DesuperheatingError(
                "Required flow length of desuperheating region "
                "could not be determined"
            ) from None
        logger.debug(
            f"Required desuperheating flow length = "
            f"{L_flow_desuper.to('mm'):~P.3f}"
        )
        L_flow = L_flow_subcool + L_flow_cond + L_flow_desuper
        N_rows = int(np.ceil(
            L_flow.to('mm').m / self.geometry.pitch_lon.to('mm').m
        ))
        self._set_flow_length(N_rows)

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
        return L_flow, N_rows

    def _set_air_in(self, air_in: HumidAir) -> None:
        self.air_in = air_in
        self.subcooling_region.air_in = air_in

    def _set_air_m_dot(self, air_m_dot: Quantity) -> None:
        self.air_m_dot = air_m_dot
        self.subcooling_region.air_m_dot = air_m_dot
        self.condensing_region.air_m_dot = air_m_dot
        self.desuperheating_region.air_m_dot = air_m_dot

    def _set_rfg_in(
        self,
        rfg_in: FluidState,
        dT_sc: Quantity | None = None
    ) -> None:
        self.rfg_in = rfg_in
        self.desuperheating_region.rfg_in = rfg_in
        self.P_cnd = self.rfg_in.P
        Rfg = rfg_in.fluid
        self.desuperheating_region.rfg_out = Rfg(P=self.P_cnd, x=Q_(1, 'frac'))
        self.condensing_region.rfg_in = self.desuperheating_region.rfg_out
        self.condensing_region.rfg_out = Rfg(P=self.P_cnd, x=Q_(0, 'frac'))
        self.T_cnd = self.condensing_region.rfg_out.T
        self.subcooling_region.rfg_in = self.condensing_region.rfg_out
        if isinstance(dT_sc, Quantity):
            self.dT_sc = dT_sc
            self.subcooling_region.dT_sc = dT_sc
            self.subcooling_region.rfg_out = Rfg(
                P=self.P_cnd,
                T=self.T_cnd - dT_sc.to('K')
            )
            self.rfg_out = self.subcooling_region.rfg_out

    def _set_rfg_m_dot(self, rfg_m_dot: Quantity) -> None:
        self.rfg_m_dot = rfg_m_dot
        self.desuperheating_region.rfg_m_dot = rfg_m_dot
        self.condensing_region.rfg_m_dot = rfg_m_dot
        self.subcooling_region.rfg_m_dot = rfg_m_dot

    def _set_flow_length(self, N_rows: int) -> None:
        self.L_flow = N_rows * self.geometry.pitch_lon
        self.geometry.length = self.L_flow

    @property
    def N_rows(self) -> int:
        return self.geometry.num_rows
