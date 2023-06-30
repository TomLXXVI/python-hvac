"""Rating of a plain finned tube cross flow heat exchanger with single-phase
fluids.

References
----------
[1] Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger Design.
    John Wiley & Sons. (§ 9.2.1 and § 9.3).
"""
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.forced_convection.external_flow.general import prandtl_number
from hvac.heat_transfer.finned_surface.fins import Fin, PlainContinuousFin
from hvac.heat_transfer.heat_exchanger.eps_ntu import CounterFlowHeatExchanger
from hvac.heat_transfer.heat_exchanger.misc import correct_friction_factor
from ..core.geometry import TubeBankInside, PlainFinStaggeredTBO
from ..correlations.plain_flat_fin.heat_transfer import j_Wang_and_Chi
from ..correlations.plain_flat_fin.friction_factor import f_Wang_and_Chi
from .circular_fin_tube_heat_exchanger import (
    CircularFinTubeCrossFlowHeatExchanger as CFT_CRF_HEX,
    TTubeBankOutside
)

Q_ = Quantity


class PlainFinTubeCrossFlowHeatExchanger(CFT_CRF_HEX):
    """Represent a single-pass, staggered, fin-tube heat exchanger with
    triangular equilateral pitch and plain flat fins in a crossflow arrangement.

    Parameters
    ----------
    L1:
        Tube length in the direction of inside flow, i.e. the length of the
        tube available for heat transfer with the outside flow.
    L2:
        Length of the tube bank parallel to the direction of external flow.
    L3:
        Length of the tube bank perpendicular to the direction of external
        flow.
    S_t:
        Lateral or transverse pitch, i.e. distance between tubes of the
        same row.
    S_l:
        Longitudinal pitch, i.e. distance between tubes of two adjacent tube
        rows.
    t_header:
        Header thickness; L1 + 2 * t_header is the tube length that needs
        to be taken into account for determining pressure drop along one
        tube.
    D_i:
        Inside diameter of the tubes.
    D_o:
        Outside diameter of the tubes.
    D_r:
        Effective diameter of the root of a circular fin. Depending upon the
        manufacturing techniques, it may be the tube outside diameter or tube
        outside diameter plus two times the fin collar thickness.
    t_f:
        Thickness of a circular fin.
    N_f:
        Number of fins per unit length (= fin density)
    k_f:
        Thermal conductivity of the external fin material.
    """
    def __init__(
        self,
        L1: Quantity,
        L2: Quantity,
        L3: Quantity,
        S_t: Quantity,
        S_l: Quantity,
        D_i: Quantity,
        D_o: Quantity,
        D_r: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        k_f: Quantity = Q_(237, 'W / (m * K)'),
        t_header: Quantity = Q_(0.0, 'mm')
    ) -> None:
        super().__init__(
            L1, L2, L3, S_t, S_l, D_i, D_o,
            D_r, None, t_f, N_f, k_f, t_header
        )

    @staticmethod
    def _create_internal_geometry(
        S_t: Quantity,
        S_l: Quantity,
        D_i: Quantity,
        L1: Quantity,
        L2: Quantity,
        L3: Quantity,
        t_header: Quantity = Q_(0.0, 'mm')
    ) -> TubeBankInside:
        geometry_int = TubeBankInside(
            S_t=S_t, S_l=S_l, D_i=D_i,
            L1=L1, L2=L2, L3=L3,
            t_header=t_header.to('m')
        )
        return geometry_int

    @staticmethod
    def _create_external_geometry(
        S_t: Quantity,
        S_l: Quantity,
        D_o: Quantity,
        D_r: Quantity,
        D_f: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        L1: Quantity,
        L2: Quantity,
        L3: Quantity
    ) -> TTubeBankOutside:
        geometry_ext = PlainFinStaggeredTBO(
            S_t=S_t, S_l=S_l, D_o=D_o, D_r=D_r,
            t_f=t_f, N_f=N_f, L1=L1, L2=L2, L3=L3
        )
        return geometry_ext

    def _get_ext_heat_trf_coeff(self, fluid_ext: FluidState) -> Quantity:
        Re_D_r_ext = self._get_ext_reynolds_number(
            D=self.geometry_ext.D_o,
            fluid_ext=fluid_ext
        )
        Pr_ext = prandtl_number(
            rho=fluid_ext.rho,
            mu=fluid_ext.mu,
            k=fluid_ext.k,
            cp=fluid_ext.cp
        )
        j_ext = j_Wang_and_Chi(
            Re_D_r=Re_D_r_ext,
            D_r=self.geometry_ext.D_r.to('m').m,
            D_h=self.geometry_ext.D_h.to('m').m,
            S_t=self.geometry_ext.S_t.to('m').m,
            S_l=self.geometry_ext.S_l.to('m').m,
            N_f=self.geometry_ext.N_f.to('1 / m').m,
            N_r=self.geometry_ext.N_r
        )
        Re_D_h_ext = self._get_ext_reynolds_number(
            D=self.geometry_ext.D_h,
            fluid_ext=fluid_ext
        )
        Nu_ext = j_ext * Re_D_h_ext * Pr_ext ** (1 / 3)
        h_ext = Nu_ext * fluid_ext.k / self.geometry_ext.D_h
        return h_ext.to('W / (m ** 2 * K)')

    def _get_fin(self) -> Fin:
        M = L = self.geometry_ext.S_t / 2
        fin = PlainContinuousFin(
            r_i=self.geometry_ext.D_r / 2,
            t=self.geometry_ext.t_f,
            k=self.k_f,
            M=M,
            L=L
        )
        return fin

    def _get_fanning_friction_factor(
        self,
        fluid_int: FluidState,
        h_int: Quantity,
        fluid_ext: FluidState,
        h_ext: Quantity,
        D_h: Quantity | None = None
    ) -> float:
        # calculate Fanning friction factor
        f = f_Wang_and_Chi(
            Re_D_r=self._get_ext_reynolds_number(
                self.geometry_ext.D_o,
                fluid_ext
            ),
            D_r=self.geometry_ext.D_r.to('m').m,
            S_t=self.geometry_ext.S_t.to('m').m,
            S_l=self.geometry_ext.S_l.to('m').m,
            N_f=self.geometry_ext.N_f.to('1 / m').m,
            N_r=self.geometry_ext.N_r
        )
        # correct friction factor for variable fluid properties effects
        f_corr = correct_friction_factor(
            f_cp=f,
            T_w=self._get_wall_temperature(fluid_int, h_int, fluid_ext, h_ext),
            flow_regime=self._get_ext_flow_regime(fluid_ext),
            thermal_regime=self._get_ext_thermal_regime(),
            fluid=fluid_ext
        )
        return f_corr


class PlainFinTubeCounterFlowHeatExchanger(PlainFinTubeCrossFlowHeatExchanger):

    def _get_hex_performance(
        self,
        fluid_int: FluidState,
        fluid_ext: FluidState,
        UA: Quantity
    ) -> tuple:
        # get outlet temperatures on internal and external side, heat transfer
        # rate, effectiveness and overall conductance of heat exchanger
        if self._get_fluid('hot') is self.Fluid_int:
            # internal fluid is hot fluid
            counterflow_hex = CounterFlowHeatExchanger(
                C_cold=self.m_dot_ext * fluid_ext.cp,
                C_hot=self.m_dot_int * fluid_int.cp,
                T_cold_in=self.T_ext_in,
                T_hot_in=self.T_int_in,
                UA=UA
            )
            T_int_out = counterflow_hex.T_hot_out
            T_ext_out = counterflow_hex.T_cold_out
        else:
            # internal flow is cold fluid
            counterflow_hex = CounterFlowHeatExchanger(
                C_cold=self.m_dot_int * fluid_int.cp,
                C_hot=self.m_dot_ext * fluid_ext.cp,
                T_cold_in=self.T_int_in,
                T_hot_in=self.T_ext_in,
                UA=UA
            )
            T_int_out = counterflow_hex.T_cold_out
            T_ext_out = counterflow_hex.T_hot_out
        return (
            T_int_out,
            T_ext_out,
            counterflow_hex.Q,
            counterflow_hex.eps,
            counterflow_hex.UA
        )

    def _get_int_flow_length(self) -> Quantity:
        return self.geometry_int.L2
