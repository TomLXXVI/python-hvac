"""
Sizing of a plain flat finned tube counterflow or crossflow heat exchanger
for single-phase fluids.

References
----------
[1] Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger Design.
    John Wiley & Sons. (§ 9.2.2 and § 9.3).
"""
import numpy as np
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.forced_convection.external_flow.general import prandtl_number
from hvac.heat_transfer.finned_surface.fins import PlainContinuousFin
from hvac.heat_transfer.heat_exchanger import eps_ntu
from ..core.geometry import PlainFinStaggeredTBO
from ..correlations.plain_flat_fin.heat_transfer import j_Wang_and_Chi
from ..correlations.plain_flat_fin.friction_factor import f_Wang_and_Chi
from .circular_fin_tube_heat_exchanger import (
    CircularFinTubeCounterFlowHeatExchanger as CFT_COF_HEX,
    TTubeBankOutside
)

Q_ = Quantity


class PlainFinTubeCounterFlowHeatExchanger(CFT_COF_HEX):
    """Represent a single-pass, staggered, fin-tube heat exchanger with
    equilateral triangular pitch and with plain flat fins in a counterflow
    arrangement.
    """
    def __init__(
        self,
        D_i: Quantity,
        D_o: Quantity,
        S_t: Quantity,
        S_l: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        alpha_int: Quantity,
        alpha_ext: Quantity,
        sigma_ext: float,
        L1: Quantity | None = None,
        L3: Quantity | None = None,
        k_f: Quantity = Q_(237, 'W / (m * K)')  # aluminium
    ) -> None:
        """Creates a `CircularFinTubeCounterFlowHeatExchanger` instance.

        Parameters
        ----------
        D_i:
            Inner diameter of tubes.
        D_o:
            Outer diameter of tubes.
        S_t:
            Lateral or transverse pitch, i.e. the spacing between tubes of the
            same row.
        S_l:
            Longitudinal pitch, i.e. the spacing between tubes of two adjacent
            tube row.
        t_f:
            Thickness of a circular fin.
        N_f:
            Number of fins per unit length (= fin density)
        alpha_int:
            Ratio of internal heat transfer surface area to total volume of the
            heat exchanger.
        alpha_ext:
            Ratio of external heat transfer surface area to total volume of the
            heat exchanger.
        sigma_ext:
            Ratio of the free flow area to the frontal area of the heat
            exchanger.
        L1: optional
            Tube length in the direction of internal flow, i.e. the core length
            of the tube available for heat transfer with the external flow or
            the width of the heat exchanger core. If `L1` is left to None,
            parameter `L3` must be specified.
        L3: optional
            Length of the heat exchanger core perpendicular to the direction of
            external flow, i.e. the height of the heat exchanger core. If `L3`
            is left to None, parameter `L1` must be specified.
        k_f: default 237 W/(m.K) for aluminium
            Thermal conductivity of the external fin material.
        """
        super().__init__(
            D_i, D_o, S_t, S_l, None, t_f, N_f, alpha_int,
            alpha_ext, sigma_ext, L1, L3, k_f
        )

    @staticmethod
    def _create_external_geometry(
        S_t: Quantity,
        S_l: Quantity,
        D_o: Quantity,
        D_f: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        alpha_ext: Quantity | None,
        sigma_ext: float | None,
        L1: Quantity | None,
        L2: Quantity | None,
        L3: Quantity | None
    ) -> TTubeBankOutside:
        geometry_ext = PlainFinStaggeredTBO(
            S_t=S_t.to('m'),
            S_l=S_l.to('m'),
            D_o=D_o.to('m'),
            D_r=D_o.to('m'),
            t_f=t_f.to('m'),
            N_f=N_f.to('1 / m'),
            alpha=alpha_ext.to('1 / m') if alpha_ext is not None else None,
            sigma=sigma_ext if sigma_ext is not None else None,
            L1=L1.to('m') if L1 is not None else None,
            L2=L2.to('m') if L2 is not None else None,
            L3=L3.to('m') if L3 is not None else None
        )
        return geometry_ext

    def _get_jf_ratio_ext(
        self,
        Re_D_o_ext: float | None,
        N_r: float | None = None
    ) -> float:
        def __colburn_j_factor(Re: float) -> float:
            return j_Wang_and_Chi(
                Re_D_r=Re,
                D_r=self.geometry_ext.D_r.to('m').m,
                D_h=self.geometry_ext.D_h.to('m').m,
                S_t=self.geometry_ext.S_t.to('m').m,
                S_l=self.geometry_ext.S_l.to('m').m,
                N_f=self.geometry_ext.N_f.to('1 / m').m,
                N_r=N_r
            )

        def __friction_factor(Re: float) -> float:
            f = f_Wang_and_Chi(
                Re_D_r=Re,
                D_r=self.geometry_ext.D_r.to('m').m,
                S_t=self.geometry_ext.S_t.to('m').m,
                S_l=self.geometry_ext.S_l.to('m').m,
                N_f=self.geometry_ext.N_f.to('1 / m').m,
                N_r=N_r
            )
            return f

        if Re_D_o_ext is None:
            Re_rng = np.arange(100, 50_000, 100)
            j_rng = np.array([__colburn_j_factor(Re) for Re in Re_rng])
            f_rng = np.array([__friction_factor(Re) for Re in Re_rng])
            jf_ratio_rng = j_rng / f_rng
            jf_ratio = float(np.mean(jf_ratio_rng))
        else:
            j = __colburn_j_factor(Re_D_o_ext)
            f = __friction_factor(Re_D_o_ext)
            jf_ratio = j / f
        return jf_ratio

    def _get_ext_heat_trf_coeff(
        self,
        fluid_ext: FluidState,
        G_ext: Quantity,
        N_r: float | None = None
    ) -> Quantity:
        Re_ext = self._get_reynolds_number_ext(
            fluid_ext, G_ext,
            self.geometry_ext.D_o
        )
        Pr_ext = prandtl_number(
            rho=fluid_ext.rho,
            mu=fluid_ext.mu,
            k=fluid_ext.k,
            cp=fluid_ext.cp
        )
        j_ext = j_Wang_and_Chi(
            Re_D_r=Re_ext,
            D_r=self.geometry_ext.D_r.to('m').m,
            D_h=self.geometry_ext.D_h.to('m').m,
            S_t=self.geometry_ext.S_t.to('m').m,
            S_l=self.geometry_ext.S_l.to('m').m,
            N_f=self.geometry_ext.N_f.to('1 / m').m,
            N_r=N_r
        )
        D_h_ext = self.geometry_ext.D_h
        Re_D_h_ext = self._get_reynolds_number_ext(fluid_ext, G_ext, D_h_ext)
        Nu_ext = j_ext * Re_D_h_ext * Pr_ext ** (1 / 3)
        h_ext = Nu_ext * fluid_ext.k / D_h_ext
        return h_ext.to('W / (m ** 2 * K)')

    def _get_ext_friction_factor(
        self,
        fluid_ext: FluidState,
        G_ext: Quantity,
        N_r: float | None = None
    ) -> float:
        Re_ext_D_r = self._get_reynolds_number_ext(
            fluid_ext, G_ext,
            self.geometry_ext.D_r
        )
        f = f_Wang_and_Chi(
            Re_D_r=Re_ext_D_r,
            D_r=self.geometry_ext.D_r.to('m').m,
            S_t=self.geometry_ext.S_t.to('m').m,
            S_l=self.geometry_ext.S_l.to('m').m,
            N_f=self.geometry_ext.N_f.to('1 / m').m,
            N_r=N_r
        )
        return f

    def _get_overall_fin_efficiency(self, h_ext) -> float:
        # get overall efficiency of finned surface
        M = L = self.geometry_ext.S_t / 2
        fin = PlainContinuousFin(
            r_i=self.geometry_ext.D_r / 2,
            t=self.geometry_ext.t_f,
            k=self.k_f,
            M=M,
            L=L
        )
        fin.h_avg = h_ext
        eta_fin = fin.efficiency
        r = self.geometry_ext.A_f_to_A
        eta_ext = 1 - r * (1 - eta_fin)
        return eta_ext


class PlainFinTubeCrossFlowHeatExchanger(PlainFinTubeCounterFlowHeatExchanger):
    """Represent a single-pass, staggered, fin-tube heat exchanger with
    equilateral triangular pitch and with circular fins in a crossflow
    arrangement.
    """
    def __init__(
        self,
        D_i: Quantity,
        D_o: Quantity,
        S_t: Quantity,
        S_l: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        alpha_int: Quantity,
        alpha_ext: Quantity,
        sigma_ext: float,
        k_f: Quantity = Q_(237, 'W / (m * K)')  # aluminium
    ) -> None:
        """Creates a `CircularFinTubeCrossFlowHeatExchanger` instance.

        Parameters
        ----------
        D_i:
            Inner diameter of tubes.
        D_o:
            Outer diameter of tubes.
        S_t:
            Lateral or transverse pitch, i.e. the spacing between tubes of the
            same row.
        S_l:
            Longitudinal pitch, i.e. the spacing between tubes of two adjacent
            tube row.
        t_f:
            Thickness of a circular fin.
        N_f:
            Number of fins per unit length (= fin density)
        alpha_int:
            Ratio of internal heat transfer surface area to total volume of the
            heat exchanger.
        alpha_ext:
            Ratio of external heat transfer surface area to total volume of the
            heat exchanger.
        sigma_ext:
            Ratio of the free flow area to the frontal area of the heat
            exchanger.
        k_f: default 237 W/(m.K) for aluminium
            Thermal conductivity of the external fin material.
        """
        super().__init__(
            D_i=D_i,
            D_o=D_o,
            S_t=S_t,
            S_l=S_l,
            t_f=t_f,
            N_f=N_f,
            alpha_int=alpha_int,
            alpha_ext=alpha_ext,
            sigma_ext=sigma_ext,
            k_f=k_f
        )

    @staticmethod
    def _get_int_flow_length(L1: Quantity, L2: Quantity) -> Quantity:
        return L1

    def _get_heat_transfer_params(
        self,
        fluid_int: FluidState,
        fluid_ext: FluidState,
    ) -> eps_ntu.CrossFlowHeatExchanger:
        if self._get_fluid('hot') is self.Fluid_int:
            # internal fluid is hot fluid
            C_cold = fluid_ext.cp * self.m_dot_ext
            T_cold_in = self.T_ext_in
            C_hot = fluid_int.cp * self.m_dot_int
            T_hot_in = self.T_int_in
        else:
            # internal fluid is cold fluid
            C_cold = fluid_int.cp * self.m_dot_int
            T_cold_in = self.T_int_in
            C_hot = fluid_ext.cp * self.m_dot_ext
            T_hot_in = self.T_ext_in
        if (self.Q is None) and (self.eps is not None):
            self.Q = self._get_heat_transfer_rate(
                eps=self.eps,
                C_min=min(C_cold, C_hot)
            )
        if (self.Q is None) and (self.eps is None):
            raise ValueError('parameter `Q` and `eps` cannot be both None')
        crf_hex = eps_ntu.CrossFlowHeatExchanger(
            type_=eps_ntu.CrossFlowHeatExchanger.Type.C_min_MIXED_C_max_UNMIXED,
            C_cold=C_cold,
            C_hot=C_hot,
            T_cold_in=T_cold_in,
            T_hot_in=T_hot_in,
            Q=self.Q
        )
        return crf_hex

    def _get_core_dimensions(
        self,
        NTU: float,
        C_min: Quantity,
        U_ext: Quantity,
        G_ext: Quantity,
        G_int: Quantity | None = None
    ) -> tuple:
        # determine external heat transfer area
        A_ext = self._get_heat_transfer_area_ext(NTU, C_min, U_ext)
        # determine internal heat transfer area
        A_int = self._get_heat_transfer_area_int(A_ext)
        # determine external free flow area
        A_o_ext = self._get_free_flow_area_ext(G_ext)
        # determine internal free flow area
        A_o_int = self._get_free_flow_area_int(G_int)
        # determine external frontal area
        A_fr_ext = self._get_frontal_area_ext(A_o_ext)
        # determine internal frontal area
        A_fr_int = self._get_frontal_area_int(A_o_int)
        # determine heat exchanger dimensions
        L2 = self._get_core_length_ext(A_ext, A_o_ext)
        L1 = self._get_core_length_int(A_int, A_o_int)
        L3 = self._get_no_flow_height(L2, A_fr_int, L1, A_fr_ext)
        return A_int, A_ext, L1, L2, L3
