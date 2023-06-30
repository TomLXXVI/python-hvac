"""Sizing of a circular finned tube counterflow or crossflow heat exchanger for
single-phase fluids.

References
----------
[1] Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger Design.
    John Wiley & Sons. (§ 9.2.2 and § 9.3).
"""
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid, FluidState
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube
from hvac.heat_transfer.forced_convection.external_flow.general import prandtl_number
from hvac.heat_transfer.finned_surface.fins import CircularRectangularFin
from hvac.heat_transfer.heat_exchanger import eps_ntu
from hvac.heat_transfer.heat_exchanger.misc import correct_friction_factor
from ..core.geometry import TubeBankInside, TTubeBankOutside, CircularFinStaggeredTBO
from ..correlations.circular_fin.heat_transfer import (
    j_Briggs_and_Young,
    Nu_Briggs_and_Young
)
from ..correlations.circular_fin.friction_factor import (
    f_Chai,
    convert_friction_factor
)


Q_ = Quantity


class CircularFinTubeCounterFlowHeatExchanger:
    """Represent a single-pass, staggered, fin-tube heat exchanger with
    equilateral triangular pitch and with circular fins in a counterflow
    arrangement.
    """
    def __init__(
        self,
        D_i: Quantity,
        D_o: Quantity,
        S_t: Quantity,
        S_l: Quantity,
        D_f: Quantity,
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
        D_f:
            Fin tip diameter.
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
        self.geometry_int = self._create_internal_geometry(
            S_t, S_l, D_i, alpha_int
        )
        self.geometry_ext = self._create_external_geometry(
            S_t, S_l, D_o, D_f, t_f,
            N_f, alpha_ext, sigma_ext,
            None, None, None
        )
        self.L1 = L1
        self.L3 = L3
        self.k_f = k_f
        self.eta_ext = 0.8  # initial guess
        self._tube_int = CircularTube(
            Di=self.geometry_int.D_i,
            L=Q_(0.0, 'm'),
            fluid=None
        )
        self.Fluid_int: Fluid | None = None
        self.Fluid_ext: Fluid | None = None
        self.m_dot_int: Quantity | None = None
        self.m_dot_ext: Quantity | None = None
        self.T_int_in: Quantity | None = None
        self.T_ext_in: Quantity | None = None
        self.P_int_in: Quantity | None = None
        self.P_ext_in: Quantity | None = None
        self.dP_int: Quantity | None = None
        self.dP_ext: Quantity | None = None
        self.eps: float | None = None
        self.Q: Quantity | None = None
        self.fluid_int_in: FluidState | None = None
        self.fluid_int_out: FluidState | None = None
        self.fluid_ext_in: FluidState | None = None
        self.fluid_ext_out: FluidState | None = None
        self.ht_params: eps_ntu.CounterFlowHeatExchanger | None = None

    def set_design_conditions(
        self,
        Fluid_int: Fluid,
        Fluid_ext: Fluid,
        m_dot_int: Quantity,
        m_dot_ext: Quantity,
        T_int_in: Quantity,
        T_ext_in: Quantity,
        P_int_in: Quantity,
        P_ext_in: Quantity,
        dP_int: Quantity,
        dP_ext: Quantity,
        eps: Quantity | None = None,
        Q: Quantity | None = None
    ) -> None:
        """Set the design conditions for which the heat exchanger will be
        sized.

        Parameters
        ----------
        Fluid_int:
            The type of internal fluid.
        Fluid_ext:
            The type of external fluid.
        m_dot_int:
            Mass flow rate of internal fluid.
        m_dot_ext:
            Mass flow rate of external fluid.
        T_int_in:
            Inlet temperature of internal fluid.
        T_ext_in:
            Inlet temperature of external fluid.
        P_int_in:
            Inlet pressure of internal fluid.
        P_ext_in:
            Outlet pressure of external fluid.
        dP_int:
            Design pressure drop between inlet and outlet at internal side.
        dP_ext:
            Design pressure drop between inlet and outlet at external side.
        eps: optional
            Desired heat exchanger effectiveness. If `eps` is left to None,
            parameter `Q` must be set.
        Q: optional
            Desired heat transfer rate between internal and external fluids.
            If `Q` is left to None, parameter `eps` must be set.

        Returns
        -------
        None
        """
        self.Fluid_int = Fluid_int
        self.Fluid_ext = Fluid_ext
        self.m_dot_int = m_dot_int.to('kg / s')
        self.m_dot_ext = m_dot_ext.to('kg / s')
        self.T_int_in = T_int_in.to('K')
        self.T_ext_in = T_ext_in.to('K')
        self.P_int_in = P_int_in.to('Pa')
        self.P_ext_in = P_ext_in.to('Pa')
        self.dP_int = dP_int.to('Pa')
        self.dP_ext = dP_ext.to('Pa')
        self.eps = eps
        self.Q = Q.to('W') if Q is not None else None

    # noinspection PyUnresolvedReferences
    def size(
        self,
        i_max: int = 5,
        tol_dP: Quantity = Q_(1, 'Pa')
    ) -> tuple[TubeBankInside, CircularFinStaggeredTBO]:
        """Runs the sizing routine, which is an iterative procedure.

        Parameters
        ----------
        i_max: default 5
            Maximum number of iterations to perform. If no final solution
            is found after `i_max` iterations, a `ValueError` exception is
            raised.
        tol_dP: default 1 Pa
            Allowable deviation between the current and previous calculated
            pressure drops on internal and external side of the heat exchanger
            for which a final solution is accepted.

        Returns
        -------
        The internal and external geometry and dimensions of the heat exchanger,
        being objects of class `TubeBankInside`, respectively
        `CircularFinStaggeredTBO` (see tube_fin_heat_exchanger.geometry.py).

        Attributes
        ----------
        After calling `size` the following attributes are also available on the
        instance:

        fluid_int_in: FluidState
            The state of the internal fluid at the inlet.
        fluid_int_out: FluidState
            The state of the internal fluid at the outlet.
        fluid_ext_in: FluidState
            The state of the external fluid at the inlet.
        fluid_ext_out: FluidState
            The state of the external fluid at the outlet.
        ht_params: eps_ntu.CounterFlowHeatExchanger | eps_ntu.CrossFlowHeatExchanger
            Instance holding the heat transfer parameters from the eps-ntu
            method (NTU, eps, Q, ...).
        """
        dP_int, dP_ext = self.dP_int, self.dP_ext
        eta_ext = self.eta_ext
        Re_int, Re_ext = None, None
        N_r = 1.0
        i = 0
        while i <= i_max:
            # determine fluid properties
            fluid_int, fluid_ext = self._determine_fluid_properties(
                dP_int,
                dP_ext
            )
            self._tube_int.fluid = fluid_int
            # determine outlet temperatures
            # T_int_out, T_ext_out = self._get_fluid_outlet_temperatures(
            #     fluid_int.cp,
            #     fluid_ext.cp
            # )
            # determine heat transfer parameters (with eps-NTU method)
            ht_params = self._get_heat_transfer_params(fluid_int, fluid_ext)
            NTU, C_r, C_min = ht_params.NTU, ht_params.C_r, ht_params.C_min
            if self._get_fluid('hot') is self.Fluid_int:
                T_int_out = ht_params.T_hot_out
                T_ext_out = ht_params.T_cold_out
            else:
                T_int_out = ht_params.T_cold_out
                T_ext_out = ht_params.T_hot_out
            # determine mass flow velocity on internal side
            ntu_int = self._get_ntu_int(fluid_int, fluid_ext, NTU, C_r)
            jf_ratio_int = self._get_jf_ratio_int(fluid_int, Re_int)
            G_int = self._get_mass_velocity_int(
                fluid_int,
                T_int_out,
                ntu_int,
                jf_ratio_int
            )
            # determine mass flow velocity on external side
            ntu_ext = self._get_ntu_ext(fluid_int, fluid_ext, NTU, C_r)
            jf_ratio_ext = self._get_jf_ratio_ext(Re_ext, N_r)
            G_ext = self._get_mass_velocity_ext(
                fluid_ext,
                T_ext_out,
                ntu_ext,
                eta_ext,
                jf_ratio_ext
            )
            # determine reynolds number for next iteration
            Re_int = self._get_reynolds_number_int(fluid_int, G_int)
            Re_ext = self._get_reynolds_number_ext(
                fluid_ext, G_ext,
                self.geometry_ext.D_h
            )
            # determine unit conductance referred to external heat transfer area
            h_int = self._get_int_heat_trf_coeff(fluid_int, G_int)
            h_ext = self._get_ext_heat_trf_coeff(fluid_ext, G_ext, N_r)
            eta_ext = self._get_overall_fin_efficiency(h_ext)
            U_ext = self._get_unit_conductance_ext(h_int, h_ext, eta_ext)
            # get core dimensions
            A_int, A_ext, L1, L2, L3 = self._get_core_dimensions(
                NTU, C_min, U_ext, G_ext, G_int
            )

            N_r = (L2.to('m') / self.geometry_ext.S_l).to('m / m').m
            # calculate internal pressure drop
            dP_int_new = self._get_pressure_drop_int(
                L=self._get_int_flow_length(L1, L2),
                D_h=self.geometry_int.D_h,
                h_int=h_int,
                A_int=A_int,
                fluid_int=fluid_int,
                G_int=G_int,
                h_ext=h_ext,
                A_ext=A_ext,
                fluid_ext=fluid_ext
            )
            # calculate external pressure drop
            dP_ext_new = self._get_pressure_drop_ext(
                L=L2,
                D_h=self.geometry_ext.D_h,
                sigma=self.geometry_ext.sigma,
                h_int=h_int,
                A_int=A_int,
                fluid_int=fluid_int,
                h_ext=h_ext,
                A_ext=A_ext,
                fluid_ext=fluid_ext,
                G_ext=G_ext,
                T_ext_out=T_ext_out,
                P_ext_out=self.P_ext_in - self.dP_ext,
                N_r=N_r
            )
            # check stop criterium
            dev_dP_int = abs(dP_int - dP_int_new)
            dev_dP_ext = abs(dP_ext - dP_ext_new)
            if (dev_dP_int <= tol_dP) and (dev_dP_ext <= tol_dP):
                # define the complete internal heat exchanger geometry and
                # dimensions
                self.geometry_int = TubeBankInside(
                    S_t=self.geometry_int.S_t,
                    S_l=self.geometry_int.S_l,
                    D_i=self.geometry_int.D_i,
                    L1=L1.to('m'),
                    L2=L2.to('m'),
                    L3=L3.to('m')
                )
                # define the complete external heat exchanger geometry and
                # dimensions
                self.geometry_ext = self._create_external_geometry(
                    S_t=self.geometry_ext.S_t,
                    S_l=self.geometry_ext.S_l,
                    D_o=self.geometry_ext.D_o,
                    D_f=(
                        self.geometry_ext.D_f
                        if isinstance(self.geometry_ext, CircularFinStaggeredTBO)
                        else None
                    ),
                    t_f=self.geometry_ext.t_f,
                    N_f=self.geometry_ext.N_f,
                    alpha_ext=None,
                    sigma_ext=None,
                    L1=L1.to('m'),
                    L2=L2.to('m'),
                    L3=L3.to('m')
                )
                # set the final calculated pressure drops
                self.dP_int = dP_int_new
                self.dP_ext = dP_ext_new
                # determine the final fluid properties and outlet temperatures
                fluid_int, fluid_ext = self._determine_fluid_properties(
                    self.dP_int,
                    self.dP_ext
                )
                # T_int_out, T_ext_out = self._get_fluid_outlet_temperatures(
                #     fluid_int.cp,
                #     fluid_ext.cp
                # )
                # calculate the final heat transfer parameters
                self.ht_params = self._get_heat_transfer_params(fluid_int, fluid_ext)
                if self._get_fluid('hot') is self.Fluid_int:
                    T_int_out = ht_params.T_hot_out
                    T_ext_out = ht_params.T_cold_out
                else:
                    T_int_out = ht_params.T_cold_out
                    T_ext_out = ht_params.T_hot_out
                # set the state of the internal and external fluid at the inlet
                # and outlet of the heat exchanger
                self.fluid_int_in = self.Fluid_int(
                    T=self.T_int_in,
                    P=self.P_int_in
                )
                self.fluid_int_out = self.Fluid_int(
                    T=T_int_out,
                    P=self.P_int_in - self.dP_int
                )
                self.fluid_ext_in = self.Fluid_ext(
                    T=self.T_ext_in,
                    P=self.P_ext_in
                )
                self.fluid_ext_out = self.Fluid_ext(
                    T=T_ext_out,
                    P=self.P_ext_in - self.dP_ext
                )
                # return the internal and external heat exchanger geometry and
                # dimensions
                return self.geometry_int, self.geometry_ext
            dP_int = dP_int_new
            dP_ext = dP_ext_new
            i += 1
        else:
            raise ValueError(f"no solution after {i_max} iterations")

    def __call__(
        self,
        Fluid_int: Fluid,
        Fluid_ext: Fluid,
        m_dot_int: Quantity,
        m_dot_ext: Quantity,
        T_int_in: Quantity,
        T_ext_in: Quantity,
        P_int_in: Quantity,
        P_ext_in: Quantity,
        dP_int: Quantity,
        dP_ext: Quantity,
        eps: Quantity | None = None,
        Q: Quantity | None = None,
        tol_dP: Quantity = Q_(1, 'Pa'),
        i_max: int = 5
    ) -> tuple[TubeBankInside, CircularFinStaggeredTBO]:
        self.set_design_conditions(
            Fluid_int, Fluid_ext, m_dot_int, m_dot_ext, T_int_in, T_ext_in,
            P_int_in, P_ext_in, dP_int, dP_ext, eps, Q
        )
        geometry_int, geometry_ext = self.size(i_max, tol_dP)
        return geometry_int, geometry_ext

    @staticmethod
    def _create_internal_geometry(
        S_t: Quantity,
        S_l: Quantity,
        D_i: Quantity,
        alpha_int: Quantity
    ) -> TubeBankInside:
        geometry_int = TubeBankInside(
            S_t=S_t.to('m'),
            S_l=S_l.to('m'),
            D_i=D_i.to('m'),
            alpha=alpha_int.to('1 / m')
        )
        return geometry_int

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
        geometry_ext = CircularFinStaggeredTBO(
            S_t=S_t.to('m'),
            S_l=S_l.to('m'),
            D_o=D_o.to('m'),
            D_r=D_o.to('m'),
            D_f=D_f.to('m'),
            t_f=t_f.to('m'),
            N_f=N_f.to('1 / m'),
            alpha=alpha_ext.to('1 / m') if alpha_ext is not None else None,
            sigma=sigma_ext if sigma_ext is not None else None,
            L1=L1.to('m') if L1 is not None else None,
            L2=L2.to('m') if L2 is not None else None,
            L3=L3.to('m') if L3 is not None else None
        )
        return geometry_ext

    # noinspection PyUnusedLocal
    @staticmethod
    def _get_int_flow_length(L1: Quantity, L2: Quantity) -> Quantity:
        return L2

    def _get_fluid(self, which: str) -> Fluid:
        if self.T_int_in > self.T_ext_in:
            return self.Fluid_int if which == 'hot' else self.Fluid_ext
        else:
            return self.Fluid_ext if which == 'hot' else self.Fluid_int

    def _get_heat_transfer_rate(
        self,
        eps: Quantity,
        C_min: Quantity,
    ) -> Quantity:
        if self._get_fluid('hot') is self.Fluid_int:
            dT = self.T_int_in - self.T_ext_in
        else:
            dT = self.T_ext_in - self.T_int_in
        Q = eps * C_min * dT
        return Q

    def _get_fluid_outlet_temperatures(
        self,
        cp_int: Quantity,
        cp_ext: Quantity,
    ) -> tuple[Quantity, Quantity]:
        C_int = cp_int * self.m_dot_int
        C_ext = cp_ext * self.m_dot_ext
        if self.eps is not None:
            self.Q = self._get_heat_transfer_rate(
                eps=self.eps,
                C_min=min(C_int, C_ext)
            )
        if self.Q is None:
            raise ValueError('parameter `Q` and `eps` cannot be both None')
        if self._get_fluid('hot') is self.Fluid_int:
            T_int_out = self.T_int_in.to('K') - self.Q / C_int
            T_ext_out = self.T_ext_in.to('K') + self.Q / C_ext
        else:
            T_int_out = self.T_int_in.to('K') + self.Q / C_int
            T_ext_out = self.T_ext_in.to('K') - self.Q / C_ext
        return T_int_out, T_ext_out

    def _get_fluid_bulk_temperatures(
        self,
        cp_int: Quantity,
        cp_ext: Quantity,
        T_int_out: Quantity,
        T_ext_out: Quantity
    ) -> tuple[Quantity, Quantity]:
        C_int = cp_int * self.m_dot_int
        C_ext = cp_ext * self.m_dot_ext
        C_max = max(C_int, C_ext)
        C_min = min(C_int, C_ext)
        C_r = C_min / C_max
        if C_r >= 0.5:
            T_int_bulk = (self.T_int_in.to('K') + T_int_out.to('K')) / 2
            T_ext_bulk = (self.T_ext_in.to('K') + T_ext_out.to('K')) / 2
            return T_int_bulk, T_ext_bulk
        else:
            if self._get_fluid('hot') is self.Fluid_int:
                # fluid on internal side is hot fluid
                if C_max == C_int:
                    # internal fluid has the smallest temperature change
                    T_int_bulk = (self.T_int_in.to('K') + T_int_out.to('K')) / 2
                    DT_max = T_int_bulk - self.T_ext_in.to('K')
                    DT_min = T_int_bulk - T_ext_out.to('K')
                    LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
                    T_ext_bulk = T_int_bulk - LMTD
                    return T_int_bulk, T_ext_bulk
                else:  # C_max == C_ext
                    # external fluid has the smallest temperature change
                    T_ext_bulk = (self.T_ext_in.to('K') + T_ext_out.to('K')) / 2
                    DT_max = self.T_int_in.to('K') - T_ext_bulk
                    DT_min = T_int_out.to('K') - T_ext_bulk
                    LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
                    T_int_bulk = T_ext_bulk + LMTD
                    return T_int_bulk, T_ext_bulk
            else:
                # fluid on external side is hot fluid
                if C_max == C_int:
                    # internal fluid has the smallest temperature change
                    T_int_bulk = (self.T_int_in.to('K') + T_int_out.to('K')) / 2
                    DT_max = self.T_ext_in.to('K') - T_int_bulk
                    DT_min = T_ext_out.to('K') - T_int_bulk
                    LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
                    T_ext_bulk = T_int_bulk + LMTD
                    return T_int_bulk, T_ext_bulk
                else:  # C_max == C_ext
                    # external fluid has the smallest temperature change
                    T_ext_bulk = (self.T_ext_in.to('K') + T_ext_out.to('K')) / 2
                    DT_max = T_ext_bulk - self.T_int_in.to('K')
                    DT_min = T_ext_bulk - T_int_out.to('K')
                    LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
                    T_int_bulk = T_ext_bulk - LMTD
                    return T_int_bulk, T_ext_bulk

    @staticmethod
    def _get_fluid_state(
        T_bulk: Quantity,
        P: Quantity,
        fluid: Fluid
    ) -> FluidState:
        fluid_state = fluid(T=T_bulk, P=P)
        return fluid_state

    def _determine_fluid_properties(
        self,
        dP_int: Quantity,
        dP_ext: Quantity,
    ) -> tuple[FluidState, FluidState]:
        P_int_avg = (2 * self.P_int_in - dP_int) / 2
        P_ext_avg = (2 * self.P_ext_in - dP_ext) / 2
        fluid_int = self.Fluid_int(T=self.T_int_in, P=P_int_avg)
        fluid_ext = self.Fluid_ext(T=self.T_ext_in, P=P_ext_avg)
        i_max = 5
        i = 0
        tol = Q_(0.1, 'J / (kg * K)')
        while i < i_max:
            T_int_out, T_ext_out = self._get_fluid_outlet_temperatures(
                cp_int=fluid_int.cp,
                cp_ext=fluid_ext.cp
            )
            T_int_bulk, T_ext_bulk = self._get_fluid_bulk_temperatures(
                cp_int=fluid_int.cp,
                cp_ext=fluid_ext.cp,
                T_int_out=T_int_out,
                T_ext_out=T_ext_out
            )
            fluid_int_new = self._get_fluid_state(
                T_bulk=T_int_bulk,
                P=P_int_avg,
                fluid=self.Fluid_int
            )
            fluid_ext_new = self._get_fluid_state(
                T_bulk=T_ext_bulk,
                P=P_ext_avg,
                fluid=self.Fluid_ext
            )
            delta_cp_int = abs(fluid_int_new.cp - fluid_int.cp)
            delta_cp_ext = abs(fluid_ext_new.cp - fluid_ext.cp)
            if (delta_cp_int <= tol) and (delta_cp_ext <= tol):
                return fluid_int_new, fluid_ext_new
            fluid_int = fluid_int_new
            fluid_ext = fluid_ext_new
            i += 1
        else:
            raise ValueError("Fluid properties could not be determined.")

    def _get_heat_transfer_params(
        self,
        fluid_int: FluidState,
        fluid_ext: FluidState,
    ) -> eps_ntu.CounterFlowHeatExchanger:
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
        if self.eps is not None:
            self.Q = self._get_heat_transfer_rate(
                eps=self.eps,
                C_min=min(C_cold, C_hot)
            )
        if self.Q is None:
            raise ValueError('parameter `Q` and `eps` cannot be both None')
        cof_hex = eps_ntu.CounterFlowHeatExchanger(
            C_cold=C_cold,
            C_hot=C_hot,
            T_cold_in=T_cold_in,
            T_hot_in=T_hot_in,
            Q=self.Q
        )
        return cof_hex

    @staticmethod
    def _get_fluid_phase(fluid: FluidState) -> str:
        if 'liquid' in fluid.phase:
            return 'liquid'
        elif 'gas' in fluid.phase:
            return 'gas'
        else:
            raise ValueError('the fluid should be in a liquid or gaseous phase.')

    def _get_ntu_int(
        self,
        fluid_int: FluidState,
        fluid_ext: FluidState,
        NTU: float,
        C_r: float
    ) -> float:
        fluid_int_phase = self._get_fluid_phase(fluid_int)
        fluid_ext_phase = self._get_fluid_phase(fluid_ext)
        if fluid_int_phase == fluid_ext_phase:
            return 2 * NTU
        elif fluid_int_phase == 'liquid':
            return 10 * C_r * NTU
        elif fluid_int_phase == 'gas':
            return 1.11 * NTU
        else:
            raise ValueError('internal fluid is neither a liquid nor a gas')

    def _get_ntu_ext(
        self,
        fluid_int: FluidState,
        fluid_ext: FluidState,
        NTU: float,
        C_r: float
    ) -> float:
        fluid_int_phase = self._get_fluid_phase(fluid_int)
        fluid_ext_phase = self._get_fluid_phase(fluid_ext)
        if fluid_int_phase == fluid_ext_phase:
            return 2 * NTU
        elif fluid_ext_phase == 'liquid':
            return 10 * C_r * NTU
        elif fluid_ext_phase == 'gas':
            return 1.11 * NTU
        else:
            raise ValueError('external fluid is neither a liquid nor a gas')

    def _get_m_dot_tube(
        self,
        fluid_int: FluidState,
        Re_int: float
    ) -> Quantity:
        u_tube = Re_int * fluid_int.mu / (self._tube_int.Dh * fluid_int.rho)
        m_dot_tube = fluid_int.rho * self._tube_int.Ac * u_tube
        return m_dot_tube.to('kg / s')

    def _get_jf_ratio_int(self, fluid_int: FluidState, Re_int: float | None) -> float:

        def __jf_ratio(Re: float) -> float:
            m_dot_tube = self._get_m_dot_tube(fluid_int, Re)
            self._tube_int.m_dot = m_dot_tube
            Nu = self._tube_int.fully_dev_nusselt_number()
            if isinstance(Nu, tuple):
                Nu = Nu[0]  # laminar flow with uniform heat flux
            Pr = self._tube_int.prandtl_number()
            j = Nu * (Pr ** -(1 / 3)) / Re
            f = self._tube_int.fully_dev_friction_factor()
            jf_ratio = j / f
            return jf_ratio

        if Re_int is None:
            Re_rng = np.arange(100, 50_000, 100)
            jf_ratio_rng = np.array([__jf_ratio(Re) for Re in Re_rng])
            jf_ratio = float(np.mean(jf_ratio_rng))
        else:
            jf_ratio = __jf_ratio(Re_int)
        return jf_ratio

    def _get_jf_ratio_ext(
        self,
        Re_D_o_ext: float | None,
        N_r: float | None = None
    ) -> float:
        def __colburn_j_factor(Re_D_o: float) -> float:
            return j_Briggs_and_Young(
                Re_D_o=Re_D_o,
                N_f=self.geometry_ext.N_f.to('1 / m').m,
                t_f=self.geometry_ext.t_f.to('m').m,
                h_f=(self.geometry_ext.D_f - self.geometry_ext.D_r).to('m').m / 2
            )

        def __friction_factor(Re_D_o: float) -> float:
            f_tb = f_Chai(
                Re_D_o=Re_D_o,
                S_t=self.geometry_ext.S_t.to('m').m,
                S_l=self.geometry_ext.S_l.to('m').m,
                D_r=self.geometry_ext.D_r.to('m').m,
                D_f=self.geometry_ext.D_f.to('m').m,
                t_f=self.geometry_ext.t_f.to('m').m,
                N_f=self.geometry_ext.N_f.to('1 / m').m
            )
            f = convert_friction_factor(
                f_tb=f_tb,
                S_l=self.geometry_ext.S_l.to('m').m,
                D_h=self.geometry_ext.D_h.to('m').m
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

    def _get_mass_velocity_int(
        self,
        fluid_int: FluidState,
        T_int_out: Quantity,
        ntu_int: float,
        jf_ratio: float
    ) -> Quantity:
        Pr_int = prandtl_number(
            rho=fluid_int.rho,
            mu=fluid_int.mu,
            k=fluid_int.k,
            cp=fluid_int.cp
        )
        rho_int_in = self.Fluid_int(
            T=self.T_int_in,
            P=self.P_int_in
        ).rho.to('kg / m ** 3')
        rho_int_out = self.Fluid_int(
            T=T_int_out,
            P=self.P_int_in - self.dP_int
        ).rho.to('kg / m ** 3')
        rho_int_avg = 2 * rho_int_in * rho_int_out / (rho_int_in + rho_int_out)
        G = (
            2 * self.dP_int * jf_ratio /
            ((1 / rho_int_avg) * (Pr_int ** (2 / 3)) * ntu_int)
        ) ** 0.5
        return G.to('kg / (m ** 2 * s)')

    def _get_mass_velocity_ext(
        self,
        fluid_ext: FluidState,
        T_ext_out: Quantity,
        ntu_ext: float,
        eta_ext: float,
        jf_ratio: float
    ) -> Quantity:
        Pr_ext = prandtl_number(
            rho=fluid_ext.rho,
            mu=fluid_ext.mu,
            k=fluid_ext.k,
            cp=fluid_ext.cp
        )
        rho_ext_in = self.Fluid_ext(
            T=self.T_ext_in,
            P=self.P_ext_in
        ).rho.to('kg / m ** 3')
        rho_ext_out = self.Fluid_ext(
            T=T_ext_out,
            P=self.P_ext_in - self.dP_ext
        ).rho.to('kg / m ** 3')
        rho_ext_avg = 2 * rho_ext_in * rho_ext_out / (rho_ext_in + rho_ext_out)
        G = (
            2 * eta_ext * self.dP_ext * jf_ratio /
            ((1 / rho_ext_avg) * (Pr_ext ** (2 / 3)) * ntu_ext)
        ) ** 0.5
        return G.to('kg / (m ** 2 * s)')

    def _get_reynolds_number_int(self, fluid_int: FluidState, G_int: Quantity) -> float:
        G_int.ito('kg / (m ** 2 * s)')
        self._tube_int.Dh.ito('m')
        fluid_int.mu.ito('Pa * s')
        Re_int = (G_int * self._tube_int.Dh / fluid_int.mu).m
        return Re_int

    def _get_int_heat_trf_coeff(
        self,
        fluid_int: FluidState,
        G_int: Quantity
    ) -> Quantity:
        Re_int = self._get_reynolds_number_int(fluid_int, G_int)
        self._tube_int.m_dot = self._get_m_dot_tube(fluid_int, Re_int)
        h_int = self._tube_int.fully_dev_heat_transfer_coefficient()
        if isinstance(h_int, tuple):
            h_int = h_int[0]  # laminar flow with uniform heat flux
        return h_int

    def _get_int_friction_factor(
        self,
        fluid_int: FluidState,
        G_int: Quantity
    ) -> float:
        Re_int = self._get_reynolds_number_int(fluid_int, G_int)
        self._tube_int.m_dot = self._get_m_dot_tube(fluid_int, Re_int)
        f_int = self._tube_int.fully_dev_friction_factor() / 4  # Fanning friction factor
        return f_int

    @staticmethod
    def _get_reynolds_number_ext(
        fluid_ext: FluidState,
        G_ext: Quantity,
        D: Quantity
    ) -> float:
        G_ext.ito('kg / (m ** 2 * s)')
        fluid_ext.mu.ito('Pa * s')
        D.ito('m')
        Re_ext = (G_ext * D / fluid_ext.mu).m
        return Re_ext

    def _get_ext_heat_trf_coeff(
        self,
        fluid_ext: FluidState,
        G_ext: Quantity,
        N_r: float | None = None
    ) -> Quantity:
        Re_ext_D_o = self._get_reynolds_number_ext(
            fluid_ext, G_ext,
            self.geometry_ext.D_o
        )
        Pr_ext = prandtl_number(
            rho=fluid_ext.rho,
            mu=fluid_ext.mu,
            k=fluid_ext.k,
            cp=fluid_ext.cp
        )
        Nu_ext = Nu_Briggs_and_Young(
            Re_D_o=Re_ext_D_o,
            Pr=Pr_ext,
            N_f=self.geometry_ext.N_f.to('1 / m').m,
            t_f=self.geometry_ext.t_f.to('m').m,
            h_f=(self.geometry_ext.D_f - self.geometry_ext.D_r).to('m').m / 2
        )
        h_ext = Nu_ext * fluid_ext.k / self.geometry_ext.D_h
        return h_ext.to('W / (m ** 2 * K)')

    def _get_ext_friction_factor(
        self,
        fluid_ext: FluidState,
        G_ext: Quantity,
        N_r: float | None = None
    ) -> float:
        Re_ext_D_o = self._get_reynolds_number_ext(
            fluid_ext, G_ext,
            self.geometry_ext.D_o
        )
        f_tb = f_Chai(
            Re_D_o=Re_ext_D_o,
            S_t=self.geometry_ext.S_t.to('m').m,
            S_l=self.geometry_ext.S_l.to('m').m,
            D_r=self.geometry_ext.D_r.to('m').m,
            D_f=self.geometry_ext.D_f.to('m').m,
            t_f=self.geometry_ext.t_f.to('m').m,
            N_f=self.geometry_ext.N_f.to('1 / m').m
        )
        f = convert_friction_factor(
            f_tb,
            self.geometry_ext.S_l.to('m').m,
            self.geometry_ext.D_h.to('m').m
        )
        return f

    def _get_overall_fin_efficiency(self, h_ext) -> float:
        # get overall efficiency of finned surface
        fin = CircularRectangularFin(
            r_i=self.geometry_ext.D_r / 2,
            r_o=self.geometry_ext.D_f / 2,
            t=self.geometry_ext.t_f,
            k=self.k_f
        )
        fin.h_avg = h_ext
        eta_fin = fin.efficiency
        r = self.geometry_ext.A_f_to_A
        eta_ext = 1 - r * (1 - eta_fin)
        return eta_ext

    def _get_unit_conductance_ext(
        self,
        h_int: Quantity,
        h_ext: Quantity,
        eta_ext: float
    ) -> Quantity:
        # get unit conductance referred to external side
        R_int = 1 / h_int
        R_ext = 1 / (eta_ext * h_ext)
        alpha_ext = self.geometry_ext.alpha.to('1 / m')
        alpha_int = self.geometry_int.alpha.to('1 / m')
        R_t = (alpha_ext / alpha_int) * R_int + R_ext
        U_ext = 1 / R_t
        return U_ext.to('W / (m ** 2 * K)')

    @staticmethod
    def _get_heat_transfer_area_ext(NTU: float, C_min: float, U_ext: Quantity) -> Quantity:
        A_ext = NTU * C_min / U_ext
        return A_ext.to('m ** 2')

    def _get_heat_transfer_area_int(self, A_ext: Quantity) -> Quantity:
        alpha_ext = self.geometry_ext.alpha.to('1 / m')
        alpha_int = self.geometry_int.alpha.to('1 / m')
        A_int = (alpha_int / alpha_ext) * A_ext
        return A_int.to('m ** 2')

    def _get_free_flow_area_ext(self, G_ext: Quantity) -> Quantity:
        A_o_ext = self.m_dot_ext / G_ext
        return A_o_ext.to('m ** 2')

    def _get_free_flow_area_int(self, G_int: Quantity) -> Quantity:
        A_o_int = self.m_dot_int / G_int
        return A_o_int.to('m ** 2')

    def _get_frontal_area_ext(self, A_o_ext: Quantity) -> Quantity:
        A_fr_ext = A_o_ext / self.geometry_ext.sigma
        return A_fr_ext.to('m ** 2')

    def _get_frontal_area_int(self, A_o_int: Quantity) -> Quantity:
        A_fr_int = A_o_int / self.geometry_int.sigma
        return A_fr_int.to('m ** 2')

    def _get_core_length_ext(self, A_ext: Quantity, A_o_ext: Quantity) -> Quantity:
        D_h_ext = self.geometry_ext.D_h
        L2 = D_h_ext * A_ext / (4 * A_o_ext)
        return L2.to('m')

    def _get_core_length_int(self, A_int: Quantity, A_o_int: Quantity) -> Quantity:
        L1 = self.geometry_int.D_h * A_int / (4 * A_o_int)
        return L1.to('m')

    @staticmethod
    def _get_no_flow_height(
        L2: Quantity,
        A_fr_int: Quantity,
        L1: Quantity,
        A_fr_ext: Quantity
    ) -> Quantity:
        L3_1 = A_fr_int / L2
        L3_2 = A_fr_ext / L1
        L3 = (L3_1 + L3_2) / 2
        return L3

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
        # determine frontal area
        A_fr = max(A_fr_ext, A_fr_int)
        # determine core dimensions
        L1 = self.L1  # width
        L2 = self._get_core_length_ext(A_ext, A_o_ext)  # depth
        L3 = self.L3  # height
        if L1 is not None:
            # determine no flow (stack) height
            L3 = A_fr / L1.to('m')
        if L3 is not None:
            # determine internal core length
            L1 = A_fr / L3.to('m')
        return A_int, A_ext, L1, L2, L3

    @staticmethod
    def _get_int_resistance(h_int: Quantity, A_int: Quantity) -> Quantity:
        # get thermal resistance on internal side of heat exchanger
        h_int.ito('W / (m ** 2 * K)')
        A_int.ito('m ** 2')
        R_int = 1 / (h_int * A_int)  # K / W
        return R_int

    @staticmethod
    def _get_ext_resistance(h_ext: Quantity, eta_ext: float, A_ext: Quantity) -> Quantity:
        # get thermal resistance on external side of heat exchanger
        h_ext.ito('W / (m ** 2 * K)')
        A_ext.ito('m ** 2')
        R_ext = 1 / (eta_ext * h_ext * A_ext)  # K / W
        return R_ext

    def _get_wall_temperature(
        self,
        h_int: Quantity,
        A_int: Quantity,
        fluid_int: FluidState,
        h_ext: Quantity,
        A_ext: Quantity,
        fluid_ext: FluidState
    ) -> Quantity:
        # get the average wall temperature
        R_int = self._get_int_resistance(h_int, A_int)
        eta_ext = self._get_overall_fin_efficiency(h_ext)
        R_ext = self._get_ext_resistance(h_ext, eta_ext, A_ext)
        T_int = fluid_int.T.to('K')
        T_ext = fluid_ext.T.to('K')
        n = T_int / R_int + T_ext / R_ext
        d = 1 / R_int + 1 / R_ext
        T_w = n / d
        return T_w

    def _get_int_thermal_regime(self) -> str:
        # get thermal regime at internal side
        if self.T_int_in > self.T_ext_in:
            # internal fluid is the hot fluid giving off heat to cold fluid on
            # external side.
            return "cooling"
        else:
            return "heating"

    def _get_pressure_drop_int(
        self,
        L: Quantity,
        D_h: Quantity,
        h_int: Quantity,
        A_int: Quantity,
        fluid_int: FluidState,
        G_int: Quantity,
        h_ext: Quantity,
        A_ext: Quantity,
        fluid_ext: FluidState
    ) -> Quantity:
        Re_int = self._get_reynolds_number_int(fluid_int, G_int)
        self._tube_int.m_dot = self._get_m_dot_tube(fluid_int, Re_int)
        f = self._get_int_friction_factor(fluid_int, G_int)
        f_corr = correct_friction_factor(
            f_cp=f,
            T_w=self._get_wall_temperature(
                h_int, A_int, fluid_int,
                h_ext, A_ext, fluid_ext
            ),
            flow_regime=self._tube_int.get_flow_condition(),
            thermal_regime=self._get_int_thermal_regime(),
            fluid=fluid_int
        )
        u_m = self._tube_int.mean_velocity()
        rho_int = fluid_int.rho
        dP_int = f_corr * (L / D_h) * rho_int * u_m ** 2 / 2
        return dP_int.to('Pa')

    def _get_ext_flow_regime(self, fluid_ext: FluidState, G_ext: Quantity) -> str:
        Re_ext_D_h = self._get_reynolds_number_ext(
            fluid_ext, G_ext,
            self.geometry_ext.D_h
        )
        if Re_ext_D_h > 2300:
            return 'turbulent'
        else:
            return 'laminar'

    def _get_ext_thermal_regime(self) -> str:
        # get thermal regime at external side
        if self.T_int_in > self.T_ext_in:
            # internal fluid is the hot fluid giving off heat to cold fluid on
            # external side.
            return "heating"
        else:
            return "cooling"

    def _get_pressure_drop_ext(
        self,
        L: Quantity,
        D_h: Quantity,
        sigma: Quantity,
        h_int: Quantity,
        A_int: Quantity,
        fluid_int: FluidState,
        h_ext: Quantity,
        A_ext: Quantity,
        fluid_ext: FluidState,
        G_ext: Quantity,
        T_ext_out: Quantity,
        P_ext_out: Quantity,
        N_r: float | None = None
    ) -> Quantity:
        f = self._get_ext_friction_factor(fluid_ext, G_ext, N_r)
        f_corr = correct_friction_factor(
            f_cp=f,
            T_w=self._get_wall_temperature(
                h_int, A_int, fluid_int,
                h_ext, A_ext, fluid_ext
            ),
            flow_regime=self._get_ext_flow_regime(fluid_ext, G_ext),
            thermal_regime=self._get_ext_thermal_regime(),
            fluid=fluid_ext
        )
        fluid_ext_in = self._get_fluid_state(self.T_ext_in, self.P_ext_in, self.Fluid_ext)
        fluid_ext_out = self._get_fluid_state(T_ext_out, P_ext_out, self.Fluid_ext)
        rho_ext_avg = (
            2 * fluid_ext_in.rho * fluid_ext_out.rho
            / (fluid_ext_in.rho + fluid_ext_out.rho)
        )
        k = (G_ext ** 2) / (2 * fluid_ext_in.rho)
        a = f_corr * (4 * L / D_h) * (fluid_ext_in.rho / rho_ext_avg)
        b = (1 + sigma ** 2) * (fluid_ext_in.rho / fluid_ext_out.rho - 1)
        dP_ext = k * (a + b)
        return dP_ext.to('Pa')


class CircularFinTubeCrossFlowHeatExchanger(CircularFinTubeCounterFlowHeatExchanger):
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
        D_f: Quantity,
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
        D_f:
            Fin tip diameter.
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
            D_f=D_f,
            t_f=t_f,
            N_f=N_f,
            alpha_int=alpha_int,
            alpha_ext=alpha_ext,
            sigma_ext=sigma_ext,
            L1=None,
            L3=None,
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
