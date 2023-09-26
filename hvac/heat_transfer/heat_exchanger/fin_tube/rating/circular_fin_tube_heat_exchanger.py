"""
Rating of a circular finned tube crossflow heat exchanger for single-phase
fluids.

References
----------
[1] Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger Design.
    John Wiley & Sons. (§ 9.2.1 and § 9.3).
"""
import math
from hvac import Quantity
from hvac.fluids import Fluid, FluidState
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube
from hvac.heat_transfer.forced_convection.external_flow.general import prandtl_number
from hvac.heat_transfer.finned_surface.fins import (
    Fin,
    CircularRectangularFin
)
from hvac.heat_transfer.heat_exchanger.eps_ntu import (
    CrossFlowHeatExchanger,
    CounterFlowHeatExchanger
)
from hvac.heat_transfer.heat_exchanger.misc import (
    correct_nusselt_number,
    correct_friction_factor
)
from ..core.geometry import (
    TubeBankInside,
    TTubeBankOutside,
    CircularFinStaggeredTBO
)
from ..correlations.circular_fin.heat_transfer import Nu_Briggs_and_Young
from ..correlations.circular_fin.friction_factor import (
    f_Chai,
    convert_friction_factor
)


Q_ = Quantity


class CircularFinTubeCrossFlowHeatExchanger:
    """Represent a single-pass, staggered, fin-tube heat exchanger with
    equilateral triangular pitch and circular fins in a cross-flow arrangement.

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
    D_f:
        Fin tip diameter.
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
        D_f: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        k_f: Quantity = Q_(237, 'W / (m * K)'),
        t_header: Quantity = Q_(0.0, 'mm')
    ) -> None:
        # geometrical properties and dimensions of heat exchanger are set on
        # instantiation
        self.geometry_int = self._create_internal_geometry(
            S_t, S_l, D_i,
            L1, L2, L3, t_header
        )
        self.geometry_ext: CircularFinStaggeredTBO = self._create_external_geometry(
            S_t, S_l, D_o, D_r, D_f, t_f, N_f,
            L1, L2, L3
        )
        self.k_f = k_f
        self._tube_int = CircularTube(
            Di=self.geometry_int.D_i,
            L=self.geometry_int.L1,
            fluid=None,
            e=Q_(0.0, 'mm')
        )

        # operating conditions of heat exchanger are set after instantiation; see
        # method `set_operating_conditions`.
        self.m_dot_int: Quantity | None = None
        self.m_dot_ext: Quantity | None = None
        self.T_int_in: Quantity | None = None
        self.T_ext_in: Quantity | None = None
        self.P_int_in: Quantity | None = None
        self.P_ext_in: Quantity | None = None
        self.Fluid_int: Fluid | None = None
        self.Fluid_ext: Fluid | None = None
        self.eps_ini: float | None = None

        # performance, results available after calling `get_performance`
        self.T_int_out: Quantity | None = None
        self.T_ext_out: Quantity | None = None
        self.Q: Quantity | None = None
        self.eps: float | None = None
        self.UA: Quantity | None = None
        self.dP_int: Quantity | None = None
        self.dP_ext: Quantity | None = None
        self.v_fa: Quantity | None = None

    def set_operating_conditions(
        self,
        m_dot_int: Quantity,
        m_dot_ext: Quantity,
        T_int_in: Quantity,
        T_ext_in: Quantity,
        P_int_in: Quantity,
        P_ext_in: Quantity,
        Fluid_int: Fluid,
        Fluid_ext: Fluid,
        eps_ini: float
    ) -> None:
        """Set the operating conditions for the heat exchanger.

        Parameters
        ----------
        m_dot_int:
            Mass flow rate of fluid on internal side of heat exchanger.
        m_dot_ext:
            Mass flow rate of fluid on external side of heat exchanger.
        T_int_in:
            Inlet temperature of internal fluid.
        T_ext_in:
            Inlet temperature of external fluid.
        P_int_in:
            Inlet pressure of internal fluid.
        P_ext_in:
            Inlet pressure of external fluid.
        Fluid_int:
            Type of fluid on internal side.
        Fluid_ext
            Type of fluid on external side.
        eps_ini:
            initial guess for the heat exchanger effectiveness (for most
            single-pass crossflow exchangers assume 50-75 %)

        Returns
        -------
        None
        """
        self.m_dot_int = m_dot_int.to('kg / s')
        self.m_dot_ext = m_dot_ext.to('kg / s')
        self.T_int_in = T_int_in.to('K')
        self.T_ext_in = T_ext_in.to('K')
        self.P_int_in = P_int_in.to('Pa')
        self.P_ext_in = P_ext_in.to('Pa')
        self.Fluid_int = Fluid_int
        self.Fluid_ext = Fluid_ext
        self.eps_ini = eps_ini

    def rate(
        self,
        i_max: int = 5,
        tol_T_out: Quantity = Q_(0.1, 'K')
    ) -> tuple[Quantity, Quantity, Quantity]:
        """Returns for the given set of operating conditions the heat
        transfer rate between the internal and external fluid in the heat
        exchanger and their outlet temperatures. The solution is determined
        by iterative calculations.

        Parameters
        ----------
        i_max:
            Maximum number of iterations.
        tol_T_out:
            Allowable temperature difference between the current and previous
            calculated outlet temperatures at which the iterative loop will
            be exited.

        Raises
        ------
        ValueError
            If no solution within the given tolerance is found after the
            maximum number of iterations.

        Attributes
        ----------
        Additional parameters can be retrieved when the rating routine has
        finished.
        self.eps: float
            Effectiveness of the heat exchanger.
        self.UA: Quantity
            Overall conductance of the heat exchanger.
        self.dP_int: Quantity
            Pressure drop on the internal side of the heat exchanger.
        self.dP_ext: Quantity
            Pressure drop on the external side of the heat exchanger.
        self.v_fa: Quantity
            Face velocity of external fluid.
        """
        # determine an initial guess for the outlet temperatures
        P_int_out, P_ext_out = self.P_int_in, self.P_ext_in
        cp_int = self.Fluid_int(T=self.T_int_in, P=self.P_int_in).cp
        cp_ext = self.Fluid_ext(T=self.T_ext_in, P=self.P_ext_in).cp
        T_int_out, T_ext_out = self._get_fluid_outlet_temperatures(
            cp_int,
            cp_ext,
            self.eps_ini
        )
        i = 0
        while i < i_max:
            fluid_int, fluid_ext = self._determine_fluid_properties(
                T_int_out, P_int_out,
                T_ext_out, P_ext_out
            )
            h_int, h_ext = self._get_heat_transfer_coeffs(fluid_int, fluid_ext)
            UA = self._get_overall_conductance(h_int, h_ext)
            T_int_out_new, T_ext_out_new, *ht_params = self._get_hex_performance(
                fluid_int,
                fluid_ext,
                UA
            )
            dP_int = self._get_int_pressure_drop(
                fluid_int, h_int,
                fluid_ext, h_ext
            )
            dP_ext = self._get_ext_pressure_drop(
                fluid_int, h_int,
                fluid_ext, h_ext, T_ext_out_new
            )
            dev_T_int_out = abs(T_int_out_new.to('K') - T_int_out.to('K'))
            dev_T_ext_out = abs(T_ext_out_new.to('K') - T_ext_out.to('K'))
            if (dev_T_int_out <= tol_T_out) and (dev_T_ext_out <= tol_T_out):
                self.T_int_out = T_int_out_new
                self.T_ext_out = T_ext_out_new
                self.Q = ht_params[0]
                self.eps = ht_params[1]
                self.UA = ht_params[2]
                self.dP_int = dP_int
                self.dP_ext = dP_ext
                self.v_fa = self._get_face_velocity()
                return self.T_int_out, self.T_ext_out, self.Q
            P_int_out = self.P_int_in - dP_int
            P_ext_out = self.P_ext_in - dP_ext
            T_int_out = T_int_out_new
            T_ext_out = T_ext_out_new
            i += 1
        else:
            raise ValueError(
                "no solution found for outlet temperatures within "
                f"tolerance {tol_T_out.to('K'):~P} after {i_max} iterations"
            )

    def __call__(
        self,
        m_dot_int: Quantity,
        m_dot_ext: Quantity,
        T_int_in: Quantity,
        T_ext_in: Quantity,
        P_int_in: Quantity,
        P_ext_in: Quantity,
        Fluid_int: Fluid,
        Fluid_ext: Fluid,
        eps_ini: float,
        i_max: int = 5,
        tol_T_out: Quantity = Q_(0.1, 'K')
    ) -> tuple[Quantity, Quantity, Quantity]:
        """Combines the methods `set_operating_conditions` and `rate` in a
        single call.
        """
        self.set_operating_conditions(
            m_dot_int, m_dot_ext, T_int_in, T_ext_in, P_int_in, P_ext_in,
            Fluid_int, Fluid_ext, eps_ini)
        T_int_out, T_ext_out, Q = self.rate(i_max, tol_T_out)
        return T_int_out, T_ext_out, Q

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
            S_t=S_t, S_l=S_l,
            D_i=D_i,
            L1=L1, L2=L2, L3=L3,
            t_header=t_header
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
        geometry_ext = CircularFinStaggeredTBO(
            S_t=S_t, S_l=S_l,
            D_o=D_o, D_r=D_r, D_f=D_f,
            t_f=t_f, N_f=N_f,
            L1=L1, L2=L2, L3=L3
        )
        return geometry_ext

    def _get_fluid(self, which: str) -> Fluid:
        """Returns the type of fluid indicated by `which`.
        Parameter `which` can either be 'hot' or 'cold'.
        In case `which` is 'hot', the hot fluid will be returned, i.e. the fluid
        on the internal or external side that gives off heat and reduces in
        temperature between inlet and outlet.
        In case `which` is 'cold', the cold fluid will be returned, i.e. the fluid
        on the internal or external side that receives heat and raises in
        temperature between inlet and outlet.
        """
        if self.T_int_in > self.T_ext_in:
            return self.Fluid_int if which == 'hot' else self.Fluid_ext
        else:
            return self.Fluid_ext if which == 'hot' else self.Fluid_int

    def _get_heat_transfer_rate(
        self,
        eps: Quantity,
        C_min: Quantity,
    ) -> Quantity:
        """Returns the heat transfer rate between the internal and external
        fluid, knowing the effectiveness of the heat exchanger and the
        minimum capacitance rate.
        """
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
        eps: float
    ) -> tuple[Quantity, Quantity]:
        """Returns the internal fluid outlet temperature and the external fluid
        outlet temperature."""
        C_int = cp_int * self.m_dot_int
        C_ext = cp_ext * self.m_dot_ext
        Q = self._get_heat_transfer_rate(
            eps=eps,
            C_min=min(C_int, C_ext)
        )
        if self._get_fluid('hot') is self.Fluid_int:
            T_int_out = self.T_int_in.to('K') - Q / C_int
            T_ext_out = self.T_ext_in.to('K') + Q / C_ext
        else:
            T_int_out = self.T_int_in.to('K') + Q / C_int
            T_ext_out = self.T_ext_in.to('K') - Q / C_ext
        return T_int_out, T_ext_out

    def _get_fluid_bulk_temperatures(
        self,
        cp_int: Quantity,
        cp_ext: Quantity,
        T_int_out: Quantity,
        T_ext_out: Quantity
    ) -> tuple[Quantity, Quantity]:
        """Returns the internal fluid mean or bulk temperature and the external
        fluid mean or bulk temperature."""
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
                    LMTD = (DT_max - DT_min) / math.log(DT_max / DT_min)
                    T_ext_bulk = T_int_bulk - LMTD
                    return T_int_bulk, T_ext_bulk
                else:  # C_max == C_ext
                    # external fluid has the smallest temperature change
                    T_ext_bulk = (self.T_ext_in.to('K') + T_ext_out.to('K')) / 2
                    DT_max = self.T_int_in.to('K') - T_ext_bulk
                    DT_min = T_int_out.to('K') - T_ext_bulk
                    LMTD = (DT_max - DT_min) / math.log(DT_max / DT_min)
                    T_int_bulk = T_ext_bulk + LMTD
                    return T_int_bulk, T_ext_bulk
            else:
                # fluid on external side is hot fluid
                if C_max == C_int:
                    # internal fluid has the smallest temperature change
                    T_int_bulk = (self.T_int_in.to('K') + T_int_out.to('K')) / 2
                    DT_max = self.T_ext_in.to('K') - T_int_bulk
                    DT_min = T_ext_out.to('K') - T_int_bulk
                    LMTD = (DT_max - DT_min) / math.log(DT_max / DT_min)
                    T_ext_bulk = T_int_bulk + LMTD
                    return T_int_bulk, T_ext_bulk
                else:  # C_max == C_ext
                    # external fluid has the smallest temperature change
                    T_ext_bulk = (self.T_ext_in.to('K') + T_ext_out.to('K')) / 2
                    DT_max = T_ext_bulk - self.T_int_in.to('K')
                    DT_min = T_ext_bulk - T_int_out.to('K')
                    LMTD = (DT_max - DT_min) / math.log(DT_max / DT_min)
                    T_int_bulk = T_ext_bulk - LMTD
                    return T_int_bulk, T_ext_bulk

    @staticmethod
    def _get_fluid_state(
        T_bulk: Quantity,
        P: Quantity,
        fluid: Fluid
    ) -> FluidState:
        """Returns the state of the given type of fluid when the bulk
        temperature and pressure are given.
        """
        fluid_state = fluid(T=T_bulk, P=P)
        return fluid_state

    def _determine_fluid_properties(
        self,
        T_int_out: Quantity,
        P_int_out: Quantity,
        T_ext_out: Quantity,
        P_ext_out: Quantity
    ) -> tuple[FluidState, FluidState]:
        """Returns the mean or bulk state of the fluid on the internal and
        external side of the heat exchanger.
        """
        T_int_avg = (self.T_int_in.to('K') + T_int_out.to('K')) / 2
        T_ext_avg = (self.T_ext_in.to('K') + T_ext_out.to('K')) / 2
        P_int_avg = (self.P_int_in + P_int_out) / 2
        P_ext_avg = (self.P_ext_in + P_ext_out) / 2
        fluid_int = self.Fluid_int(T=T_int_avg, P=P_int_avg)
        fluid_ext = self.Fluid_ext(T=T_ext_avg, P=P_ext_avg)
        i_max = 5
        i = 0
        tol = Q_(0.1, 'J / (kg * K)')
        while i < i_max:
            T_int_bulk, T_ext_bulk = self._get_fluid_bulk_temperatures(
                cp_int=fluid_int.cp,
                cp_ext=fluid_ext.cp,
                T_int_out=T_int_out,
                T_ext_out=T_ext_out
            )
            fluid_int_new = self._get_fluid_state(
                T_bulk=T_int_bulk,
                P=self.P_int_in,
                fluid=self.Fluid_int
            )
            fluid_ext_new = self._get_fluid_state(
                T_bulk=T_ext_bulk,
                P=self.P_ext_in,
                fluid=self.Fluid_ext
            )
            delta_cp_int = abs(fluid_int_new.cp - fluid_int.cp)
            delta_cp_ext = abs(fluid_ext_new.cp - fluid_ext.cp)
            if (delta_cp_int <= tol) and (delta_cp_ext <= tol):
                self._tube_int.fluid = fluid_int_new  # assign internal fluid state to the tube
                return fluid_int_new, fluid_ext_new
            fluid_int = fluid_int_new
            fluid_ext = fluid_ext_new
            i += 1
        else:
            raise ValueError("fluid properties could not be determined")

    def _get_m_dot_tube(self, fluid_in: FluidState) -> Quantity:
        """Returns the mass flow rate of internal fluid in a single tube."""
        # determine mean velocity in tubes
        rho_int = fluid_in.rho.to('kg / m ** 3')
        m_dot_int = self.m_dot_int.to('kg / s')
        A_o_int = self.geometry_int.A_o
        u_m_int = m_dot_int / (A_o_int * rho_int)  # m / s
        # determine mass flow rate in tubes
        A_tube = self._tube_int.Ac.to('m ** 2')
        m_dot_tube = rho_int * A_tube * u_m_int  # kg / s
        return m_dot_tube

    def _get_ext_mass_velocity(self) -> Quantity:
        """Returns the mass velocity of the external fluid."""
        A_o_ext = self.geometry_ext.A_o.to('m ** 2')
        m_dot_ext = self.m_dot_ext.to('kg / s')
        G_ext = m_dot_ext / A_o_ext
        return G_ext

    def _get_ext_reynolds_number(self, D: Quantity, fluid_ext: FluidState) -> float:
        """Returns the Reynolds number of the external fluid based on diameter
        `D` (either the hydraulic diameter, or the outside diameter of the
        tubes).
        """
        G_ext = self._get_ext_mass_velocity()
        mu_ext = fluid_ext.mu.to('Pa * s')
        D.ito('m')
        Re_ext = (G_ext * D / mu_ext).m
        return Re_ext

    def _get_int_heat_trf_coeff(self, fluid_int: FluidState) -> Quantity:
        """Returns heat transfer coefficient on internal side."""
        self._tube_int.m_dot = self._get_m_dot_tube(fluid_int)
        h_int = self._tube_int.avg_heat_transfer_coefficient(self.geometry_int.L1)
        # this assumes the same heat transfer coefficient in all tubes
        if isinstance(h_int, tuple):
            h_int = h_int[0]  # laminar flow with constant heat flux
        return h_int.to('W / (m ** 2 * K)')

    def _get_ext_heat_trf_coeff(self, fluid_ext: FluidState) -> Quantity:
        """Returns heat transfer coefficient on external side."""
        Re_D_ext = self._get_ext_reynolds_number(
            D=self.geometry_ext.D_o,
            fluid_ext=fluid_ext
        )
        Pr_ext = prandtl_number(
            rho=fluid_ext.rho,
            mu=fluid_ext.mu,
            k=fluid_ext.k,
            cp=fluid_ext.cp
        )
        Nu_ext = Nu_Briggs_and_Young(
            Re_D_o=Re_D_ext,
            Pr=Pr_ext,
            N_f=self.geometry_ext.N_f.to('1 / m').m,
            t_f=self.geometry_ext.t_f.to('m').m,
            h_f=(self.geometry_ext.D_f - self.geometry_ext.D_r).to('m').m / 2
        )
        h_ext = Nu_ext * fluid_ext.k / self.geometry_ext.D_h
        return h_ext.to('W / (m ** 2 * K)')

    @staticmethod
    def _get_int_resistance(h_int: Quantity, A_int: Quantity) -> Quantity:
        """Returns thermal resistance on internal side of heat exchanger."""
        h_int.ito('W / (m ** 2 * K)')
        A_int.ito('m ** 2')
        R_int = 1 / (h_int * A_int)  # K / W
        return R_int

    def _get_fin(self) -> Fin:
        fin = CircularRectangularFin(
            r_i=self.geometry_ext.D_r / 2,
            r_o=self.geometry_ext.D_f / 2,
            t=self.geometry_ext.t_f,
            k=self.k_f
        )
        return fin

    def _get_overall_fin_efficiency(self, h_ext: Quantity) -> float:
        """Returns overall efficiency of finned surface."""
        fin = self._get_fin()
        fin.h_avg = h_ext
        eta_fin = fin.efficiency
        A_f = self.geometry_ext.A_f  # m ** 2
        A_ext = self.geometry_ext.A  # m ** 2
        r = (A_f / A_ext).to('m ** 2 / m ** 2').m
        eta_ext = 1 - r * (1 - eta_fin)
        return eta_ext

    @staticmethod
    def _get_ext_resistance(
        h_ext: Quantity,
        eta_ext: float,
        A_ext: Quantity
    ) -> Quantity:
        """Returns thermal resistance on external side of heat exchanger."""
        h_ext.ito('W / (m ** 2 * K)')
        A_ext.ito('m ** 2')
        R_ext = 1 / (eta_ext * h_ext * A_ext)  # K / W
        return R_ext

    def _get_wall_temperature(
        self,
        fluid_int: FluidState,
        h_int: Quantity,
        fluid_ext: FluidState,
        h_ext: Quantity
    ) -> Quantity:
        """Returns the average wall temperature of the tubes."""
        R_int = self._get_int_resistance(
            h_int,
            self.geometry_int.A
        )
        eta_ext = self._get_overall_fin_efficiency(h_ext)
        R_ext = self._get_ext_resistance(
            h_ext,
            eta_ext,
            self.geometry_ext.A
        )
        T_int_bulk = fluid_int.T.to('K')
        T_ext_bulk = fluid_ext.T.to('K')
        n = T_int_bulk / R_int + T_ext_bulk / R_ext
        d = 1 / R_int + 1 / R_ext
        T_w = n / d
        return T_w.to('K')

    def _get_int_thermal_regime(self) -> str:
        """Returns thermal regime of the internal fluid."""
        if self.T_int_in > self.T_ext_in:
            # internal fluid is the hot fluid giving off heat to cold fluid on
            # external side.
            return "cooling"
        else:
            return "heating"

    @staticmethod
    def _get_nusselt_number(h: Quantity, D: Quantity, k: Quantity) -> float:
        h.ito('W / (m ** 2 * K)')
        D.ito('m')
        k.ito('W / (m * K)')
        Nu = (h * D / k).m
        return Nu

    def _correct_int_heat_trf_coeff(
        self,
        fluid_int: FluidState,
        h_int: Quantity,
        T_w: Quantity
    ) -> Quantity:
        """Correct internal heat transfer coefficient for variable fluid
        properties effects.
        """
        # determine Nusselt number from h_int
        Nu_cp = self._get_nusselt_number(h_int, self._tube_int.Dh, fluid_int.k)
        # correct Nusselt number
        Nu = correct_nusselt_number(
            Nu_cp=Nu_cp,
            T_w=T_w.to('K'),
            flow_regime=self._tube_int.get_flow_condition(),
            thermal_regime=self._get_int_thermal_regime(),
            fluid=fluid_int
        )
        # calculate corrected internal heat transfer coefficient
        h_int = Nu * fluid_int.k / self._tube_int.Dh
        return h_int.to('W / (m ** 2 * K)')

    def _get_ext_flow_regime(self, fluid_ext: FluidState) -> str:
        """Returns flow regime of external fluid."""
        Re_ext = self._get_ext_reynolds_number(
            self.geometry_ext.D_h,
            fluid_ext
        )
        if Re_ext > 2300:
            return 'turbulent'
        else:
            return 'laminar'

    def _get_ext_thermal_regime(self) -> str:
        """Returns thermal regime of external fluid."""
        if self.T_int_in > self.T_ext_in:
            # internal fluid is the hot fluid giving off heat to cold fluid on
            # external side.
            return "heating"
        else:
            return "cooling"

    def _correct_ext_heat_trf_coeff(
        self,
        fluid_ext: FluidState,
        h_ext: Quantity,
        T_w: Quantity
    ) -> Quantity:
        """Correct external heat transfer coefficient for variable fluid
        properties effects.
        """
        D_h = self.geometry_ext.D_h
        Nu_cp = self._get_nusselt_number(h_ext, D_h, fluid_ext.k)
        # correct Nusselt number
        Nu = correct_nusselt_number(
            Nu_cp=Nu_cp,
            T_w=T_w.to('K'),
            flow_regime=self._get_ext_flow_regime(fluid_ext),
            thermal_regime=self._get_ext_thermal_regime(),
            fluid=fluid_ext
        )
        # calculate corrected internal heat transfer coefficient
        h_ext = Nu * fluid_ext.k / D_h
        return h_ext.to('W / (m ** 2 * K)')

    def _get_heat_transfer_coeffs(
        self,
        fluid_int: FluidState,
        fluid_ext: FluidState
    ) -> tuple[Quantity, Quantity]:
        """Returns internal and external heat transfer coefficients."""
        h_int = self._get_int_heat_trf_coeff(fluid_int)
        h_ext = self._get_ext_heat_trf_coeff(fluid_ext)
        # correct heat transfer coefficients for variable fluid properties effects
        T_w = self._get_wall_temperature(fluid_int, h_int, fluid_ext, h_ext)
        h_int = self._correct_int_heat_trf_coeff(fluid_int, h_int, T_w)
        h_ext = self._correct_ext_heat_trf_coeff(fluid_ext, h_ext, T_w)
        return h_int, h_ext

    def _get_overall_conductance(
        self,
        h_int: Quantity,
        h_ext: Quantity
    ) -> Quantity:
        """Returns overall heat transfer conductance of the heat exchanger."""
        A_int = self.geometry_int.A
        A_ext = self.geometry_ext.A
        # calculate finned surface efficiency
        eta_ext = self._get_overall_fin_efficiency(h_ext)
        # calculate overall thermal conductance
        R_int = self._get_int_resistance(h_int, A_int)
        R_foul_int = Q_(0.0, 'K / W')  # we ignore fouling on internal side
        R_cond = Q_(0.0, 'K / W')  # we ignore the resistance to conduction
        R_ext = self._get_ext_resistance(h_ext, eta_ext, A_ext)
        R_foul_ext = Q_(0.0, 'K / W')  # we ignore fouling on external side
        R_tot = R_int + R_foul_int + R_cond + R_ext + R_foul_ext
        UA = 1 / R_tot
        return UA

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
            crossflow_hex = CrossFlowHeatExchanger(
                type_=CrossFlowHeatExchanger.Type.C_min_MIXED_C_max_UNMIXED,
                C_cold=self.m_dot_ext * fluid_ext.cp,
                C_hot=self.m_dot_int * fluid_int.cp,
                T_cold_in=self.T_ext_in,
                T_hot_in=self.T_int_in,
                UA=UA
            )
            T_int_out = crossflow_hex.T_hot_out
            T_ext_out = crossflow_hex.T_cold_out
        else:
            # internal flow is cold fluid
            crossflow_hex = CrossFlowHeatExchanger(
                type_=CrossFlowHeatExchanger.Type.C_min_MIXED_C_max_UNMIXED,
                C_cold=self.m_dot_int * fluid_int.cp,
                C_hot=self.m_dot_ext * fluid_ext.cp,
                T_cold_in=self.T_int_in,
                T_hot_in=self.T_ext_in,
                UA=UA
            )
            T_int_out = crossflow_hex.T_cold_out
            T_ext_out = crossflow_hex.T_hot_out
        return (
            T_int_out,
            T_ext_out,
            crossflow_hex.Q,
            crossflow_hex.eps,
            crossflow_hex.UA
        )

    def _get_int_flow_length(self) -> Quantity:
        """Returns the internal flow path length (pressure drop length) in case
        of a crossflow heat exchanger."""
        return self.geometry_int.L1

    def _get_int_pressure_drop(
        self,
        fluid_int: FluidState,
        h_int: Quantity,
        fluid_ext: FluidState,
        h_ext: Quantity
    ) -> Quantity:
        """Returns pressure drop on internal side."""
        f_cp = self._tube_int.friction_factor() / 4  # Fanning friction factor
        # correct friction factor for variable fluid properties effects
        f_corr = correct_friction_factor(
            f_cp=f_cp,
            T_w=self._get_wall_temperature(fluid_int, h_int, fluid_ext, h_ext),
            flow_regime=self._tube_int.get_flow_condition(),
            thermal_regime=self._get_int_thermal_regime(),
            fluid=fluid_int
        )
        # calculate pressure drop
        L_dP = self._get_int_flow_length()
        D_h = self._tube_int.Dh
        u_m = self._tube_int.mean_velocity()
        rho = fluid_int.rho
        dP_int = f_corr * (L_dP / D_h) * rho * u_m ** 2 / 2
        return dP_int.to('Pa')

    def _get_fanning_friction_factor(
        self,
        fluid_int: FluidState,
        h_int: Quantity,
        fluid_ext: FluidState,
        h_ext: Quantity,
        D_h: Quantity | None = None
    ) -> float:
        # calculate tube bank friction factor
        f_tb = f_Chai(
            Re_D_o=self._get_ext_reynolds_number(
                Q_(self.geometry_ext.D_o, 'm'),
                fluid_ext
            ),
            S_t=self.geometry_ext.S_t.to('m').m,
            S_l=self.geometry_ext.S_l.to('m').m,
            D_r=self.geometry_ext.D_r.to('m').m,
            D_f=self.geometry_ext.D_f.to('m').m,
            t_f=self.geometry_ext.t_f.to('m').m,
            N_f=self.geometry_ext.N_f.to('1 / m').m
        )
        # convert tube bank friction factor to Fanning friction factor
        f = convert_friction_factor(
            f_tb,
            self.geometry_ext.S_l.to('m').m,
            D_h.to('m').m
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

    def _get_ext_pressure_drop(
        self,
        fluid_int: FluidState,
        h_int: Quantity,
        fluid_ext: FluidState,
        h_ext: Quantity,
        T_ext_out: Quantity
    ) -> Quantity:
        """Returns pressure drop on external side."""
        G_ext = self._get_ext_mass_velocity()  # kg / (s * m ** 2)
        D_h = self.geometry_ext.D_h.to('m')
        L = self.geometry_ext.L2.to('m')
        sigma = self.geometry_ext.sigma
        rho_ext_in = self.Fluid_ext(
            T=self.T_ext_in,
            P=self.P_ext_in
        ).rho.to('kg / m ** 3')
        P_ext_out = self.P_ext_in  # initial guess
        tol_dP_ext = Q_(1.0, 'Pa')
        i_max = 5
        i = 0
        while i <= i_max:
            rho_ext_out = self.Fluid_ext(
                T=T_ext_out,
                P=P_ext_out
            ).rho.to('kg / m ** 3')
            rho_ext_avg = 2 * rho_ext_in * rho_ext_out / (rho_ext_in + rho_ext_out)
            f_corr = self._get_fanning_friction_factor(
                fluid_int, h_int,
                fluid_ext, h_ext,
                D_h
            )
            k = (G_ext ** 2) / (2 * rho_ext_in)
            a = f_corr * (4 * L / D_h) * (rho_ext_in / rho_ext_avg)
            b = (1 + sigma ** 2) * (rho_ext_in / rho_ext_out - 1)
            dP_ext = k * (a + b)
            P_ext_out_new = self.P_ext_in - dP_ext
            dev_dP_ext = abs(P_ext_out_new - P_ext_out)
            if dev_dP_ext.to('Pa') <= tol_dP_ext:
                return dP_ext
            P_ext_out = P_ext_out_new
            i += 1
        else:
            raise ValueError('no solution for external pressure drop')

    def _get_face_velocity(self) -> Quantity:
        """Get face velocity at external side inlet of heat exchanger."""
        m_dot_ext = self.m_dot_ext.to('kg / s')
        rho_ext_in = self.Fluid_ext(
            T=self.T_ext_in,
            P=self.P_ext_in
        ).rho.to('kg / m ** 3')
        v_fa = m_dot_ext / (rho_ext_in * self.geometry_ext.A_fr)
        return Q_(v_fa, 'm / s')


class CircularFinTubeCounterFlowHeatExchanger(CircularFinTubeCrossFlowHeatExchanger):

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
        """Returns the internal flow path length (pressure drop length) in case
        of a counterflow heat exchanger."""
        return self.geometry_int.L2
