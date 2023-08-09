import numpy as np
from scipy import optimize
from hvac import Quantity
from hvac.fluids import Fluid, FluidState, HumidAir, CP_HUMID_AIR, CoolPropError
from hvac.heat_transfer.heat_exchanger.fin_tube import core
from hvac.heat_transfer.heat_exchanger.eps_ntu import CounterFlowHeatExchanger

Q_ = Quantity
HexCore = core.plain_fin_tube.PlainFinTubeHeatExchangerCore


class PlainFinTubeCounterflowCondensingCondenser:
    """Model of the condensing part of the condenser."""

    def __init__(
        self,
        L1: Quantity,
        L3: Quantity,
        S_t: Quantity,
        S_l: Quantity,
        D_i: Quantity,
        D_o: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        k_f: Quantity = Q_(237, 'W / (m * K)')
    ) -> None:
        """Creates a `PlainFinTubeCounterFlowCondensingCondenser` instance by
        defining the known dimensions of the heat exchanger core and its
        geometrical characteristics.

        Parameters
        ----------
        L1:
            Tube length in the direction of inside flow, i.e. the length of the
            tube available for heat transfer with the outside flow.
        L3:
            Length of the tube bank perpendicular to the direction of external
            flow.
        S_t:
            Lateral or transverse pitch, i.e. distance between tubes of the
            same row.
        S_l:
            Longitudinal pitch, i.e. distance between tubes of two adjacent tube
            rows.
        D_i:
            Inside diameter of the tubes.
        D_o:
            Outside diameter of the tubes.
        t_f:
            Thickness of a circular fin.
        N_f:
            Number of fins per unit length (= fin density)
        k_f:
            Thermal conductivity of the external fin material.
        """
        # Create heat exchanger core:
        self._hex_core = HexCore(
            L1, L3,
            S_t, S_l,
            D_i, D_o,
            t_f, N_f, k_f,
            condensing=True
        )

        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.rfg_sat_vap_in: FluidState | None = None
        self.rfg_sat_liq_out: FluidState | None = None
        self.rfg_out: FluidState | None = None
        self.Q_dot: Quantity | None = None
        self.eps: Quantity | None = None
        self.dT_air: Quantity | None = None
        self.L2_condensing: Quantity | None = None
        self.Rfg: Fluid | None = None
        self.P_rfg: Quantity | None = None

    def set_fixed_operating_conditions(
        self,
        m_dot_air: Quantity,
        Rfg: Fluid,
        P_rfg: Quantity,
        m_dot_rfg: Quantity
    ) -> None:
        """Sets the known, fixed operating conditions on the condensing region
        of the condenser.

        Parameters
        ----------
        m_dot_air:
            Mass flow rate of air.
        Rfg:
            Type of refrigerant.
        P_rfg:
            Condenser pressure.
        m_dot_rfg:
            Mass flow rate of refrigerant.

        Returns
        -------
        None
        """
        self.air_in = None
        self.air_out = None
        self.rfg_out = None
        self.eps = None
        self.L2_condensing = None
        self._hex_core.m_dot_ext = m_dot_air
        self._hex_core.m_dot_int = m_dot_rfg
        self.rfg_sat_vap_in = Rfg(P=P_rfg, x=Q_(1.0, 'frac'))
        self.rfg_sat_liq_out = Rfg(P=P_rfg, x=Q_(0.0, 'frac'))
        self.Q_dot = self._hex_core.m_dot_int * (self.rfg_sat_vap_in.h - self.rfg_sat_liq_out.h)
        self.dT_air = self.Q_dot / (self._hex_core.m_dot_ext * CP_HUMID_AIR)
        self.Rfg = Rfg
        self.P_rfg = P_rfg

    def set_air_in(self, air_in: HumidAir) -> None:
        """Sets the state of air leaving the subcooling part and entering the
        condensing part of the condenser. To be used after solving for
        the subcooling part first.
        """
        self.air_in = air_in
        T_air_out = self.air_in.Tdb + self.dT_air
        # When entering air state is set, the state of air leaving the
        # condensing part can also be determined:
        self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)

    def _get_mean_fluid_states(self) -> tuple[FluidState, HumidAir]:
        """Returns the mean states of refrigerant and air across the
        condensing region of the condenser.
        """
        T_rfg_mean = self.rfg_sat_vap_in.T.to('K')
        DT_max = T_rfg_mean - self.air_in.Tdb.to('K')
        DT_min = max(T_rfg_mean - self.air_out.Tdb.to('K'), Q_(1.e-12, 'K'))
        LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
        T_air_mean = T_rfg_mean - LMTD
        x_rfg = max(
            self.rfg_out.x if self.rfg_out is not None else Q_(0.5, 'frac'),
            Q_(1.e-12, 'frac')
        )
        try:
            rfg_mean = self.Rfg(T=T_rfg_mean, x=x_rfg)
        except CoolPropError:
            # if `self.Rfg` is a pseudo-pure fluid, `x` cannot be an input
            # variable
            x_target = x_rfg.to('frac').m

            def _eq(h: float) -> float:
                state = self.Rfg(
                    P=self.rfg_sat_liq_out.P,
                    h=Q_(h, 'kJ / kg')
                )
                x = state.x.to('frac').m
                return x - x_target

            h_ini = self.rfg_sat_liq_out.h.to('kJ / kg').m
            h_fin = self.rfg_sat_vap_in.h.to('kJ / kg').m
            sol = optimize.root_scalar(_eq, bracket=[h_ini, h_fin])
            h = Q_(sol.root, 'kJ / kg')
            rfg_mean = self.Rfg(
                P=self.rfg_sat_liq_out.P,
                h=h
            )
        air_mean = HumidAir(Tdb=T_air_mean, W=self.air_in.W)
        return rfg_mean, air_mean

    def determine_flow_length(self, L2_ini: Quantity) -> Quantity:
        """Finds the flow length needed to condense refrigerant from saturated
        vapor to saturated liquid.
        """
        def _get_new_flow_length(A_int: Quantity) -> Quantity:
            n = A_int / (np.pi * self._hex_core.int.geo.D_i * self._hex_core.int.geo.L1)
            n -= 0.5
            d = self._hex_core.int.geo.L3 / (self._hex_core.int.geo.S_t * self._hex_core.int.geo.S_l)
            d -= 1 / (2 * self._hex_core.int.geo.S_l)
            L2 = n / d
            return L2.to('mm')

        def _eq(unknowns):
            L2 = Q_(unknowns[0], 'mm')
            self._hex_core.L2 = L2
            h_int = self._hex_core.int.h
            T_wall = self._hex_core.T_wall
            R_int = (rfg_mean.T - T_wall) / self.Q_dot
            A_int = 1 / (R_int * h_int)
            L2_new = _get_new_flow_length(A_int.to('m ** 2'))
            dev = L2_new - L2
            return np.array([dev.to('mm').m])

        rfg_mean, air_mean = self._get_mean_fluid_states()
        self._hex_core.int.fluid_mean = rfg_mean
        self._hex_core.ext.fluid_mean = air_mean
        roots = optimize.fsolve(_eq, np.array([L2_ini.to('mm').m]), xtol=0.001)
        L2 = Q_(roots[0], 'mm')
        return L2

    def solve(self, L2: Quantity, i_max: int = 100, tol: float = 0.01) -> None:
        """Solves for the heat transfer performance of the condensing part
        of the condenser when the flow length of the condensing part is
        given.

        Parameters
        ----------
        L2:
            Flow length of condensing part.
        i_max:
            Maximum number of iterations.
        tol:
            Acceptable deviation between last and previous calculated value
            of the heat transfer effectiveness `eps` of the condensing part.

        Returns
        -------
        None
        """
        self._hex_core.L2 = L2
        # Guess specific heat of air:
        cp_air = CP_HUMID_AIR
        C_min = C_air = self._hex_core.m_dot_ext * cp_air
        # Guess an initial value for the heat transfer effectiveness:
        self.eps = 0.75
        # Calculate heat transfer rate:
        Q_max = C_min * (self.rfg_sat_vap_in.T - self.air_in.Tdb)
        self.Q_dot = self.eps * Q_max
        # Calculate state of air and refrigerant leaving condensing part:
        T_air_out = self.air_in.Tdb + self.Q_dot / C_air
        h_rfg_out = self.rfg_sat_vap_in.h - self.Q_dot / self._hex_core.m_dot_int
        if T_air_out > self.rfg_sat_vap_in.T:
            T_air_out = self.rfg_sat_vap_in.T
            self.Q_dot = C_air * (T_air_out - self.air_in.Tdb)
            self.eps = self.Q_dot / Q_max
            h_rfg_out = self.rfg_sat_vap_in.h - self.Q_dot / self._hex_core.m_dot_int
            self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
            self.rfg_out = self.Rfg(h=h_rfg_out, P=self.P_rfg)
        if h_rfg_out < self.rfg_sat_liq_out.h:
            h_rfg_out = self.rfg_sat_liq_out.h
            self.Q_dot = self._hex_core.m_dot_int * (self.rfg_sat_vap_in.h - h_rfg_out)
            self.eps = self.Q_dot / Q_max
            T_air_out = self.air_in.Tdb + self.Q_dot / C_air
            self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
            self.rfg_out = self.Rfg(h=h_rfg_out, P=self.P_rfg)
            return
        for i in range(i_max):
            # Calculate mean state of air and refrigerant between in- and
            # outlet of condensing region:
            rfg_mean, air_mean = self._get_mean_fluid_states()
            # Update heat exchanger core with mean states for calculating overall
            # conductance:
            self._hex_core.int.fluid_mean = rfg_mean
            self._hex_core.ext.fluid_mean = air_mean
            # Solve heat exchanger equation:
            cof_hex = CounterFlowHeatExchanger(
                C_cold=self._hex_core.m_dot_ext * CP_HUMID_AIR,
                C_hot=Q_(float('inf'), 'W / K'),
                T_cold_in=self.air_in.Tdb,
                T_hot_in=self.rfg_sat_vap_in.T,
                UA=self._hex_core.UA
            )
            # Get new value of heat transfer effectiveness:
            eps_new = cof_hex.eps
            # Calculate new value of heat transfer rate:
            self.Q_dot = eps_new * Q_max
            # Get new state of leaving air and refrigerant:
            T_air_out = cof_hex.T_cold_out
            h_rfg_out = self.rfg_sat_vap_in.h - self.Q_dot / self._hex_core.m_dot_int
            self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
            self.rfg_out = self.Rfg(h=h_rfg_out, P=self.P_rfg)
            # Check if deviation between new and previous value is small enough:
            if abs(eps_new - self.eps) <= tol:
                return
            # Otherwise, repeat loop with new value:
            self.eps = eps_new
        else:
            raise ValueError(
                f'No acceptable solution found after {i_max} iterations.'
            )

    @property
    def dP_air(self) -> Quantity:
        """Air-side pressure drop across condensing part of condenser."""
        dP_air = self._hex_core.ext.get_pressure_drop(self.air_in, self.air_out)
        return dP_air
