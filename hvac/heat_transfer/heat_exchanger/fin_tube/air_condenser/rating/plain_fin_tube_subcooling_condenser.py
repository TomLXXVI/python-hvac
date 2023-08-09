import numpy as np
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid, FluidState, CP_HUMID_AIR
from hvac.heat_transfer.heat_exchanger.eps_ntu import CounterFlowHeatExchanger
from hvac.heat_transfer.heat_exchanger.fin_tube import core


Q_ = Quantity
HexCore = core.plain_fin_tube.PlainFinTubeHeatExchangerCore


class PlainFinTubeCounterFlowSubcoolingCondenser:

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
        """Creates a `PlainFinTubeCounterFlowSubcoolingCondenser` instance
        by defining the frontal area dimensions and geometrical characteristics
        of the condenser heat transfer surfaces.

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
        self.hex_core = HexCore(
            L1, L3,
            S_t, S_l,
            D_i, D_o,
            t_f, N_f, k_f
        )
        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.rfg_sat_liq_in: FluidState | None = None
        self.rfg_out: FluidState | None = None
        self.Q_dot: Quantity | None = None
        self.eps: float | None = None
        self.Rfg: Fluid | None = None
        self.P_rfg: Quantity | None = None

    def set_fixed_operating_conditions(
        self,
        air_in: HumidAir,
        m_dot_air: Quantity,
        Rfg: Fluid,
        P_rfg: Quantity,
        m_dot_rfg: Quantity
    ) -> None:
        """Sets the known, fixed operating conditions on the subcooling region
        of the condenser.

        Parameters
        ----------
        air_in:
            State of air entering the condenser.
        m_dot_air:
            Mass flow rate of air.
        Rfg:
            Type of refrigerant.
        P_rfg:
            Condenser pressure (assumed to be constant throughout the condenser).
        m_dot_rfg:
            Mass flow rate of refrigerant.

        Returns
        -------
        None
        """
        self.air_out = None
        self.rfg_out = None
        self.Q_dot = None
        self.eps = None
        self.hex_core.m_dot_ext = m_dot_air
        self.hex_core.m_dot_int = m_dot_rfg
        self.air_in = air_in
        self.rfg_sat_liq_in = Rfg(P=P_rfg, x=Q_(0.0, 'frac'))
        self.Rfg = Rfg
        self.P_rfg = P_rfg

    def _get_mean_fluid_states(
        self,
        T_rfg_out: Quantity,
        T_air_out: Quantity,
        C_rfg: Quantity,
        C_air: Quantity
    ) -> tuple[FluidState, HumidAir]:
        """Returns the mean states of refrigerant and air across the
        desuperheating region of the condenser.
        """
        C_min = min(C_air, C_rfg)
        C_max = max(C_air, C_rfg)
        C_r = C_min / C_max
        if C_r >= 0.5:
            T_rfg_mean = (self.rfg_sat_liq_in.T.to('K') + T_rfg_out.to('K')) / 2
            T_air_mean = (self.air_in.Tdb.to('K') + T_air_out.to('K')) / 2
            rfg_mean = self.Rfg(T=T_rfg_mean, P=self.rfg_sat_liq_in.P)  # we assume constant condenser pressure
            air_mean = HumidAir(Tdb=T_air_mean, W=self.air_in.W)
            return rfg_mean, air_mean
        else:
            if C_max == C_rfg:
                # Refrigerant has the smallest temperature change.
                T_rfg_mean = (self.rfg_sat_liq_in.T.to('K') + T_rfg_out.to('K')) / 2
                DT_max = T_rfg_mean - self.air_in.Tdb.to('K')
                DT_min = max(Q_(1.e-12, 'K'), T_rfg_mean - T_air_out.to('K'))
                LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
                T_air_mean = T_rfg_mean - LMTD
                rfg_mean = self.Rfg(T=T_rfg_mean, P=self.rfg_sat_liq_in.P)
                air_mean = HumidAir(Tdb=T_air_mean, W=self.air_in.W)
                return rfg_mean, air_mean
            else:  # C_max == C_ext
                # Air has the smallest temperature change.
                T_air_mean = (self.air_in.Tdb.to('K') + T_air_out.to('K')) / 2
                DT_max = self.rfg_sat_liq_in.T.to('K') - T_air_mean
                DT_min = max(Q_(1.e-12, 'K'), T_rfg_out.to('K') - T_air_mean)
                LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
                T_rfg_mean = T_air_mean + LMTD
                rfg_mean = self.Rfg(T=T_rfg_mean, P=self.rfg_sat_liq_in.P)
                air_mean = HumidAir(Tdb=T_air_mean, W=self.air_in.W)
                return rfg_mean, air_mean

    def solve(self, L2: Quantity, i_max: int = 100, tol: float = 0.01) -> None:
        """Solves for the heat transfer performance of the
        subcooling part of the condenser for a given flow length L2 of the
        subcooling part.

        Parameters
        ----------
        L2:
            Flow length of subcooling part.
        i_max:
            Maximum number of iterations.
        tol:
            Acceptable deviation between last and previous calculated value
            of the heat transfer effectiveness `eps` of the condensing part.

        Returns
        -------
        None
        """
        self.hex_core.L2 = L2
        # Guess specific heats of air and refrigerant:
        cp_air = CP_HUMID_AIR
        cp_rfg = self.rfg_sat_liq_in.cp
        # Calculate capacitance rates:
        C_air = cp_air * self.hex_core.m_dot_ext
        C_rfg = cp_rfg * self.hex_core.m_dot_int
        C_min = min(C_air, C_rfg)
        # Guess initial value for the heat transfer effectiveness:
        self.eps = 0.75
        # Calculate heat transfer rate:
        Q_max = C_min * (self.rfg_sat_liq_in.T - self.air_in.Tdb)
        self.Q_dot = self.eps * Q_max
        # Calculate state of leaving air and leaving refrigerant:
        T_air_out = self.air_in.Tdb + self.Q_dot / C_air
        T_rfg_out = self.rfg_sat_liq_in.T - self.Q_dot / C_rfg
        if T_air_out > self.rfg_sat_liq_in.T:
            T_air_out = self.rfg_sat_liq_in.T
            self.Q_dot = C_air * (T_air_out - self.air_in.Tdb)
            self.eps = self.Q_dot / Q_max
            T_rfg_out = self.rfg_sat_liq_in.T - self.Q_dot / C_rfg
        if T_rfg_out < self.air_in.Tdb:
            T_rfg_out = self.air_in.Tdb
            self.Q_dot = C_rfg * (self.rfg_sat_liq_in.T - T_rfg_out)
            T_air_out = self.air_in.Tdb + self.Q_dot / C_air
            self.eps = self.Q_dot / Q_max
        self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
        self.rfg_out = self.Rfg(T=T_rfg_out, P=self.P_rfg)
        for i in range(i_max):
            # Calculate mean state of air and refrigerant between inlet and
            # outlet of subcooling region:
            rfg_mean, air_mean = self._get_mean_fluid_states(
                T_rfg_out, T_air_out,
                C_rfg, C_air
            )
            # Update heat exchanger core with mean states for calculating
            # the overall heat transfer conductance:
            self.hex_core.int.fluid_mean = rfg_mean
            self.hex_core.ext.fluid_mean = air_mean
            # Solve heat exchanger equation:
            cof_hex = CounterFlowHeatExchanger(
                C_cold=C_air,
                C_hot=C_rfg,
                T_cold_in=self.air_in.Tdb,
                T_hot_in=self.rfg_sat_liq_in.T,
                UA=self.hex_core.UA
            )
            # Get new value for heat transfer effectiveness:
            eps_new = cof_hex.eps
            # Get new value for heat transfer rate:
            self.Q_dot = cof_hex.Q
            # Calculate state of leaving air and leaving refrigerant:
            T_air_out = cof_hex.T_cold_out
            T_rfg_out = cof_hex.T_hot_out
            self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
            self.rfg_out = self.Rfg(T=T_rfg_out, P=self.P_rfg)
            # Recalculate capacitance rates:
            C_air = self.Q_dot / (T_air_out - self.air_in.Tdb)
            C_rfg = self.Q_dot / (self.rfg_sat_liq_in.T - T_rfg_out)
            # Check if deviation between new and previous value is small enough:
            if abs(eps_new - self.eps) <= tol:
                return
            # otherwise, repeat loop with new value
            self.eps = eps_new
        else:
            raise ValueError(
                f'No acceptable solution found after {i_max} iterations.'
            )

    @property
    def dP_air(self) -> Quantity:
        """Air-side pressure drop across condensing part of condenser."""
        return self.hex_core.ext.get_pressure_drop(self.air_in, self.air_out)

    @property
    def dT_sc(self) -> Quantity:
        """Degree of subcooling of refrigerant leaving condenser."""
        return self.rfg_sat_liq_in.T - self.rfg_out.T
