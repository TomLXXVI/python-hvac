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
        self._hex_core = HexCore(
            L1, L3,
            S_t, S_l,
            D_i, D_o,
            t_f, N_f, k_f
        )
        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.rfg_sat_liq_in: FluidState | None = None
        self.rfg_out: FluidState | None = None
        self.Q: Quantity | None = None
        self.eps: float | None = None
        self.Rfg: Fluid | None = None
        self.P_rfg: Quantity | None = None
        self.dT_sco: Quantity | None = None

    def set_fixed_operating_conditions(
        self,
        air_in: HumidAir,
        m_dot_air: Quantity,
        Rfg: Fluid,
        P_rfg: Quantity,
        m_dot_rfg: Quantity
    ) -> None:
        """Sets the fixed operating conditions on the subcooling region of
        the condenser.

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
        self._hex_core.m_dot_ext = m_dot_air
        self._hex_core.m_dot_int = m_dot_rfg
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
                # refrigerant has the smallest temperature change
                T_rfg_mean = (self.rfg_sat_liq_in.T.to('K') + T_rfg_out.to('K')) / 2
                DT_max = T_rfg_mean - self.air_in.Tdb.to('K')
                DT_min = max(Q_(1.e-12, 'K'), T_rfg_mean - T_air_out.to('K'))
                LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
                T_air_mean = T_rfg_mean - LMTD
                rfg_mean = self.Rfg(T=T_rfg_mean, P=self.rfg_sat_liq_in.P)
                air_mean = HumidAir(Tdb=T_air_mean, W=self.air_in.W)
                return rfg_mean, air_mean
            else:  # C_max == C_ext
                # air has the smallest temperature change
                T_air_mean = (self.air_in.Tdb.to('K') + T_air_out.to('K')) / 2
                DT_max = self.rfg_sat_liq_in.T.to('K') - T_air_mean
                DT_min = max(Q_(1.e-12, 'K'), T_rfg_out.to('K') - T_air_mean)
                LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
                T_rfg_mean = T_air_mean + LMTD
                rfg_mean = self.Rfg(T=T_rfg_mean, P=self.rfg_sat_liq_in.P)
                air_mean = HumidAir(Tdb=T_air_mean, W=self.air_in.W)
                return rfg_mean, air_mean

    def solve(self, L2: Quantity, i_max: int = 10, tol_eps: float = 0.01):
        """Determines by iteration the heat transfer performance of the
        subcooling part of the condenser for a given flow length L2 of the
        subcooling part.

        Returns
        -------
        rfg_out:
            State of refrigerant at outlet of condenser.
        air_out:
            State of air at outlet of subcooling part = inlet of condensing
            part.
        Q:
            Heat transfer rate in subcooling part.
        eps:
            Heat transfer effectiveness of subcooling part.
        """
        self._hex_core.L2 = L2
        # guess specific heats of air and refrigerant
        cp_air = CP_HUMID_AIR
        cp_rfg = self.rfg_sat_liq_in.cp
        # calculate capacitance rates
        C_air = cp_air * self._hex_core.m_dot_ext
        C_rfg = cp_rfg * self._hex_core.m_dot_int
        C_min = min(C_air, C_rfg)
        # guess a value of the heat transfer effectiveness
        eps = 0.75
        # calculate heat transfer rate and outlet states of air and refrigerant
        Q = eps * C_min * (self.rfg_sat_liq_in.T - self.air_in.Tdb)
        T_air_out = min(self.air_in.Tdb + Q / C_air, self.rfg_sat_liq_in.T)
        T_rfg_out = max(self.rfg_sat_liq_in.T - Q / C_rfg, self.air_in.Tdb)
        i = 0
        while i < i_max:
            # calculate the mean state of air and refrigerant between inlet and
            # outlet of subcooling region
            rfg_mean, air_mean = self._get_mean_fluid_states(T_rfg_out, T_air_out, C_rfg, C_air)
            # update heat exchanger core with mean states for calculating overall
            # conductance
            self._hex_core.int.fluid_mean = rfg_mean
            self._hex_core.ext.fluid_mean = air_mean
            # solve the heat exchanger equations
            cof_hex = CounterFlowHeatExchanger(
                C_cold=C_air,
                C_hot=C_rfg,
                T_cold_in=self.air_in.Tdb,
                T_hot_in=self.rfg_sat_liq_in.T,
                UA=self._hex_core.UA
            )
            # get new value for effectiveness and heat transfer rate
            eps_new = cof_hex.eps
            Q = cof_hex.Q
            # get new values for the outlet state of air and refrigerant
            T_air_out = cof_hex.T_cold_out
            T_rfg_out = cof_hex.T_hot_out
            # recalculate the capacitance rates
            C_air = Q / (T_air_out - self.air_in.Tdb)
            C_rfg = Q / (self.rfg_sat_liq_in.T - T_rfg_out)
            # check if the change of effectiveness has become small enough
            dev_eps = abs(eps_new - eps)
            if dev_eps <= tol_eps:
                self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
                self.rfg_out = self.Rfg(T=T_rfg_out, P=self.P_rfg)
                self.dT_sco = self.rfg_sat_liq_in.T - T_rfg_out
                self.Q = Q
                self.eps = eps
                return self.rfg_out, self.air_out, self.Q, self.eps
            # if not, do next iteration
            eps = eps_new
            i += 1
        else:
            raise ValueError(
                f'no acceptable solution found after {i_max} iterations'
            )
