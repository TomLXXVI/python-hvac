import numpy as np
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid, FluidState, CoolPropError
from hvac.heat_transfer.heat_exchanger.eps_ntu import CounterFlowHeatExchanger
from hvac.heat_transfer.heat_exchanger.fin_tube import core

from hvac.logging import ModuleLogger
logger = ModuleLogger.get_logger(__name__)

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
        rfg_out: FluidState,
        air_out: HumidAir,
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
            h_rfg_mean = (self.rfg_sat_liq_in.h + rfg_out.h) / 2
            h_air_mean = (self.air_in.h + air_out.h) / 2
            rfg_mean = self.Rfg(h=h_rfg_mean, P=self.rfg_sat_liq_in.P)  # we assume constant condenser pressure
            air_mean = HumidAir(h=h_air_mean, W=self.air_in.W)
            return rfg_mean, air_mean
        else:
            Dh = (
                self.rfg_sat_liq_in.h - air_out.h,
                rfg_out.h - self.air_in.h
            )
            Dh_max = max(Dh)
            Dh_min = min(Dh)
            if Dh_min < 0:
                logger.debug('Dh_min < 0 -> set to 1e-12 J/kg')
                Dh_min = Q_(1e-12, 'J / kg')
            LMED = (Dh_max - Dh_min) / np.log(Dh_max / Dh_min)
            if C_max == C_rfg:
                # Refrigerant has the smallest temperature change.
                h_rfg_mean = (self.rfg_sat_liq_in.h + rfg_out.h) / 2
                h_air_mean = h_rfg_mean - LMED
            else:  # C_max == C_ext
                # Air has the smallest temperature change.
                h_air_mean = (self.air_in.h + air_out.h) / 2
                h_rfg_mean = h_air_mean + LMED
            try:
                rfg_mean = self.Rfg(h=h_rfg_mean, P=self.rfg_sat_liq_in.P)
                # if T_rfg_mean is very near to saturation temperature, CoolProp
                # throws an exception to indicate that the full state cannot be
                # determined.
            except CoolPropError:
                rfg_mean = self.rfg_sat_liq_in
            air_mean = HumidAir(h=h_air_mean, W=self.air_in.W)
            return rfg_mean, air_mean

    def solve(self, L2: Quantity, i_max: int = 100, tol: float = 0.01) -> None:
        """Solves for the heat transfer performance of the subcooling part of
        the condenser for a given flow length L2 of the subcooling part.

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
        cp_air = self.air_in.cp
        cp_rfg = self.rfg_sat_liq_in.cp
        # Calculate capacitance rates:
        C_air = cp_air * self.hex_core.m_dot_ext
        C_rfg = cp_rfg * self.hex_core.m_dot_int
        C_min = min(C_air, C_rfg)
        # Guess initial value for the heat transfer effectiveness:
        self.eps = 0.25
        # Calculate heat transfer rate:
        Q_max = C_min * (self.rfg_sat_liq_in.T - self.air_in.Tdb)
        self.Q_dot = self.eps * Q_max
        # Calculate state of leaving air and leaving refrigerant:
        h_air_out = self.air_in.h + self.Q_dot / self.hex_core.ext.m_dot
        h_rfg_out = self.rfg_sat_liq_in.h - self.Q_dot / self.hex_core.int.m_dot
        self.air_out = HumidAir(h=h_air_out, W=self.air_in.W)
        self.rfg_out = self.Rfg(h=h_rfg_out, P=self.P_rfg)
        for i in range(i_max):
            # Calculate mean state of air and refrigerant between inlet and
            # outlet of subcooling region:
            rfg_mean, air_mean = self._get_mean_fluid_states(
                self.rfg_out, self.air_out,
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
            # Get new value of heat transfer effectiveness:
            eps_new = cof_hex.eps
            # Get new value of heat transfer rate:
            self.Q_dot = cof_hex.Q
            # Calculate state of leaving air and leaving refrigerant:
            self.air_out = HumidAir(
                Tdb=cof_hex.T_cold_out,
                W=self.air_in.W
            )
            try:
                self.rfg_out = self.Rfg(
                    T=cof_hex.T_hot_out,
                    P=self.P_rfg
                )
            except CoolPropError:
                # if rfg_out is very near to saturation, CoolProp
                # throws an exception to indicate that the full state cannot be
                # determined.
                self.rfg_out = self.rfg_sat_liq_in
            # Recalculate capacitance rates:
            dT_air = self.air_out.Tdb - self.air_in.Tdb
            dh_air = self.air_out.h - self.air_in.h
            try:
                cp_air_avg = dh_air / dT_air
            except ZeroDivisionError:
                cp_air_avg = self.air_in.cp
            dT_rfg = self.rfg_sat_liq_in.T - self.rfg_out.T
            dh_rfg = self.rfg_sat_liq_in.h - self.rfg_out.h
            try:
                cp_rfg_avg = dh_rfg / dT_rfg
            except ZeroDivisionError:
                cp_rfg_avg = self.rfg_sat_liq_in.cp
            C_air = self.hex_core.ext.m_dot * cp_air_avg
            C_rfg = self.hex_core.int.m_dot * cp_rfg_avg
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
