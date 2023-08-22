import numpy as np
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState, Fluid
from hvac.heat_transfer.heat_exchanger.eps_ntu import CounterFlowHeatExchanger
from hvac.heat_transfer.heat_exchanger.fin_tube import core

from hvac.logging import ModuleLogger
logger = ModuleLogger.get_logger(__name__)

Q_ = Quantity
HexCore = core.plain_fin_tube.PlainFinTubeHeatExchangerCore


class PlainFinTubeCounterFlowSuperheatEvaporator:
    """Model of a single-pass plain fin tube counterflow heat exchanger with
    humid air flowing on the external side and with superheated refrigerant on the
    internal side. The states of the refrigerant at the evaporator's inlet and
    exit are fixed: at the inlet the refrigerant is saturated vapor and at the
    outlet it is superheated to the degree of superheat set on the expansion
    device.
    """
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
        """Creates a `PlainFinTubeCounterFlowSuperheatEvaporator` instance by
        defining the known dimensions of the heat exchanger core and its
        geometrical characteristics. The frontal area dimensions L1 and L3
        are known. The aim of this class is to determine the flow length L2 for
        which the refrigerant will be superheated to the predetermined degree
        of superheat.

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
        self.rfg_sat_in: FluidState | None = None
        self.Rfg: Fluid | None = None
        self._dT_rfg_sh: Quantity | None = None
        self.rfg_out: FluidState | None = None
        self.Q_dot: Quantity | None = None
        self.air_out: HumidAir | None = None

    def set_fixed_operating_conditions(
        self,
        air_in: HumidAir,
        m_dot_air: Quantity,
        Rfg: Fluid,
        P_rfg_sat: Quantity,
        dT_rfg_sh: Quantity
    ) -> None:
        """Sets the known operating conditions of the evaporator with refrigerant
        being superheated inside.

        Parameters
        ----------
        air_in:
            State of humid air at the evaporator's air inlet.
        m_dot_air:
            Mass flow rate of humid air (referred to mass of dry-air fraction).
        Rfg:
            Type of refrigerant.
        P_rfg_sat:
            Refrigerant saturation pressure (= evaporation pressure).
        dT_rfg_sh:
            Degree of superheat (set on expansion device).

        Notes
        -----
        The states of refrigerant at the inlet and outlet are known. At the
        inlet it is saturated vapor. Only the saturation temperature must be
        specified to fix the state of saturated vapor (vapor quality = 100 %).
        At the outlet the state of the refrigerant is fixed by the given degree
        of superheat, assuming that the evaporation pressure remains constant
        throughout the evaporator.
        """
        self.Q_dot = None
        self.air_out = None
        self.air_in = air_in
        self.hex_core.m_dot_ext = m_dot_air.to('kg / s')
        self.Rfg = Rfg
        # Refrigerant state at the inlet of the superheated region = saturated vapor
        self.rfg_sat_in = self.Rfg(P=P_rfg_sat, x=Q_(1.0, 'frac'))
        self._dT_rfg_sh = dT_rfg_sh.to('K')
        # Determine refrigerant outlet state with known degree of superheat:
        self.rfg_out = self.Rfg(
            T=self.rfg_sat_in.T.to('K') + self._dT_rfg_sh,
            P=self.rfg_sat_in.P  # ignore any pressure drop on refrigerant side
        )

    def set_mass_flow_rate_refrigerant(self, m_dot_rfg: Quantity) -> None:
        """Sets the provisional mass flow rate of refrigerant through the
        superheated region of the evaporator.
        """
        self.hex_core.m_dot_int = m_dot_rfg.to('kg / s')
        # Determine heat transfer rate in superheated region of evaporator:
        self.Q_dot = self.hex_core.m_dot_int * (self.rfg_out.h - self.rfg_sat_in.h)
        # Determine air outlet state:
        h_a_out = self.air_in.h - self.Q_dot / self.hex_core.m_dot_ext
        W_a_out = self.air_in.W  # assume heat transfer is only sensible: humidity ratio = constant
        self.air_out = HumidAir(h=h_a_out, W=W_a_out)

    def determine_flow_length(self, L2_ini: Quantity) -> Quantity:
        """Determines the flow length that is needed to superheat the
        refrigerant from saturated vapor to the given degree of superheat,
        knowing the mass flow rate of refrigerant.
        """
        # Finds the flow length for which the heat transfer rate equals the
        # heat transfer rate that follows from the fixed refrigerant inlet
        # and outlet states and the mass flow rate of refrigerant (see
        # method `set_operating_conditions`).
        def eq(L2: float) -> float:
            self.hex_core.L2 = Q_(L2, 'm')
            cof_hex = CounterFlowHeatExchanger(
                C_cold=rfg_mean.cp * self.hex_core.m_dot_int,
                C_hot=air_mean.cp * self.hex_core.m_dot_ext,
                T_cold_in=self.rfg_sat_in.T,
                T_hot_in=self.air_in.Tdb,
                UA=self.hex_core.UA
            )
            dQ = cof_hex.Q - self.Q_dot
            return dQ.to('W').m

        air_mean, rfg_mean = self._determine_mean_fluid_states()
        self.hex_core.ext.fluid_mean = air_mean
        self.hex_core.int.fluid_mean = rfg_mean
        try:
            L2 = root_scalar(eq, bracket=[1.e-6, L2_ini.to('m').m]).root
        except ValueError:
            raise ValueError(
                "Unable to reach set degree of superheat."
            )
        return Q_(L2, 'm')

    def _determine_mean_fluid_states(self) -> tuple[HumidAir, FluidState]:
        # Determine specific heat of air:
        air_avg = HumidAir(
            Tdb=(self.air_in.Tdb.to('K') + self.air_out.Tdb.to('K')) / 2,
            W=self.air_in.W
        )
        cp_air = air_avg.cp
        # Determine specific heat of refrigerant:
        rfg_avg = self.Rfg(
            T=(self.rfg_sat_in.T + self.rfg_out.T) / 2,
            P=self.rfg_sat_in.P
        )
        cp_rfg = rfg_avg.cp
        # Determine capacitance rates:
        C_air = cp_air * self.hex_core.m_dot_ext
        C_rfg = cp_rfg * self.hex_core.m_dot_int
        C_max = max(C_air, C_rfg)
        C_min = min(C_air, C_rfg)
        C_r = C_min / C_max
        # Determine mean air and mean refrigerant states:
        if C_r >= 0.5:
            T_air_m = (self.air_in.Tdb.to('K') + self.air_out.Tdb.to('K')) / 2
            T_rfg_m = (self.rfg_sat_in.T.to('K') + self.rfg_out.T.to('K')) / 2
            air_mean = HumidAir(Tdb=T_air_m, W=self.air_in.W)
            rfg_mean = self.Rfg(T=T_rfg_m, P=self.rfg_sat_in.P)
            return air_mean, rfg_mean
        else:
            DT = (
                self.air_in.Tdb.to('K') - self.rfg_out.T.to('K'),
                self.air_out.Tdb.to('K') - self.rfg_sat_in.T.to('K')
            )
            DT_max = max(DT)
            DT_min = min(DT)
            if DT_min < 0:
                logger.debug('DT_min < 0 -> set to 1e-12 K')
                DT_min = Q_(1e-12, 'K')
            LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
            if C_max == C_rfg:
                # Refrigerant has the smallest temperature change.
                T_rfg_m = (self.rfg_sat_in.T.to('K') + self.rfg_out.T.to('K')) / 2
                T_air_m = T_rfg_m + LMTD
            else:  # C_max == C_air
                # Air has the smallest temperature change.
                T_air_m = (self.air_in.Tdb.to('K') + self.air_out.Tdb.to('K')) / 2
                T_rfg_m = T_air_m - LMTD
            air_mean = HumidAir(Tdb=T_air_m, W=self.air_in.W)
            rfg_mean = self.Rfg(T=T_rfg_m, P=self.rfg_sat_in.P)
            return air_mean, rfg_mean

    @property
    def dP_air(self) -> Quantity:
        """Returns air-side pressure drop across superheating region of
        evaporator.
        """
        dP_air = self.hex_core.ext.get_pressure_drop(self.air_in, self.air_out)
        return dP_air
