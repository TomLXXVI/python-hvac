from collections import namedtuple
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import FluidState, HumidAir, CP_HUMID_AIR
from .plain_fin_tube_subcooling_condenser import PlainFinTubeCounterFlowSubcoolingCondenser
from .plain_fin_tube_condensing_condenser import PlainFinTubeCounterflowCondensingCondenser
from .plain_fin_tube_desuperheat_condenser import PlainFinTubeCounterFlowDesuperheatCondenser


Q_ = Quantity


Result = namedtuple(
    'Result',
    [
        'rfg_out',
        'air_out',
        'Q',
        'eps',
        'dP_air',
        'dT_sc',
        'L2_desuperheating',
        'L2_condensing',
        'L2_subcooling'
    ]
)


class PlainFinTubeCounterFlowCondenser:
    """Model of a single-pass, counterflow air condenser with an equilateral
    triangular staggered tube array and plain flat fins.
    Aim of this model is to find the heat transfer performance when the state
    and mass flow rate of refrigerant and of air at both entries of the condenser
    are known.
    """
    def __init__(
        self,
        L1: Quantity,
        L3: Quantity,
        N_r: int,
        S_t: Quantity,
        S_l: Quantity,
        D_i: Quantity,
        D_o: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        k_f: Quantity = Q_(237, 'W / (m * K)')
    ) -> None:
        """Creates a `PlainFinTubeCounterFlowCondenser` instance by
        defining the dimensions of the heat exchanger core and its geometrical
        characteristics.

        Parameters
        ----------
        L1:
            Width of the heat exchanger core.
        L3:
            Height of the heat exchanger core. The product of L1 and L3
            determines the frontal area for the air flow through the
            heat exchanger.
        N_r:
            Number of rows.
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
        # we distinguish a subcooling part, a condensing part and a desuperheating
        # part of the condenser.
        self.L2 = N_r * S_l
        kwargs = {
            'L1': L1, 'L3': L3,
            'S_t': S_t, 'S_l': S_l,
            'D_i': D_i, 'D_o': D_o,
            't_f': t_f, 'N_f': N_f,
            'k_f': k_f
        }
        self.subcooling_part = PlainFinTubeCounterFlowSubcoolingCondenser(**kwargs)
        self.condensing_part = PlainFinTubeCounterflowCondensingCondenser(**kwargs)
        self.desuperheating_part = PlainFinTubeCounterFlowDesuperheatCondenser(**kwargs)

        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.dP_air: Quantity | None = None

        self.rfg_in: FluidState | None = None
        self.rfg_out: FluidState | None = None

        self.m_dot_air: Quantity | None = None
        self.m_dot_rfg: Quantity | None = None

        self.Q_dsh: Quantity | None = None
        self.Q_cnd: Quantity | None = None
        self.Q_sco: Quantity | None = None
        self.Q: Quantity | None = None

        self.L2_dsh: Quantity | None = None
        self.L2_cnd: Quantity | None = None
        self.L2_sco: Quantity | None = None

    def set_operating_conditions(
        self,
        air_in: HumidAir,
        m_dot_air: Quantity,
        rfg_in: FluidState,
        m_dot_rfg: Quantity
    ) -> None:
        """Sets the operating conditions of the condenser.

        Parameters
        ----------
        air_in:
            State of air entering the condenser.
        m_dot_air:
            Mass flow rate of air through the condenser.
        rfg_in:
            State of refrigerant entering the condenser.
        m_dot_rfg:
            Mass flow rate of refrigerant through the condenser, which under
            steady-state conditions equals the mass flow rate through the
            evaporator.
        """
        self.m_dot_air = m_dot_air.to('kg / s')
        self.m_dot_rfg = m_dot_rfg.to('kg / s')
        self.air_in = air_in
        self.rfg_in = rfg_in
        self.desuperheating_part.set_fixed_operating_conditions(
            m_dot_air=self.m_dot_air,
            rfg_in=self.rfg_in,
            m_dot_rfg=self.m_dot_rfg
        )
        self.condensing_part.set_fixed_operating_conditions(
            m_dot_air=self.m_dot_air,
            Rfg=self.rfg_in.fluid,
            P_rfg=self.rfg_in.P,
            m_dot_rfg=self.m_dot_rfg
        )
        self.subcooling_part.set_fixed_operating_conditions(
            air_in=self.air_in,
            m_dot_air=self.m_dot_air,
            Rfg=self.rfg_in.fluid,
            P_rfg=self.rfg_in.P,
            m_dot_rfg=self.m_dot_rfg
        )

    def rate(
        self,
        i_max: int = 10,
        tol_eps: float = 0.01
    ) -> Result:
        """Determines by iteration for a given state of entering air and
        refrigerant, the state of the air and refrigerant leaving the condenser.

        Parameters
        ----------
        i_max:
            Maximum number of iterations.
        tol_eps:
            Tolerance for the effectiveness of the subcooling part of the
            condenser.

        Raises
        ------
        ValueError if the number of iterations reached the maximum number and
        no solution was found within the given tolerance.

        Returns
        -------
        namedtuple Result:
            rfg_out: FluidState
                State of leaving refrigerant.
            air_out: HumidAir
                State of leaving air.
            Q: Quantity
                Total heat rejection rate of the condenser.
            dT_sc: Quantity
                Degree of subcooling.
            L2_desuperheating: Quantity
                Flow length along the desuperheating region of the condenser.
            L2_condensing: Quantity
                Flow length along the condensing region of the condenser.
            L2_subcooling: Quantity
                Flow length along the subcooling region of the condenser.
        """
        # First we solve the subcooling part of the condenser using a guessed
        # value of the subcooling flow length. We know the refrigerant inlet
        # state to the subcooling part (saturated liquid). We also know the
        # inlet state of air. With this it is possible to determine the air
        # state at the inlet of the condensing part of the condenser and the
        # state of refrigerant leaving the condenser.
        # We know the inlet (saturated vapor) and outlet state (saturated liquid)
        # of refrigerant at the condensing part of the condenser. As mass flow
        # rates are given, we can immediately determine the heat transfer rate
        # and the temperature rise of the air in the condensing part (this is
        # already done when we called the `set_fixed_operation_conditions`
        # methods in the `set_operations_conditions` method of
        # `PlainFinTubeCounterFlowCondenser`).
        # For the desuperheating part, the same applies. We know the inlet
        # state of refrigerant to the condenser, and we know that the
        # desuperheating region ends where refrigerant has turned into saturated
        # vapor.
        # Knowing the air outlet state of the subcooling part (= air inlet
        # state of the condensing part) and the temperature rise across the
        # condensing part, we can determine the air outlet state of the
        # condensing part (= air inlet state of the desuperheating part).
        # Knowing also the temperature rise of air across the desuperheating
        # part, the state of air leaving the condenser can be determined.
        # Next, we can determine the flow length of the condensing part and
        # desuperheating part. The sum of the subcooling, condensing and
        # desuperheating flow lengths should be equal to the total flow length
        # of the condenser. So, we need to the find the subcooling flow length
        # for which the sum of flow lengths equals the total flow length of the
        # condenser.

        def eq(L2_subcooling: float) -> float:
            # solve the subcooling part
            L2_subcooling = Q_(L2_subcooling, 'm')
            self.subcooling_part.solve(L2_subcooling, i_max, tol_eps)
            self.rfg_out = self.subcooling_part.rfg_out
            # set the state of air entering the condensing part
            T_air_in_cnd = self.subcooling_part.air_out.Tdb
            air_in_cnd = HumidAir(Tdb=T_air_in_cnd, W=self.air_in.W)
            self.condensing_part.set_air_in(air_in_cnd)
            # set the state of air entering the desuperheating part
            T_air_in_dsh = T_air_in_cnd + self.condensing_part.dT_air
            air_in_dsh = HumidAir(Tdb=T_air_in_dsh, W=self.air_in.W)
            self.desuperheating_part.set_air_in(air_in_dsh)
            # if T_air_in_dsh > self.condensing_part.rfg_sat_vap_in.T:
            #     raise ValueError(
            #         "Air temperature exceeds condenser temperature."
            #         "Raise condenser pressure."
            #     )
            # set the state of air leaving the condenser
            T_air_out = T_air_in_dsh + self.desuperheating_part.dT_air
            self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
            # get the total heat rejection rate of the condenser
            self.Q_dsh = self.desuperheating_part.Q
            self.Q_cnd = self.condensing_part.Q
            self.Q_sco = self.subcooling_part.Q
            self.Q = self.Q_dsh + self.Q_cnd + self.Q_sco
            # retrieve the flow lengths of the subcooling, condensing and
            # desuperheating part
            self.L2_sco = L2_subcooling
            self.L2_cnd = self.condensing_part.determine_flow_length(L2_ini=self.L2)
            self.L2_dsh = self.desuperheating_part.determine_flow_length(L2_ini=self.L2)
            # determine the totale flow length and compare with the actual flow
            # length of the condenser
            L2_tot = self.L2_sco + self.L2_cnd + self.L2_dsh
            dev_L2 = self.L2 - L2_tot
            return dev_L2.to('m').m

        # Find the minimum subcooling flow length for which the subcooling part
        # can be solved.
        L2_subcooling_min = -1.0
        for i in range(1, 100):
            L2_subcooling_min = i / 100 * self.L2.to('m').m
            try:
                self.subcooling_part.solve(Q_(L2_subcooling_min, 'm'), i_max, tol_eps)
            except ValueError:
                continue
            else:
                break
        # If a minimal subcooling flow length exists, try to find the subcooling
        # flow length such that the calculated total flow length equals the
        # actual flow length of the condenser.
        if L2_subcooling_min > 0.0:
            try:
                root_scalar(eq, bracket=[L2_subcooling_min, self.L2.to('m').m])
            except ValueError:
                raise ValueError(
                    'condenser flow length too short to attain subcooling '
                    'of refrigerant'
                )
            else:
                self.dP_air = (
                    self.subcooling_part.dP_air
                    + self.condensing_part.dP_air
                    + self.desuperheating_part.dP_air
                )
                return Result(
                    rfg_out=self.rfg_out,
                    air_out=self.air_out,
                    Q=self.Q,
                    eps=self._get_eps(),
                    dP_air=self.dP_air,
                    dT_sc=self.dT_sc,
                    L2_desuperheating=self.L2_dsh,
                    L2_condensing=self.L2_cnd,
                    L2_subcooling=self.L2_sco,
                )

    def _get_eps(self) -> float:
        T_cnd = self.condensing_part.rfg_sat_liq_out.T
        eps = self.Q / (self.m_dot_air * CP_HUMID_AIR * (T_cnd - self.air_in.Tdb))
        return eps

    @property
    def dT_sc(self) -> Quantity:
        return self.subcooling_part.dT_sco

    def __call__(
        self,
        m_dot_air: Quantity,
        m_dot_rfg: Quantity,
        air_in: HumidAir,
        rfg_in: FluidState,
        i_max: int = 10,
        tol_eps: float = 0.01
    ) -> Result:
        """Combines the methods `set_operating_conditions` and `rate` in a
        single method.
        """
        self.set_operating_conditions(
            air_in, m_dot_air,
            rfg_in, m_dot_rfg
        )
        res = self.rate(i_max, tol_eps)
        return res
