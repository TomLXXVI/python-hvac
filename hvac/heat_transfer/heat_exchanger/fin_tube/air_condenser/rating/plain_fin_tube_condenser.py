import warnings
from collections import namedtuple
from scipy import optimize
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
        'Q_dot',
        'eps',
        'dP_air',
        'dT_sc',
        'L2_desuperheating',
        'L2_condensing',
        'L2_subcooling'
    ]
)


class CondenserWarning(Warning):
    pass


class PlainFinTubeCounterFlowCondenser:
    """Model of a single-pass, counterflow air condenser with an equilateral
    triangular staggered tube array and plain flat fins.
    Aim of this model is to find the heat transfer performance when the state
    and mass flow rate of refrigerant and of air at both entries of the
    condenser are known.
    """
    def __init__(
        self,
        L1: Quantity, L3: Quantity, N_r: int,
        S_t: Quantity, S_l: Quantity,
        D_i: Quantity, D_o: Quantity,
        t_f: Quantity, N_f: Quantity,
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
        self.N_r = N_r
        self.L2 = N_r * S_l
        self._geo_dict = {
            'L1': L1, 'L3': L3,
            'S_t': S_t, 'S_l': S_l,
            'D_i': D_i, 'D_o': D_o,
            't_f': t_f, 'N_f': N_f,
            'k_f': k_f
        }
        # We can distinguish a subcooling part, a condensing part and a
        # desuperheating part of the condenser:
        self.subcooling_part = PlainFinTubeCounterFlowSubcoolingCondenser(**self._geo_dict)
        self.condensing_part = PlainFinTubeCounterflowCondensingCondenser(**self._geo_dict)
        self.desuperheating_part = PlainFinTubeCounterFlowDesuperheatCondenser(**self._geo_dict)

        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.rfg_in: FluidState | None = None
        self.rfg_out: FluidState | None = None
        self.m_dot_air: Quantity | None = None
        self.m_dot_rfg: Quantity | None = None
        self.Q_dot_dsp: Quantity | None = None
        self.Q_dot_cdp: Quantity | None = None
        self.Q_dot_scp: Quantity | None = None
        self.Q_dot: Quantity | None = None
        self.L2_dsp: Quantity | None = None
        self.L2_cdp: Quantity | None = None
        self.L2_scp: Quantity | None = None
        self._sub_cooling_flag: bool = True

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
        self.air_out = None
        self.rfg_out = None
        self.Q_dot_dsp = None
        self.Q_dot_cdp = None
        self.Q_dot_scp = None
        self.Q_dot = None
        self.L2_dsp = None
        self.L2_cdp = None
        self.L2_scp = None
        self._sub_cooling_flag = True
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

    def _fun_find_subcooling_length(
        self,
        L2_scp: float,
        i_max: int = 100,
        tol: float = 0.01
    ) -> float:
        # Equation with subcooling.
        # Solve subcooling part:
        L2_scp = Q_(L2_scp, 'mm')
        self.subcooling_part.solve(L2_scp, i_max, tol)
        self.rfg_out = self.subcooling_part.rfg_out
        # Set state of air entering the condensing part:
        T_cdp_air_in = self.subcooling_part.air_out.Tdb
        cdp_air_in = HumidAir(Tdb=T_cdp_air_in, W=self.air_in.W)
        self.condensing_part.set_air_in(cdp_air_in)
        # Set state of air entering the desuperheating part:
        T_dsp_air_in = T_cdp_air_in + self.condensing_part.dT_air
        dsp_air_in = HumidAir(Tdb=T_dsp_air_in, W=self.air_in.W)
        self.desuperheating_part.set_air_in(dsp_air_in)
        # Set state of air leaving the condenser:
        T_air_out = T_dsp_air_in + self.desuperheating_part.dT_air
        self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
        # Get total heat rejection rate in the condenser:
        self.Q_dot_dsp = self.desuperheating_part.Q_dot
        self.Q_dot_cdp = self.condensing_part.Q_dot
        self.Q_dot_scp = self.subcooling_part.Q_dot
        self.Q_dot = self.Q_dot_dsp + self.Q_dot_cdp + self.Q_dot_scp
        # Get flow lengths of subcooling, condensing and desuperheating
        # part:
        self.L2_scp = L2_scp
        self.L2_cdp = self.condensing_part.determine_flow_length(self.L2)
        self.L2_dsp = self.desuperheating_part.determine_flow_length(self.L2)
        # Determine total flow length and compare with actual flow length of
        # the condenser:
        L2 = self.L2_scp + self.L2_cdp + self.L2_dsp
        dev = L2 - self.L2
        return dev.to('mm').m

    def _fun_find_condensing_length(
        self,
        L_cdp: float,
        i_max: int = 100,
        tol: float = 0.01
    ) -> float:
        # Equation without subcooling.
        # Solve condensing part:
        L2_cdp = Q_(L_cdp, 'mm')
        self.condensing_part.set_air_in(self.air_in)
        self.condensing_part.solve(L2_cdp, i_max, tol)
        self.rfg_out = self.condensing_part.rfg_out
        # Set state of air entering the desuperheating part:
        T_dsp_air_in = self.condensing_part.air_out.Tdb
        dsp_air_in = HumidAir(Tdb=T_dsp_air_in, W=self.air_in.W)
        self.desuperheating_part.set_air_in(dsp_air_in)
        # Set state of air leaving the condenser:
        T_air_out = T_dsp_air_in + self.desuperheating_part.dT_air
        self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
        # Get total heat rejection rate in the condenser:
        self.Q_dot_dsp = self.desuperheating_part.Q_dot
        self.Q_dot_cdp = self.condensing_part.Q_dot
        self.Q_dot_scp = Q_(0.0, 'W')
        self.Q_dot = self.Q_dot_dsp + self.Q_dot_cdp + self.Q_dot_scp
        # Get flow lengths of subcooling, condensing and desuperheating
        # part:
        self.L2_scp = Q_(0.0, 'mm')
        self.L2_cdp = L2_cdp
        self.L2_dsp = self.desuperheating_part.determine_flow_length(self.L2)
        # Determine total flow length and compare with actual flow length of
        # the condenser:
        L2 = self.L2_scp + self.L2_cdp + self.L2_dsp
        dev = L2 - self.L2
        return dev.to('mm').m

    def rate(
        self,
        i_max: int = 100,
        tol: float = 0.01
    ) -> Result | None:
        """Determines the heat transfer performance of the condenser.

        Parameters
        ----------
        i_max:
            Maximum number of iterations.
        tol:
            Tolerance for the effectiveness of the subcooling part of the
            condenser.

        Returns
        -------
        namedtuple `Result`:
            rfg_out: FluidState
                State of leaving refrigerant.
            air_out: HumidAir
                State of leaving air.
            Q: Quantity
                Total heat rejection rate of the condenser.
            eps: Quantity
                Heat transfer effectiveness of the condenser.
            dP_air: Quantity
                Air-side pressure drop across condenser.
            dT_sc: Quantity
                Degree of subcooling.
            L2_desuperheating: Quantity
                Flow length along the desuperheating region of the condenser.
            L2_condensing: Quantity
                Flow length along the condensing region of the condenser.
            L2_subcooling: Quantity
                Flow length along the subcooling region of the condenser.
        """
        try:
            # Try to find the subcooling flow length for which the calculated
            # total condenser flow length equals the actual total flow length of
            # the condenser:
            optimize.root_scalar(
                self._fun_find_subcooling_length,
                # method='secant',
                # x0=0.50 * self.L2.to('mm').m,  # initial guess for `L_scp`
                bracket=(1e-9, self.L2.to('mm').m),
                xtol=0.1,  # mm
                maxiter=i_max
            )
        except Exception as err:
            warnings.warn(
                f"No subcooling of refrigerant in condenser: {err}",
                category=CondenserWarning
            )
            self._sub_cooling_flag = False
            # If no subcooling is possible at current operating conditions, try
            # to find the condensing flow length for which the calculated total
            # condenser flow length equals the actual flow length of the
            # condenser:
            try:
                optimize.root_scalar(
                    self._fun_find_condensing_length,
                    # method='secant',
                    # x0=0.50 * self.L2.to('mm').m,  # initial guess for `L_cdp`
                    bracket=(1e-9, self.L2.to('mm').m),
                    xtol=0.1,  # mm
                    maxiter=i_max
                )
            except Exception as err:
                raise ValueError(
                    f'Rating of condenser failed with error: "{err}"'
                ) from None

        return Result(
            rfg_out=self.rfg_out,
            air_out=self.air_out,
            Q_dot=self.Q_dot,
            eps=self.eps,
            dP_air=self.dP_air,
            dT_sc=self.dT_sc,
            L2_desuperheating=self.L2_dsp,
            L2_condensing=self.L2_cdp,
            L2_subcooling=self.L2_scp,
        )

    @property
    def eps(self) -> Quantity:
        """Heat transfer effectiveness of the condenser."""
        T_cnd = self.condensing_part.rfg_sat_liq_out.T
        eps = self.Q_dot / (self.m_dot_air * CP_HUMID_AIR * (T_cnd - self.air_in.Tdb))
        return eps.to('frac')

    @property
    def dP_air(self) -> Quantity:
        """Air-side pressure drop across the condenser."""
        if self.L2_scp.m > 0.0:
            return (
                self.subcooling_part.dP_air +
                self.condensing_part.dP_air +
                self.desuperheating_part.dP_air
            )
        else:
            return (
                self.condensing_part.dP_air +
                self.desuperheating_part.dP_air
            )

    @property
    def dT_sc(self) -> Quantity:
        """Degree of subcooling of the refrigerant leaving the condenser."""
        if self._sub_cooling_flag:
            return self.subcooling_part.dT_sc
        else:
            return Q_(0.0, 'K')

    def __call__(
        self,
        m_dot_air: Quantity,
        m_dot_rfg: Quantity,
        air_in: HumidAir,
        rfg_in: FluidState,
        i_max: int = 100,
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
