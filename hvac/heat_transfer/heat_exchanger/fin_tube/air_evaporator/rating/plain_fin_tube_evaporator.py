from collections import namedtuple
from scipy import optimize
from hvac import Quantity
from hvac.fluids import FluidState, HumidAir
from .plain_fin_tube_boiling_evaporator import PlainFinTubeCounterFlowBoilingEvaporator
from .plain_fin_tube_superheat_evaporator import PlainFinTubeCounterFlowSuperheatEvaporator


Q_ = Quantity

Result = namedtuple(
    'Result',
    [
        'm_dot_rfg',
        'rfg_out',
        'air_out',
        'Q_dot',
        'eps',
        'dP_air',
        'L2_superheat',
        'L2_boiling'
    ]
)


class PlainFinTubeCounterFlowEvaporator:
    """Model of single-phase, plain fin-tube, counterflow air evaporator.
    Aim of the model is to find the mass flow rate of refrigerant such that
    the refrigerant at the evaporator outlet has the degree of superheat set
    by the expansion device under the given operating conditions.
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
        """Creates a `PlainFinTubeCounterFlowEvaporator` instance by
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
        self.L2 = N_r * S_l
        self._geo_dict = {
            'L1': L1, 'L3': L3,
            'S_t': S_t, 'S_l': S_l,
            'D_i': D_i, 'D_o': D_o,
            't_f': t_f, 'N_f': N_f,
            'k_f': k_f
        }
        self.superheating_part = PlainFinTubeCounterFlowSuperheatEvaporator(**self._geo_dict)
        self.boiling_part = PlainFinTubeCounterFlowBoilingEvaporator(**self._geo_dict)

        self.air_in: HumidAir | None = None
        self.m_dot_air: Quantity | None = None
        self.rfg_in: FluidState | None = None
        self.dT_rfg_sh: Quantity | None = None
        self.rfg_out: FluidState | None = None
        self.m_dot_rfg_ini: Quantity | None = None
        self.L2_superheat: Quantity | None = None
        self.L2_boiling: Quantity | None = None
        self.air_out: HumidAir | None = None
        self.Q_dot: Quantity | None = None
        self.eps: Quantity | None = None
        self.m_dot_rfg: Quantity | None = None
        self.dP_air: Quantity | None = None

    def set_operating_conditions(
        self,
        air_in: HumidAir,
        m_dot_air: Quantity,
        rfg_in: FluidState,
        dT_rfg_sh: Quantity,
        m_dot_rfg_ini: Quantity
    ) -> None:
        """Sets the operating conditions of the evaporator.

        Parameters
        ----------
        air_in:
            State of humid air at the inlet of the evaporator.
        m_dot_air:
            Mass flow rate of air through the evaporator.
        rfg_in:
            State of refrigerant at the inlet of the evaporator; a saturated
            liquid-vapor mixture with a given pressure or temperature
            (= evaporator pressure/temperature) and a given vapor quality.
        dT_rfg_sh:
            Degree of superheat of the refrigerant at the evaporator outlet;
            the degree of superheat is determined by the setting on the expansion
            device.
        m_dot_rfg_ini:
            Initial guess of the mass flow rate of refrigerant (maybe
            determined from a preliminary calculation).
        """
        self.rfg_out = None
        self.L2_superheat = None
        self.L2_boiling = None
        self.air_out = None
        self.Q_dot = None
        self.eps = None
        self.m_dot_rfg = None
        self.dP_air = None
        self.air_in = air_in
        self.m_dot_air = m_dot_air.to('kg / s')
        self.rfg_in = rfg_in
        self.dT_rfg_sh = dT_rfg_sh.to('K')
        self.m_dot_rfg_ini = m_dot_rfg_ini.to('kg / s')
        # Set fixed operating conditions on superheating part:
        self.superheating_part.set_fixed_operating_conditions(
            air_in=self.air_in,
            m_dot_air=self.m_dot_air,
            Rfg=self.rfg_in.fluid,
            P_rfg_sat=self.rfg_in.P,
            dT_rfg_sh=self.dT_rfg_sh
        )
        # Set fixed operating conditions on boiling part:
        self.boiling_part.set_fixed_operating_conditions(
            m_dot_air=self.m_dot_air,
            rfg_in=self.rfg_in
        )

    def rate(
        self,
        i_max: int = 100,
        tol_m_dot_rfg: Quantity = Q_(0.001, 'kg / s')
    ) -> Result:
        """Finds by iteration the actual mass flow rate of refrigerant under
        given operating conditions, so that the refrigerant at the outlet of the
        evaporator has the degree of superheat set on the expansion device.
        Once the actual mass flow rate of refrigerant is found, the state of
        air the evaporator, the heat transfer rate, the air-side
        pressure drop and the flow lengths of the boiling and superheated region
        will also be determined.

        Parameters
        ----------
        i_max:
            Maximum number of iterations.
        tol_m_dot_rfg:
            Acceptable deviation between the last and previous calculated value
            of the refrigerant mass flow rate.

        Raises
        ------
        ValueError if no acceptable solution for the refrigerant mass flow rate
        could be found after `i_max` iterations.

        Returns
        -------
        namedtuple `Result`:
            m_dot_rfg: Quantity
                Mass flow rate of refrigerant.
            rfg_out: FluidState
                State of refrigerant leaving the evaporator.
            air_out: HumidAir
                State of air leaving the evaporator.
            Q: Quantity
                Heat transfer rate from air flow to refrigerant.
            eps: Quantity
                Heat transfer effectiveness of evaporator.
            dP_air: Quantity
                Air-side pressure drop across evaporator.
            L2_superheat: Quantity
                Flow length of superheated region.
            L2_boiling: Quantity
                Flow length of boiling region.
        """
        def _eq(m_dot_rfg: float) -> float:
            m_dot_rfg = Q_(m_dot_rfg, 'kg / hr')
            # Set mass flow rate through superheating part:
            self.superheating_part.set_mass_flow_rate_refrigerant(m_dot_rfg)
            # Determine the flow length required to superheat the refrigerant
            # vapor to the degree set on the expansion device:
            self.L2_superheat = self.superheating_part.determine_flow_length(self.L2)
            # Determine available flow length for boiling the refrigerant:
            self.L2_boiling = self.L2 - self.L2_superheat
            self.boiling_part.set_flow_length(self.L2_boiling)
            # Set state of air entering the boiling part = state of air leaving
            # the superheating part:
            self.boiling_part.set_air_in(self.superheating_part.air_out)
            # Determine the mass flow rate of refrigerant required to transform
            # the refrigerant entering the evaporator into saturated vapor:
            m_dot_rfg_new, *_ = self.boiling_part.rate(m_dot_rfg, i_max, tol_m_dot_rfg)
            # Determine deviation between new and previous mass flow rate:
            dev_m_dot_rfg = m_dot_rfg_new - m_dot_rfg
            return dev_m_dot_rfg.to('kg / hr').m

        sol = optimize.root_scalar(
            _eq,
            method='secant',
            x0=self.m_dot_rfg_ini.to('kg / hr').m,
            maxiter=i_max,
            xtol=tol_m_dot_rfg.to('kg / hr').m
        )
        self.m_dot_rfg = Q_(sol.root, 'kg / hr')
        self.rfg_out = self.superheating_part.rfg_out
        self.air_out = self.boiling_part.air_out
        self.Q_dot = self.boiling_part.Q_dot + self.superheating_part.Q_dot
        self.eps = self._get_eps()
        self.dP_air = self.boiling_part.dP_air + self.superheating_part.dP_air
        return Result(
            self.m_dot_rfg,
            self.rfg_out,
            self.air_out,
            self.Q_dot,
            self.eps,
            self.dP_air,
            self.L2_superheat,
            self.L2_boiling
        )

    def _get_eps(self) -> Quantity:
        """Returns heat transfer effectiveness of the evaporator."""
        air_sat_out = HumidAir(Tdb=self.rfg_in.T, RH=Q_(100, 'pct'))
        Q_max = self.m_dot_air * (self.air_in.h - air_sat_out.h)
        eps = self.Q_dot / Q_max
        return eps.to('frac')

    def __call__(
        self,
        air_in: HumidAir,
        m_dot_air: Quantity,
        rfg_in: FluidState,
        dT_rfg_sh: Quantity,
        m_dot_rfg_ini: Quantity,
        i_max: int = 100,
        tol_m_dot_rfg: Quantity = Q_(0.001, 'kg / s')
    ) -> Result:
        """Combines the methods `set_operating_conditions` and `rate` in a
        single method.
        """
        self.set_operating_conditions(
            air_in, m_dot_air, rfg_in,
            dT_rfg_sh, m_dot_rfg_ini
        )
        res = self.rate(i_max, tol_m_dot_rfg)
        return res
