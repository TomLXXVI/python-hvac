from collections import namedtuple
from hvac import Quantity
from hvac.fluids import FluidState, HumidAir
from .plain_fin_tube_boiling_evaporator import PlainFinTubeCounterFlowBoilingEvaporator
from .plain_fin_tube_superheat_evaporator import PlainFinTubeCounterFlowSuperheatEvaporator

Q_ = Quantity

Result = namedtuple(
    'Result',
    ['m_dot_rfg', 'rfg_out', 'air_out', 'Q', 'eps', 'dP_air', 'L2_superheat']
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
        kwargs = {
            'L1': L1, 'L3': L3,
            'S_t': S_t, 'S_l': S_l,
            'D_i': D_i, 'D_o': D_o,
            't_f': t_f, 'N_f': N_f,
            'k_f': k_f
        }
        self.evp_superheat = PlainFinTubeCounterFlowSuperheatEvaporator(**kwargs)
        self.evp_boiling = PlainFinTubeCounterFlowBoilingEvaporator(**kwargs)

        self.air_in: HumidAir | None = None
        self.m_dot_air: Quantity | None = None
        self.rfg_in: FluidState | None = None
        self.dT_rfg_sh: Quantity | None = None
        self.rfg_out: FluidState | None = None
        self.m_dot_rfg_ini: Quantity | None = None

        self.L2_superheat: Quantity | None = None
        self.air_out: HumidAir | None = None
        self.Q: Quantity | None = None
        self.eps: float | None = None
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
            Initial guess of the mass flow rate of refrigerant determined from
            a preliminary calculation.
        """
        self.air_in = air_in
        self.m_dot_air = m_dot_air.to('kg / s')
        self.rfg_in = rfg_in
        self.dT_rfg_sh = dT_rfg_sh.to('K')
        self.m_dot_rfg_ini = m_dot_rfg_ini.to('kg / s')
        self.evp_superheat.set_fixed_operating_conditions(
            air_in=self.air_in,
            m_dot_air=self.m_dot_air,
            Rfg=self.rfg_in.fluid,
            T_rfg_sat=self.rfg_in.T,
            dT_rfg_sh=self.dT_rfg_sh
        )
        self.evp_boiling.set_fixed_operating_conditions(
            m_dot_air=self.m_dot_air,
            rfg_in=self.rfg_in
        )

    def rate(
        self,
        i_max: int = 10,
        tol_m_dot_rfg: Quantity = Q_(0.01, 'kg / s')
    ) -> Result:
        """Finds by iteration the actual mass flow rate of refrigerant under the
        given operating conditions, such that the refrigerant at the outlet of the
        evaporator has the degree of superheat dictated by the setting on
        the expansion device. Once the actual mass flow rate of refrigerant has
        been determined, the state of outlet air, heat transfer rate, air-side
        pressure drop and flow length of superheated region are also determined.

        Parameters
        ----------
        i_max:
            Maximum number of iterations.
        tol_m_dot_rfg:
            Tolerated deviation between the last calculated value of refrigerant
            mass flow rate and the previous calculated value. If the last
            calculated value is within the tolerance band of the previous value,
            iteration is stopped.

        Raises
        ------
        ValueError if no acceptable solution for the refrigerant mass flow rate
        could be found after `i_max` iterations.

        Returns
        -------
        namedtuple `Result` with members:
            m_dot_rfg:
                Mass flow rate of refrigerant.
            air_out:
                State of air at evaporator outlet.
            Q:
                Heat transfer rate from air flow to refrigerant.
            dP_air:
                Air-side pressure drop
            L2_superheat:
                Flow length of superheated region.
        """
        m_dot_rfg = self.m_dot_rfg_ini
        i = 0
        while i < i_max:
            self.evp_superheat.set_mass_flow_rate_refrigerant(m_dot_rfg)
            L2_superheat = self.evp_superheat.determine_flow_length(self.L2)
            L2_boiling = self.L2 - L2_superheat
            self.evp_boiling.set_air_in(self.evp_superheat.air_out)
            self.evp_boiling.set_flow_length(L2_boiling)
            m_dot_rfg_new, *_ = self.evp_boiling.rate(m_dot_rfg, i_max, tol_m_dot_rfg)
            dev_m_dot_rfg = abs(m_dot_rfg_new - m_dot_rfg)
            if dev_m_dot_rfg <= tol_m_dot_rfg:
                self.L2_superheat = L2_superheat
                self.air_out = self.evp_boiling.air_out
                self.rfg_out = self.evp_superheat.rfg_out
                self.Q = self.evp_boiling.Q + self.evp_superheat.Q
                self.eps = self._get_eps()
                self.m_dot_rfg = m_dot_rfg_new
                self.dP_air = self.evp_boiling.dP_air + self.evp_superheat.dP_air
                return Result(
                    self.m_dot_rfg, self.rfg_out, self.air_out, self.Q,
                    self.eps, self.dP_air, self.L2_superheat
                )
            m_dot_rfg = m_dot_rfg_new
            i += 1
        else:
            raise ValueError(
                "no acceptable solution found within tolerance "
                "{tol_m_dot_rfg.to('kg / s'):~P} after {i_max} "
                "iterations"
            )

    def _get_eps(self) -> float:
        air_sat_out = HumidAir(Tdb=self.rfg_in.T, RH=Q_(100, 'pct'))
        Q_max = self.m_dot_air * (self.air_in.h - air_sat_out.h)
        eps = self.Q / Q_max
        return eps

    def __call__(
        self,
        air_in: HumidAir,
        m_dot_air: Quantity,
        rfg_in: FluidState,
        dT_rfg_sh: Quantity,
        m_dot_rfg_ini: Quantity,
        i_max: int = 10,
        tol_m_dot_rfg: Quantity = Q_(0.01, 'kg / s')
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
