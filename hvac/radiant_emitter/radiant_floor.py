from abc import ABC
import math
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import Fluid


Q_ = Quantity
Water = Fluid('Water')


class Floor(ABC):
    """Abstract base class from which specific types of floors are derived that
    can be used in the class `RadiantFloorPanel` to:
    -   get the thermal resistance of the floor covering
    -   get the heat transfer coefficient between the average water temperature
        in the floor heating circuit and the indoor air temperature.
    -   get the linear heat loss coefficient along the edges of a heated slab on
        grade
    -   get the thermal resistance of the insulation underneath the radiant
        floor panel.

    References
    ----------
    [1] Siegenthaler, J. (2022). MODERN HYDRONIC HEATING AND COOLING: FOR
    RESIDENTIAL AND LIGHT COMMERCIAL BUILDINGS. Cengage Learning.
    """
    R_fl_cover: Quantity = None
    # thermal resistance of floor covering
    heat_transfer_coeff_dict = {
        10: float('nan'),
        15: float('nan'),
        20: float('nan'),
        30: float('nan')
    }
    # the dict keys are tube spacing values in cm, the dict values will be the
    # corresponding heat transfer coefficients between average heating water
    # temperature and indoor air temperature in W/(m².K) that can be derived
    # from heat output design tables published by underfloor heating system
    # manufacturers.
    edge_loss_coeff_dict = {
        0: 0.83,
        5: 0.39,
        7.5: 0.30,
        10: 0.23
    }
    # used to calculate the downward and edge loss coefficient from a heated
    # slab-on-grade floor; the keys are R-values of edge insulation in units
    # (°F.hr.ft²)/Btu, the values are corresponding linear heat loss coefficients
    # in units Btu/(hr.ft.°F). Source [1]

    def __init__(
        self,
        S: Quantity,
        R_edge_ins: Quantity | None = None,
        R_fl_ins: Quantity | None = None
    ) -> None:
        """Creates a `Floor` instance.

        Parameters
        ----------
        S:
            Distance between floor heating tubes. (Suggested tube spacings can
            be found in Ref. [1], Chapter 10, p. 542, Figure 10-119.)
        R_edge_ins: optional
            Thermal resistance of floor edge insulation.
        R_fl_ins: optional
            Thermal resistance of insulation installed beneath the radiant
            floor panel.
        """
        self.S = S
        self.R_edge_ins = R_edge_ins
        self.R_fl_ins = R_fl_ins

        self._calculate_upward_heat_transfer_coeff()
        self._calculate_edge_loss_coeff()

    def _calculate_upward_heat_transfer_coeff(self):
        # heat transfer coefficient between average heating water temperature
        # and indoor air temperature
        S = self.S.to('cm').m
        S_lst = list(self.heat_transfer_coeff_dict.keys())
        k_lst = list(self.heat_transfer_coeff_dict.values())
        k = interp1d(S_lst, k_lst)(S)
        self.k = Q_(float(k), 'W / (m ** 2 * K)')

    def _calculate_edge_loss_coeff(self):
        self.m = None
        if self.R_edge_ins is not None:
            R_ins_edge = self.R_edge_ins.to('delta_degF * (ft ** 2) * hr / Btu').m
            R_ins_edge_lst = list(self.edge_loss_coeff_dict.keys())
            m_lst = list(self.edge_loss_coeff_dict.values())
            m = interp1d(R_ins_edge_lst, m_lst)(R_ins_edge)
            self.m = Q_(m, 'Btu / (hr * ft * delta_degF)')


class Tiles(Floor):
    R_fl_cover = Q_(0.010, 'm ** 2 * K / W')
    heat_transfer_coeff_dict = {
        10: 6.913,
        15: 6.074,
        20: 5.274,
        30: 4.000
    }
    # Note: heat transfer coefficients were derived from a floor heating
    # design table from manufacturer VASCO (www.vasco.eu).


class PVC(Floor):
    R_fl_cover = Q_(0.025, 'm ** 2 * K / W')
    heat_transfer_coeff_dict = {
        10: 6.264,
        15: 5.536,
        20: 4.874,
        30: 3.736
    }
    # Note: heat transfer coefficients were derived from a floor heating
    # design table from VASCO (www.vasco.eu).


class Parquet(Floor):
    R_fl_cover = Q_(0.050, 'm ** 2 * K / W')
    heat_transfer_coeff_dict = {
        10: 5.400,
        15: 4.864,
        20: 4.313,
        30: 3.400
    }
    # Note: heat transfer coefficients were derived from a floor heating
    # design table from VASCO (www.vasco.eu).


class ThinCarpet(Floor):
    R_fl_cover = Q_(0.075, 'm ** 2 * K / W')
    heat_transfer_coeff_dict = {
        10: 4.761,
        15: 4.313,
        20: 3.887,
        30: 3.113
    }
    # Note: heat transfer coefficients were derived from a floor heating
    # design table from VASCO (www.vasco.eu).


class ThickCarpet(Floor):
    R_fl_cover = Q_(0.150, 'm ** 2 * K / W')
    heat_transfer_coeff_dict = {
        10: 3.487,
        15: 3.239,
        20: 3.000,
        30: 2.487
    }
    # Note: heat transfer coefficients were derived from a floor heating
    # design table from VASCO (www.vasco.eu).


class RadiantFloorPanel:
    """Implements design and analysis procedures for radiant floor heating.

    References
    ----------
    [1] Siegenthaler, J. (2022). MODERN HYDRONIC HEATING AND COOLING: FOR
    RESIDENTIAL AND LIGHT COMMERCIAL BUILDINGS. Cengage Learning.
    """
    h_fl = Q_(11.36, 'W / (m ** 2 * K)')
    # combined convection and radiation heat transfer coefficient at the
    # interface between floor surface and space air. Source [1].

    def __init__(
        self,
        floor_cover: str,
        Q_dot_load_des: Quantity,
        A_fl: Quantity,
        T_int_des: Quantity,
        S: Quantity,
        T_fl_max: Quantity,
        L_cir_max: Quantity,
        L_ldr: Quantity,
        DT_w_des: Quantity = Q_(10, 'K'),
        T_ext_des: Quantity | None = None,
        T_uf: Quantity | None = None,
        P_fl: Quantity | None = None,
        R_edg_ins: Quantity | None = None,
        R_ufl_ins: Quantity | None = None
    ) -> None:
        """Creates an instance of `RadiantFloorPanel`.

        Parameters
        ----------
        floor_cover: ['tiles', 'pvc', 'parquet', 'thin_carpet', 'thick_carpet']
            Type of floor covering.
        Q_dot_load_des:
            Heat loss of the space at design conditions, i.e. design heating
            load, but excluding the heat loss from the underside of the floor.
        A_fl:
            Available floor area for underfloor heating.
        T_int_des:
            Indoor air temperature at design conditions.
        T_fl_max:
            The maximum allowable floor surface temperature.
        L_cir_max:
            The maximum allowable length of one floor heating circuit including
            leader length (supply + return); this value will depend on the
            diameter of the floor heating tubing (see [1], Chapter 10, p. 524,
            Figure 10-99).
        L_ldr:
            The leader length required for the circuit(s) to reach the manifold
            station.
        DT_w_des: default 10 K
            Temperature difference between heating water inlet and outlet
            temperature.
        T_ext_des: default None
            Outdoor air temperature at design conditions. If set, a slab-on-
            grade floor is assumed for calculating the edge and downward heat
            loss from the heated slab.
        T_uf: default None
            Temperature underneath the heated floor. If set, a suspended floor
            deck is assumed to calculate the downward heat loss from the floor.
        P_fl: default None
            Perimeter of heated floor area. Only needs to be set in case of
            a heated slab-on-grade floor (i.e. when `T_ext_des` is also set).
        R_edg_ins: default None
            Thermal resistance of floor edge insulation. Only needs to be set in
            case of a heated slab-on-grade floor (i.e. when `T_ext_des` is also set).
        R_ufl_ins: default None
            Thermal resistance of insulation installed beneath the radiant
            floor panel. Only needs to be set in case of a suspended floor deck
            (i.e. when `T_uf` is set).
        """
        self.floor = self._get_floor(floor_cover, S, R_edg_ins, R_ufl_ins)
        self.Q_dot_load_des = Q_dot_load_des
        self.A_fl = A_fl
        self.T_int_des = T_int_des
        self.T_fl_max = T_fl_max
        self.L_cir_max = L_cir_max
        self.L_ldr = L_ldr
        self.DT_w_des = DT_w_des.to('K')
        self.T_ext_des = T_ext_des
        self.T_uf = T_uf
        self.P_fl = P_fl

        # Design procedure:
        self._calculate_floor_surface_temperature()
        self._calculate_additional_heat_output()
        self._calculate_water_supply_temperature()
        self.Q_dot_fl_loss = Q_(0.0, 'W')  # downward/edge heat loss
        if self.T_ext_des is not None:
            self._calculate_slab_on_grade_heat_loss()
        elif self.T_uf is not None:
            self._calculate_suspended_floor_heat_loss()
        self._calculate_water_volume_flow_rate()
        self._calculate_number_of_circuits()

    @staticmethod
    def _get_floor(
        floor_cover: str,
        S: Quantity,
        R_edg_ins: Quantity,
        R_ufl_ins: Quantity
    ) -> Floor:
        # Creates and returns a concrete instance of type `Floor` specified by
        # parameter `floor_cover`.
        match floor_cover:
            case 'tiles':
                return Tiles(S, R_edg_ins, R_ufl_ins)
            case 'pvc':
                return PVC(S, R_edg_ins, R_ufl_ins)
            case 'parquet':
                return Parquet(S, R_edg_ins, R_ufl_ins)
            case 'thin_carpet':
                return ThinCarpet(S, R_edg_ins, R_ufl_ins)
            case 'thick_carpet':
                return ThickCarpet(S, R_edg_ins, R_ufl_ins)

    def _calculate_floor_surface_temperature(self):
        # Calculates the required heat flux emitted by the floor surface to
        # compensate for the heat loss of the space at design conditions:
        self.q_dot_des = self.Q_dot_load_des / self.A_fl
        # Calculate the required floor surface temperature to establish the
        # required heat flux q_dot_des:
        self.T_fl = self.T_int_des.to('K') + self.q_dot_des / self.h_fl

    def _calculate_additional_heat_output(self):
        # Calculate the heat flux and total heat output emitted by the floor
        # surface at the maximum allowable floor surface temperature:
        self.q_dot_max = self.h_fl * (self.T_fl_max - self.T_int_des)
        self.Q_dot_max = self.q_dot_max * self.A_fl
        # If the maximum allowable heat output is smaller than the design load,
        # additional heating will be needed at design conditions:
        self.Q_dot_add = max(self.Q_dot_load_des - self.Q_dot_max, Q_(0, 'W'))

    def _calculate_water_supply_temperature(self):
        # Calculates the required water supply temperature given:
        # - design heating load flux
        # - indoor air temperature at design conditions
        # - the selected water temperature drop along the floor heating circuit
        # - the heat transfer coefficient between the average water temperature
        #   and the indoor air temperature
        self.T_w_sup = (
            (self.T_int_des.to('K') + self.DT_w_des / 2)
            + self.q_dot_des / self.floor.k
        )

    def _calculate_slab_on_grade_heat_loss(self):
        # Calculates the downward and edge heat loss from a heated slab on grade.
        if self.floor.m is not None:
            m = self.floor.m.to('W / (m * K)')
            P_fl = self.P_fl.to('m')
            q_dot_des = self.q_dot_des.to('W / m ** 2')
            R_cov = self.floor.R_fl_cover.to('m ** 2 * K / W')
            R_con = (1 / self.h_fl).to('m ** 2 * K / W')
            R_fl = R_cov + R_con
            T_int_des = self.T_int_des.to('K')
            T_ext_des = self.T_ext_des.to('K')
            self.Q_dot_fl_loss = m * P_fl * (q_dot_des * R_fl + T_int_des - T_ext_des)

    def _calculate_suspended_floor_heat_loss(self):
        # Calculates the downward heat loss from a suspended floor deck.
        if self.floor.R_fl_ins is not None:
            A_fl = self.A_fl.to('m ** 2')
            R_fl_ins = self.floor.R_fl_ins.to('m ** 2 * K / W')
            q_dot_des = self.q_dot_des.to('W / m ** 2')
            R_cov = self.floor.R_fl_cover.to('m ** 2 * K / W')
            R_con = (1 / self.h_fl).to('m ** 2 * K / W')
            R_fl = R_cov + R_con
            T_id = self.T_int_des.to('K')
            T_uf = self.T_uf.to('K')
            self.Q_dot_fl_loss = (A_fl / R_fl_ins) * (q_dot_des * R_fl + T_id - T_uf)

    def _calculate_water_volume_flow_rate(self):
        # Calculates the required flow rate through the circuit.
        Q_dot_tot = self.Q_dot_load_des + self.Q_dot_fl_loss
        # required upward heat + downward/edge heat loss from floor slab
        T_w_avg = self.T_w_sup - self.DT_w_des / 2
        w = Water(T=T_w_avg, P=Q_(1.5, 'bar'))
        self.V_dot_w = Q_dot_tot / (w.rho * w.cp * self.DT_w_des)

    def _calculate_number_of_circuits(self):
        L_cir = self.A_fl / self.floor.S
        N_cir = math.ceil((L_cir + self.L_ldr) / self.L_cir_max)
        while N_cir > 1:
            L_cir /= N_cir
            if L_cir + self.L_ldr <= self.L_cir_max:
                break
            N_cir += 1
        self.N_cir = N_cir
        self.L_cir = L_cir
        self.V_dot_w /= N_cir

    def design(self) -> dict[str, Quantity | int]:
        """Returns a dict with the results of the design calculations:
        'T_w_sup':
            The required supply water temperature.
        'V_dot_w':
            The required water volume flow rate.
        'T_fl':
            The resulting floor surface temperature.
        'N_cir':
            The required number of floor circuits.
        'L_cir':
            The circuit tube length (in the space only).
        'Q_dot_add':
            Additional heating to support floor heating under design conditions.
        """
        return {
            'T_w_sup': self.T_w_sup,
            'V_dot_w': self.V_dot_w,
            'T_fl': self.T_fl,
            'N_cir': self.N_cir,
            'L_cir': self.L_cir,
            'Q_dot_add': self.Q_dot_add
        }

    def Q_dot(
        self,
        T_w_sup: Quantity,
        V_dot_w: Quantity,
        T_i: Quantity | None = None,
        T_fl: Quantity | None = None,
        f_Ql: Quantity = Q_(10, 'pct')
    ) -> Quantity:
        """Returns the steady-state heat output of the underfloor heating
        circuit(s) for the given combination of input conditions.

        Parameters
        ----------
        T_w_sup:
            Water supply temperature at the inlet of the circuit.
        V_dot_w:
            Volume flow rate of water through the circuit.
        T_i: default None
            Indoor air temperature of the space.
        T_fl: default None
            Average floor surface temperature.
        f_Ql: default 10 %
            Floor heating loss coefficient: the heating water looses more
            heat than heat is transferred to space. By default, it is assumed
            that 10% of the heat given off by the water does not enter the space.

        Notes
        -----
        Either `T_i`, or `T_fl` cannot be None. If both are specified, `T_fl` is
        ignored. If both are None, the method also returns None.
        """
        if T_i is not None:
            try:
                return self._Q1(T_w_sup, V_dot_w, T_i, f_Ql)
            except ValueError:
                return Q_(0.0, 'W')
        elif T_fl is not None:
            return self._Q2(T_w_sup, V_dot_w, T_fl, f_Ql)

    def _Q1(
        self,
        T_w_sup: Quantity,
        V_dot_w: Quantity,
        T_i: Quantity,
        f_Ql: Quantity = Q_(10, 'pct')
    ) -> Quantity:
        """Returns the heat output in case `T_i` is given."""
        S = self.floor.S.to('m').m
        L_cir = self.L_cir.to('m').m * self.N_cir
        k = self.floor.k.to('W / (m ** 2 * K)').m
        f_Ql = 1.0 + f_Ql.to('frac').m
        K = f_Ql * S * L_cir * k
        T_i = T_i.to('K').m
        T_w_sup = T_w_sup.to('K').m
        V_dot_w = V_dot_w.to('m ** 3 / s').m
        w = Water(T=Q_(T_w_sup, 'K'), P=Q_(1.5, 'bar'))
        rho_w = w.rho.to('kg / m ** 3').m
        cp_w = w.cp.to('J / (kg * K)').m
        Q_dot_max = rho_w * cp_w * V_dot_w * (T_w_sup - T_i)

        def eq(Q_dot: float) -> float:
            DT_max = T_w_sup - T_i
            DT_w = Q_dot / (rho_w * cp_w * V_dot_w)
            T_w_ret = T_w_sup - DT_w
            DT_min = T_w_ret - T_i
            Q_dot_new = K * (DT_max - DT_min) / np.log(DT_max / DT_min)
            return Q_dot - Q_dot_new

        sol = root_scalar(eq, bracket=[0.01, Q_dot_max-0.01])
        Q = Q_(sol.root, 'W')
        return Q

    def _Q2(
        self,
        T_w_sup: Quantity,
        V_w_dot: Quantity,
        T_fl: Quantity,
        f_Ql: Quantity = Q_(10, 'pct')
    ) -> Quantity:
        """Returns the heat output in case `T_fl` is given."""
        w = Water(T=T_w_sup, P=Q_(1.5, 'bar'))
        rw = w.rho.to('kg / m ** 3')
        cw = w.cp.to('J / (kg * K)')
        f_Ql = Q_(1.0, 'frac') + f_Ql.to('frac')
        R_wf_t = (1 / self.floor.k - 1 / self.h_fl) / self.A_fl
        # total thermal resistance in units K/W between heating water and floor
        # surface, used in the equation Q = (Tw_avg - T_fl) / R_wf_t
        n = 2 * f_Ql * rw * cw * V_w_dot
        d = f_Ql + 2 * rw * cw * V_w_dot * R_wf_t
        K = n / d
        Q_dot = K * (T_w_sup - T_fl)
        return Q_dot
