from typing import Any
import warnings
import numpy as np
from scipy import optimize
from hvac import Quantity
from hvac.fluids import Fluid, HumidAir, CoolPropWarning, FluidState
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube
from hvac.heat_transfer.finned_surface.fins import Fin, PlainContinuousFin
import hvac.heat_exchanger.recuperator.general.eps_ntu_wet as wet
import hvac.heat_exchanger.recuperator.general.eps_ntu as dry
from hvac.heat_exchanger.recuperator.general import corrections
from .geometry import ContinuousFinStaggeredTubeBank
from .correlations import PlainContinuousFinStaggeredTubeBank as correlations


warnings.filterwarnings('ignore', category=CoolPropWarning)


Q_ = Quantity
Water = Fluid('Water')


class HeatExchangerCore:
    # Continuous plain fin, staggered tube bank.

    def __init__(
        self,
        width: Quantity,
        height: Quantity,
        num_rows: int,
        pitch_trv: Quantity,
        pitch_lon: Quantity,
        d_o: Quantity,
        d_i: Quantity,
        t_fin: Quantity,
        fin_density: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
        d_r: Quantity | None = None,
        num_circuits: int | None = None
    ) -> None:
        self.geometry = ContinuousFinStaggeredTubeBank(
            width, height, num_rows, pitch_trv,
            pitch_lon, d_o, d_i, t_fin, fin_density,
            k_fin, d_r
        )
        self.internal = InternalSurface(self)
        self.external = ExternalSurface(self)
        self.num_circuits = num_circuits

        self.T_wall: Quantity | None = None
        self.h_int: Quantity | None = None
        self.R_int: Quantity | None = None
        self.h_ext: Quantity | None = None
        self.R_ext: Quantity | None = None
        self.R_tot: Quantity | None = None
        self.UA: Quantity | None = None
        self.eta_surf: float | None = None
        self.dP_int: Quantity | None = None
        self.dP_ext: Quantity | None = None

    def hex_properties(
        self,
        air_mean: HumidAir,
        air_in: HumidAir,
        air_out: HumidAir,
        water_mean: FluidState,
        air_m_dot: Quantity,
        water_m_dot: Quantity,
    ) -> dict[str, Quantity | float]:
        self.internal.m_dot = water_m_dot
        self.internal.water = water_mean
        self.external.m_dot = air_m_dot
        self.external.air = air_mean
        self.T_wall = self.__find_wall_temperature()
        self.h_int = self.internal.heat_transfer_coeff()
        self.R_int = self.internal.thermal_resistance(self.h_int)
        self.h_ext = self.external.heat_transfer_coeff(self.T_wall)
        self.R_ext = self.external.thermal_resistance(self.h_ext)
        self.R_tot = self.R_int + self.R_ext
        # Note: The thermal conduction resistance of the heat exchanger body
        # is ignored, and fouling is not being taken into account.
        self.UA = 1 / self.R_tot
        self.dP_int = self.internal.pressure_drop(self.T_wall)
        self.dP_ext = self.external.pressure_drop(air_in.rho, air_out.rho, self.T_wall)
        self.eta_surf = self.external.eta_surf(self.h_ext)
        return {
            'h_int': self.h_int,
            'h_ext': self.h_ext,
            'R_int': self.R_int,
            'R_ext': self.R_ext,
            'R_tot': self.R_tot,
            'UA': self.UA,
            'eta_surf': self.eta_surf,
            'dP_int': self.dP_int,
            'dP_ext': self.dP_ext
        }

    def __find_wall_temperature(self) -> Quantity:

        def __eq__(T_wall_: float) -> float:
            T_wall_ = Q_(T_wall_, 'K')
            h_int = self.internal.heat_transfer_coeff()
            R_int = self.internal.thermal_resistance(h_int)
            h_ext = self.external.heat_transfer_coeff(T_wall_)
            R_ext = self.external.thermal_resistance(h_ext)
            T_int_ = self.internal.water.T.to('K')
            T_ext_ = self.external.air.Tdb.to('K')
            n = T_int_ / R_int + T_ext_ / R_ext
            d = 1 / R_int + 1 / R_ext
            T_wall_new = n / d
            dev = (T_wall_new - T_wall_).to('K').m
            return dev

        T_int = self.internal.water.T.to('K').m  # cold fluid (water)
        T_ext = self.external.air.Tdb.to('K').m  # hot fluid (air)
        sol = optimize.root_scalar(
            __eq__,
            bracket=(T_int, T_ext)
        )
        T_wall = Q_(sol.root, 'K')
        return T_wall


class InternalSurface:

    def __init__(self, parent: HeatExchangerCore):
        self.parent = parent
        self.water: FluidState | None = None
        self.m_dot: Quantity | None = None

    def __m_dot_circuit(self) -> Quantity:
        A_min = self.parent.geometry.internal.A_min
        n_r = self.parent.geometry.num_rows
        d_i = self.parent.geometry.d_i
        n_c = self.parent.num_circuits
        n_t = self.parent.geometry.num_tubes_1st_row
        # Calculate the mass flow rate through a single tube (this assumes that
        # the number of parallel circuits equals the number of tubes in the
        # first row):
        A_min_row = A_min / n_r  # min. free flow area of 1 row
        G = self.m_dot / A_min_row
        A_tube = np.pi * d_i ** 2 / 4
        m_dot_circuit = G * A_tube
        # If the number of parallel water circuits has been specified, the mass
        # flow rate through a single tube is multiplied with the number of
        # "tube heights" in one circuit to get the mass flow rate per circuit:
        if isinstance(n_c, int) and n_c >= 1:
            m_dot_circuit *= (n_t / n_c)
        return m_dot_circuit

    def heat_transfer_coeff(self) -> Quantity:
        tube = CircularTube(
            Di=self.parent.geometry.d_i,
            L=self.parent.geometry.length,
            fluid=self.water
        )
        tube.m_dot = self.__m_dot_circuit()
        h = tube.avg_heat_transfer_coefficient()
        if isinstance(h, tuple):
            h = h[0]  # if laminar flow, assume constant heat flux
        return h.to('W / (m ** 2 * K)')

    def thermal_resistance(self, h: Quantity) -> Quantity:
        A_tot = self.parent.geometry.internal.A_tot
        R = 1 / (h * A_tot)
        return R

    def pressure_drop(self, T_wall: Quantity) -> Quantity:
        tube = CircularTube(
            Di=self.parent.geometry.d_i,
            L=self.parent.geometry.length,
            fluid=self.water
        )
        tube.m_dot = self.__m_dot_circuit()
        ff = tube.friction_factor() / 4
        flow_regime = tube.get_flow_condition()
        ff = corrections.correct_friction_factor(
            ff, T_wall, flow_regime,
            'heating', self.water
        )
        w = self.parent.geometry.width
        n_r = self.parent.geometry.num_rows
        n_c = self.parent.num_circuits  # number of parallel water circuits
        L = w * n_r  # tube length of one 'tube height'
        if isinstance(n_c, int) and n_c >= 1:
            n_t = self.parent.geometry.num_tubes_1st_row
            L *= (n_t / n_c)
        d_h = tube.Dh
        u_m = tube.mean_velocity()
        rho = self.water.rho
        dP = ff * (L / d_h) * rho * u_m ** 2 / 2
        return dP


class ExternalSurface:

    def __init__(self, parent: HeatExchangerCore):
        self.parent = parent
        self.fin = self.__configure_fin()
        self.air: HumidAir | None = None
        self.m_dot: Quantity | None = None

    def __configure_fin(self) -> Fin:
        M = L = self.parent.geometry.pitch_trv / 2
        fin = PlainContinuousFin(
            r_i=self.parent.geometry.d_r / 2,
            t=self.parent.geometry.t_fin,
            k=self.parent.geometry.k_fin,
            M=M,
            L=L
        )
        return fin

    def __mass_velocity(self) -> Quantity:
        A_min = self.parent.geometry.external.A_min
        G = self.m_dot.to('kg / s') / A_min.to('m ** 2')
        return G

    def __reynolds_number(self, d: Quantity) -> float:
        G = self.__mass_velocity()
        mu = self.air.mu.to('Pa * s')
        Re = (G * d.to('m') / mu).m
        return Re

    def __prandtl_number(self) -> float:
        rho = self.air.rho.to('kg / m ** 3').m
        mu = self.air.mu.to('Pa * s').m
        k = self.air.k.to('W / (m * K)').m
        cp = self.air.cp.to('J / (kg * K)').m
        Pr = (mu / rho) / (k / (rho * cp))
        return Pr

    def __flow_regime(self) -> str:
        Re_Dh = self.__reynolds_number(self.parent.geometry.external.d_h)
        if Re_Dh > 2300:
            return 'turbulent'
        else:
            return 'laminar'

    def __correct_heat_trf_coeff(self, Nu: float, T_wall: Quantity):
        Nu = corrections.correct_nusselt_number(
            Nu_cp=Nu,
            T_w=T_wall.to('K'),
            flow_regime=self.__flow_regime(),
            thermal_regime='cooling',
            fluid=self.air
        )
        h = Nu * self.air.k / self.parent.geometry.external.d_h
        return h.to('W / (m ** 2 * K)')

    def __fanning_friction_factor(self, T_wall: Quantity) -> float:
        ff = correlations.fanning_friction_factor(
            Re=self.__reynolds_number(self.parent.geometry.d_o),
            d_r=self.parent.geometry.d_r.to('m').m,
            p_trv=self.parent.geometry.pitch_trv.to('m').m,
            p_lon=self.parent.geometry.pitch_lon.to('m').m,
            n_fin=self.parent.geometry.fin_density.to('1 / m').m,
            n_rows=self.parent.geometry.num_rows
        )
        ff = corrections.correct_friction_factor(
            f_cp=ff,
            T_w=T_wall,
            flow_regime=self.__flow_regime(),
            thermal_regime='cooling',
            fluid=self.air
        )
        return ff

    def eta_surf(self, h: Quantity) -> float:
        self.fin.h_avg = h
        A_fin = self.parent.geometry.external.A_fin.to('m ** 2').m
        A_tot = self.parent.geometry.external.A_tot.to('m ** 2').m
        eta = 1 - (A_fin / A_tot) * (1 - self.fin.efficiency)
        return eta

    def heat_transfer_coeff(self, T_wall: Quantity) -> Quantity:
        Re_Do = self.__reynolds_number(self.parent.geometry.d_o)
        Pr = self.__prandtl_number()
        j = correlations.colburn_j_factor(
            Re=Re_Do,
            d_r=self.parent.geometry.d_r.to('m').m,
            d_h=self.parent.geometry.external.d_h.to('m').m,
            p_trv=self.parent.geometry.pitch_trv.to('m').m,
            p_lon=self.parent.geometry.pitch_lon.to('m').m,
            n_fin=self.parent.geometry.fin_density.to('1 / m').m,
            n_rows=self.parent.geometry.num_rows
        )
        Re_Dh = self.__reynolds_number(self.parent.geometry.external.d_h)
        Nu = j * Re_Dh * Pr ** (1 / 3)
        h = self.__correct_heat_trf_coeff(Nu, T_wall)
        return h.to('W / (m ** 2 * K)')

    def thermal_resistance(self, h: Quantity) -> Quantity:
        A_tot = self.parent.geometry.external.A_tot
        eta_surf = self.eta_surf(h)
        R = 1 / (eta_surf * h * A_tot)
        return R

    def pressure_drop(
        self,
        rho_in: Quantity,
        rho_out: Quantity,
        T_wall: Quantity
    ) -> Quantity:
        G = self.__mass_velocity()
        d_h = self.parent.geometry.external.d_h.to('m')
        L = self.parent.geometry.length.to('m')
        sigma = self.parent.geometry.external.sigma
        rho_avg = 2 * rho_in * rho_out / (rho_in + rho_out)
        ff = self.__fanning_friction_factor(T_wall)
        k = (G ** 2) / (2 * rho_in)
        a = ff * (4 * L / d_h) * (rho_in / rho_avg)
        b = (1 + sigma ** 2) * (rho_in / rho_out - 1)
        dP = k * (a + b)
        return dP.to('Pa')


class PlainFinTubeAirToWaterCounterFlowHeatExchanger:
    """
    Models a counterflow fin-tube heat exchanger with continuous, plain fins for
    cooling humid air by water. The model assumes the air-side heat transfer
    surface is fully wet.
    """
    def __init__(
        self,
        width: Quantity,
        height: Quantity,
        num_rows: int,
        pitch_trv: Quantity,
        pitch_lon: Quantity,
        d_o: Quantity,
        d_i: Quantity,
        t_fin: Quantity,
        fin_density: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
        d_r: Quantity | None = None,
        num_circuits: int | None = None
    ) -> None:
        """
        Creates a `PlainFinTubeAirToWaterCounterFlowHeatExchanger` object.

        Parameters
        ----------
        width:
            Width of heat exchanger core.
        height:
            Height of heat exchanger core.
        num_rows:
            The number of rows in the heat exchanger core.
        pitch_trv:
            Transverse pitch, i.e., the spacing between adjacent tubes in a
            single row.
        pitch_lon:
            Longitudinal pitch, i.e., the spacing between tubes in adjacent
            rows.
        d_o:
            Outside diameter of the tubes.
        d_i:
            Inside diameter of the tubes.
        t_fin:
            Fin thickness.
        fin_density:
            Number of fins per unit length.
        k_fin: optional
            Thermal conductivity of fin material. The default value applies to
            aluminum.
        d_r: optional
            The root or collar diameter of the fins (i.e., the sum of the
            outside diameter and two times the fin thickness). This is ignored
            by default; the root diameter is taken to be equal to the outside
            tube diameter (i.e., fins without a collar).
        num_circuits: optional
            The number of parallel water circuits. By default, the number of
            parallel water circuits is taken to be equal to the number of tubes
            in the first row.
        """
        self.core = HeatExchangerCore(
            width, height, num_rows, pitch_trv,
            pitch_lon, d_o, d_i, t_fin, fin_density,
            k_fin, d_r, num_circuits
        )
        self.air_in: HumidAir | None = None
        self.water_in: FluidState | None = None
        self.air_m_dot: Quantity | None = None
        self.water_m_dot: Quantity | None = None
        self.air_surface_condition: str | None = None
        self.Q_dot_max: Quantity | None = None

    def set_operating_conditions(
        self,
        air_in: HumidAir,
        water_in: FluidState,
        air_m_dot: Quantity,
        water_m_dot: Quantity,
    ) -> None:
        """Sets the operating conditions of the heat exchanger.

        Parameters
        ----------
        air_in:
            State of humid air entering the heat exchanger.
        water_in:
            State of water entering the heat exchanger.
        air_m_dot:
            Mass flow rate of air through the heat exchanger.
        water_m_dot:
            Mass flow rate of water through the heat exchanger.
        """
        self.air_in = air_in
        self.water_in = water_in
        self.air_m_dot = air_m_dot
        self.water_m_dot = water_m_dot
        air_out_sat = HumidAir(Tdb=self.water_in.T, RH=Q_(100, 'pct'))
        self.Q_dot_max = self.air_m_dot * (self.air_in.h - air_out_sat.h)

    def rate(
        self,
        eps_ini: float = 0.5,
        i_max: int = 5,
        tol: Quantity = Q_(0.1, 'K')
    ) -> dict[str, Any]:
        """Determines the performance of the heat exchanger under the current
        set of operating conditions.

        The method uses an iterative solving technique, starting with an initial
        guess for the heat exchanger effectiveness. Usually the heat exchanger
        effectiveness of an air-cooling coil lies somewhere between 0.5 and
        0.75.

        Parameters
        ----------
        eps_ini:
            An initial guess for the heat exchanger effectiveness.
        i_max:
            Maximum number of iterations to find an acceptable solution.
        tol:
            Stop criterion for the iteration. It is the allowable deviation
            between the last and previous value of the calculated air leaving
            and water leaving temperature.

        Returns
        -------
        Dictionary with keys (strings):
        'air_out': HumidAir
            State of humid air leaving the coil.
        'water_out': Fluid('Water')
            State of water leaving the coil.
        'Q_dot': Quantity
            Heat transfer rate from air to water.
        'eps': Quantity
            Heat exchanger effectiveness.
        'dP_air': Quantity
            Air-side drop of pressure.
        'dP_water': Quantity
            Water-side drop of pressure.
        'air_surface_condition': str
            The condition of the air-side heat transfer surface: either fully
            wet ('wet'), or partially wet ('partially_wet'), or fully dry.
        """
        Q_dot = eps_ini * self.Q_dot_max
        water_out = self.__init_guess_water_out(Q_dot)
        air_out = self.__init_guess_air_out(Q_dot)
        i = 0
        while i < i_max:
            air_mean = self.__air_mean(air_out, water_out)
            water_mean = self.__water_mean(water_out)
            hex_props = self.core.hex_properties(
                air_mean, self.air_in, air_out,
                water_mean,
                self.air_m_dot, self.water_m_dot
            )
            A_ext = self.core.geometry.external.A_tot.to('m ** 2')
            A_int = self.core.geometry.internal.A_tot.to('m ** 2')

            air_surf_condition = self.__get_air_side_surface_condition(
                hex_props, A_ext, A_int,
                air_mean, water_mean,
                water_out, air_out
            )
            if air_surf_condition == 'wet':
                air_out_new, T_water_out_new, heX = self.__hex_wet(
                    water_out, water_mean,
                    hex_props,
                    A_ext, A_int
                )
            elif air_surf_condition == 'dry':
                air_out_new, T_water_out_new, heX = self.__hex_dry(
                    water_mean, air_mean, hex_props
                )
            else:
                # air_surf_condition == 'partially_wet'
                # If the air-side heat transfer surface is partially dry and
                # partially wet, we take the higher heat transfer rate for a
                # totally wet and a totally dry coil to provide a good estimate
                # for a partially wet coil (J.W. Mitchell & J. E. Braun,
                # "Principles of Heating, Ventilation, and Air Conditioning in
                # Buildings" (2013), John Wiley & Sons).
                air_out_new_1, T_water_out_new_1, heX_1 = self.__hex_wet(
                    water_out, water_mean,
                    hex_props,
                    A_ext, A_int
                )
                air_out_new_2, T_water_out_new_2, heX_2 = self.__hex_dry(
                    water_mean, air_mean, hex_props
                )
                if heX_1.Q >= heX_2.Q:
                    air_out_new = air_out_new_1
                    T_water_out_new = T_water_out_new_1
                    heX = heX_1
                else:
                    air_out_new = air_out_new_2
                    T_water_out_new = T_water_out_new_2
                    heX = heX_2

            dP_air = hex_props['dP_ext']
            P_air_out_new = self.air_in.P - dP_air
            air_out_new = HumidAir(
                Tdb=air_out_new.Tdb,
                W=air_out_new.W,
                P=P_air_out_new
            )

            dP_water = hex_props['dP_int']
            P_water_out_new = self.water_in.P - dP_water
            if P_water_out_new.magnitude < 0:
                raise ValueError(
                    "Water-side pressure drop is larger than inlet water "
                    "pressure. A higher water-system pressure is required."
                )
            water_out_new = Water(
                T=T_water_out_new,
                P=P_water_out_new
            )
            dev_a = abs(air_out_new.Tdb.to('K') - air_out.Tdb.to('K'))
            dev_w = abs(water_out_new.T.to('K') - water_out.T.to('K'))
            if (dev_a <= tol.to('K')) and (dev_w <= tol.to('K')):
                return {
                    'air_out': air_out,
                    'water_out': water_out,
                    'Q_dot': heX.Q,
                    'eps': Q_(heX.eps, 'frac'),
                    'dP_air': dP_air,
                    'dP_water': dP_water,
                    'air_surface_condition': air_surf_condition
                }
            air_out = air_out_new
            water_out = water_out_new
            i += 1
        else:
            raise ValueError(
                f"No solution found within tolerance {tol:~P} "
                f"after {i_max} iterations."
            )

    def __init_guess_water_out(self, Q_dot: Quantity) -> FluidState:
        C_water = self.water_in.cp * self.water_m_dot
        T_water_out = self.water_in.T + Q_dot / C_water
        water_out = Water(T=T_water_out, P=self.water_in.P)
        return water_out

    def __init_guess_air_out(self, Q_dot: Quantity) -> HumidAir:
        h_air_out = self.air_in.h - Q_dot / self.air_m_dot
        air_out = HumidAir(h=h_air_out, RH=Q_(100, 'pct'))  # we guess saturated outlet air
        return air_out

    def __air_mean(self, air_out: HumidAir, water_out: FluidState):
        air_in_sat = HumidAir(Tdb=water_out.T, RH=Q_(100, 'pct'))
        air_out_sat = HumidAir(Tdb=self.water_in.T, RH=Q_(100, 'pct'))
        lmed = self.__lmed(air_out, air_in_sat, air_out_sat)
        h_air_avg_sat = (air_in_sat.h + air_out_sat.h) / 2
        h_air_avg = h_air_avg_sat + lmed
        lmtd = self.__lmtd(air_out, water_out)
        T_water_avg = (self.water_in.T + water_out.T) / 2
        T_air_avg = T_water_avg + lmtd
        air_mean = HumidAir(Tdb=T_air_avg, h=h_air_avg)
        return air_mean

    def __lmed(
        self,
        air_out: HumidAir,
        air_in_sat: HumidAir,
        air_out_sat: HumidAir
    ):
        dh_in = self.air_in.h - air_in_sat.h
        dh_out = air_out.h - air_out_sat.h
        dh_max = max(dh_in, dh_out)
        dh_min = min(dh_in, dh_out)
        if dh_min.m <= 0.0:
            dh_min = Q_(1.e-12, 'kJ / kg')
        lmed = (dh_max - dh_min) / np.log(dh_max / dh_min)
        return lmed

    def __lmtd(self, air_out: HumidAir, water_out: FluidState):
        dT_in = self.air_in.Tdb - water_out.T
        dT_out = air_out.Tdb - self.water_in.T
        dT_max = max(dT_in, dT_out)
        dT_min = min(dT_in, dT_out)
        if dT_min <= 0.0:
            dT_min = Q_(1.e-12, 'K')
        lmtd = (dT_max - dT_min) / np.log(dT_max / dT_min)
        return lmtd

    def __water_mean(self, water_out: FluidState):
        T_water_avg = (self.water_in.T + water_out.T) / 2
        P_water_avg = (self.water_in.P + water_out.P) / 2
        water_mean = Water(T=T_water_avg, P=P_water_avg)
        return water_mean

    def __get_air_side_surface_condition(
        self,
        hex_props: dict[str, Quantity | float],
        A_ext: Quantity,
        A_int: Quantity,
        air_mean: HumidAir,
        water_mean: FluidState,
        water_out: FluidState,
        air_out: HumidAir
    ) -> str:
        """Condensation of air starts where the surface temperature equals the
        dew point temperature of the entering air. In this function, the air
        temperature and water temperature are determined that correspond with a
        surface temperature equal to the dew point temperature of the entering
        air. Then, we can also calculate the enthalpy that belongs to this air
        temperature. If this enthalpy value is greater than the enthalpy of the
        entering air, the surface is fully wet. If this enthalpy value is
        between the enthalpy of the entering air and the leaving air, the
        surface is partially dry and partially wet. Otherwise, if this enthalpy
        value is less than the enthalpy of the leaving air, the surface is
        fully dry. (W.F. Stoecker & J.W. Jones, "Refrigeration & Air
        Conditioning" Second Ed. (1982), McGraw-Hill)
        """
        h_ext = hex_props['h_ext'].to('W / (m**2 * K)').magnitude
        h_int = hex_props['h_int'].to('W / (m**2 * K)').magnitude
        A_ext = A_ext.to('m**2').magnitude
        A_int = A_int.to('m**2').magnitude
        Tdp_air_in = self.air_in.Tdp.to('K').magnitude
        air_m_dot = self.air_m_dot.to('kg / s').magnitude
        water_m_dot = self.water_m_dot.to('kg / s').magnitude
        cp_air = air_mean.cp.to('J / (kg * K)').magnitude
        cp_water = water_mean.cp.to('J / (kg * K)').magnitude
        T_water_out = water_out.T.to('K').magnitude
        T_air_in = self.air_in.Tdb.to('K').magnitude
        air_dp = HumidAir(Tdb=self.air_in.Tdp, RH=Q_(100, 'pct'))
        h_air_dp = air_dp.h.to('J / kg').magnitude
        h_air_in = self.air_in.h.to('J / kg').magnitude
        h_air_out = air_out.h.to('J / kg').magnitude
        a11 = h_ext * (A_ext / A_int)
        a12 = h_int
        b11 = (h_int + h_ext * (A_ext / A_int)) * Tdp_air_in
        a21 = -air_m_dot * cp_air
        a22 = water_m_dot * cp_water
        b21 = water_m_dot * cp_water * T_water_out - air_m_dot * cp_air * T_air_in
        A = np.array([
            [a11, a12],
            [a21, a22]
        ])
        B = np.array([b11, b21])
        X = np.linalg.solve(A, B)
        T_air_star, T_water_star = X[0], X[1]
        c = h_ext * A_ext / (cp_air * A_int)
        h_air_star = (h_int * (Tdp_air_in - T_water_star) + c * h_air_dp) / c
        if h_air_star >= h_air_in:
            return 'wet'
        elif h_air_in > h_air_star > h_air_out:
            return 'partially_wet'
        else:
            # h_air_star <= h_air_out
            return 'dry'

    def __hex_wet(
        self,
        water_out: FluidState,
        water_mean: FluidState,
        hex_props: dict,
        A_ext: Quantity,
        A_int: Quantity
    ) -> tuple[HumidAir, Quantity, wet.CounterFlowHeatExchanger]:
        heX = wet.CounterFlowHeatExchanger(
            m_dot_r=self.water_m_dot,
            m_dot_a=self.air_m_dot,
            T_r_in=self.water_in.T,
            T_r_out=None,
            T_r_out_ini=water_out.T,
            P_r=water_mean.P,
            refrigerant=Water,
            air_in=self.air_in,
            h_ext=hex_props['h_ext'],
            h_int=hex_props['h_int'],
            eta_surf_wet=hex_props['eta_surf'],
            A_ext_to_A_int=A_ext.m / A_int.m,
            A_ext=A_ext
        )
        air_out_new = heX.air_out
        T_water_out_new = heX.T_r_out
        return air_out_new, T_water_out_new, heX

    def __hex_dry(
        self,
        water_mean: FluidState,
        air_mean: HumidAir,
        hex_props: dict
    ) -> tuple[HumidAir, Quantity, dry.CounterFlowHeatExchanger]:
        heX = dry.CounterFlowHeatExchanger(
            C_cold=water_mean.cp * self.water_m_dot,
            C_hot=air_mean.cp * self.air_m_dot,
            T_cold_in=self.water_in.T,
            T_hot_in=self.air_in.Tdb,
            UA=hex_props['UA']
        )
        air_out_new = HumidAir(Tdb=heX.T_hot_out, W=self.air_in.W)
        T_water_out_new = heX.T_cold_out
        return air_out_new, T_water_out_new, heX
