from math import pi as PI
from abc import ABC, abstractmethod
from hvac import Quantity
from hvac.fluids import Fluid, FluidState, HumidAir
from hvac.heat_transfer.forced_convection import internal_flow
from hvac.heat_transfer.condensation import flow_condensation
from hvac.heat_transfer.boiling import flow_boiling
from hvac.heat_transfer.finned_surface.fins import Fin, PlainContinuousFin
from hvac.heat_transfer.heat_exchanger.misc import correct_nusselt_number, correct_friction_factor
from .. import correlations
from . import geometry


Q_ = Quantity
ht_correlation = correlations.plain_flat_fin.heat_transfer.j_Wang_and_Chi
ff_correlation = correlations.plain_flat_fin.friction_factor.f_Wang_and_Chi


class PlainFinTubeHeatExchangerCore:

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
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
        condensing: bool = False,
        boiling: bool = False,
        **kwargs
    ) -> None:
        if condensing:
            self.internal_surface = _InternalCondensingHeatTransferSurface(
                parent=self,
                geometry=geometry.TubeBankInside(
                    S_t, S_l,
                    D_i,
                    None,
                    L1, None, L3
                )
            )
        elif boiling:
            self.internal_surface = _InternalBoilingHeatTransferSurface(
                parent=self,
                geometry=geometry.TubeBankInside(
                    S_t, S_l,
                    D_i,
                    None,
                    L1, None, L3
                )
            )
        else:
            self.internal_surface = _InternalSinglePhaseHeatTransferSurface(
                parent=self,
                geometry=geometry.TubeBankInside(
                    S_t, S_l,
                    D_i,
                    None,
                    L1, None, L3
                ),
                tube_length=kwargs.get('int_tube_length')
            )
        self.external_surface = _ExternalSinglePhaseHeatTransferSurface(
            parent=self,
            geometry=geometry.PlainFinStaggeredTBO(
                S_t, S_l,
                D_o, D_o,
                t_f, N_f,
                None, None,
                L1, None, L3
            ),
            k_f=k_fin
        )
        self._T_wall: Quantity | None = None

    @property
    def L2(self) -> Quantity:
        return self.internal_surface.geometry.L2

    @L2.setter
    def L2(self, val: Quantity) -> None:
        self.internal_surface.geometry.L2 = val
        self.external_surface.geometry.L2 = val
        self.internal_surface.tube.L = val

    @property
    def m_dot_int(self) -> Quantity:
        return self.internal_surface.m_dot

    @m_dot_int.setter
    def m_dot_int(self, val: Quantity) -> Quantity:
        self.internal_surface.m_dot = val

    @property
    def m_dot_ext(self) -> Quantity:
        return self.external_surface.m_dot

    @m_dot_ext.setter
    def m_dot_ext(self, val: Quantity) -> Quantity:
        self.external_surface.m_dot = val

    @property
    def UA(self) -> Quantity:
        """Returns overall heat transfer conductance of the heat exchanger."""
        R_int = self.internal_surface.R
        R_foul_int = Q_(0.0, 'K / W')  # we ignore fouling on internal side
        R_cond = Q_(0.0, 'K / W')  # we ignore resistance to conduction
        R_ext = self.external_surface.R
        R_foul_ext = Q_(0.0, 'K / W')  # we ignore fouling on external side
        R_tot = R_int + R_foul_int + R_cond + R_ext + R_foul_ext
        UA = 1 / R_tot
        return UA

    def get_wall_temperature(
        self,
        h_int: Quantity,
        h_ext: Quantity,
        eta_ext: float
    ) -> Quantity:
        """Returns the average wall temperature of the tubes."""
        A_int = self.internal_surface.geometry.A
        R_int = 1 / (h_int * A_int)
        A_ext = self.external_surface.geometry.A
        R_ext = 1 / (eta_ext * h_ext * A_ext)
        T_int_m = self.internal_surface.fluid_mean.T.to('K')
        if isinstance(self.external_surface.fluid_mean, HumidAir):
            T_ext_m = self.external_surface.fluid_mean.Tdb
        else:
            T_ext_m = self.external_surface.fluid_mean.T
        n = T_int_m / R_int + T_ext_m / R_ext
        d = 1 / R_int + 1 / R_ext
        T_w = n / d
        return T_w.to('K')

    @property
    def T_wall(self) -> Quantity:
        if self._T_wall is None:
            if isinstance(self.ext.fluid_mean, HumidAir):
                T_fluid_ext = self.ext.fluid_mean.Tdb.to('K')
            else:
                T_fluid_ext = self.ext.fluid_mean.T.to('K')
            T_fluid_min = min(
                self.int.fluid_mean.T.to('K'),
                T_fluid_ext
            )
            T_fluid_max = max(
                self.int.fluid_mean.T.to('K'),
                T_fluid_ext
            )
            T_w = (T_fluid_min + T_fluid_max) / 2
            i_max = 5
            i = 0
            tol_T_w = Q_(0.01, 'K')
            while i < i_max:
                if isinstance(self.int, _InternalBoilingHeatTransferSurface):
                    h_int = self.int.h
                else:
                    h_int = self.int.get_heat_trf_coeff(T_w)
                h_ext = self.ext.get_heat_trf_coeff(T_w)
                eta_ext = self.ext.get_eta(h_ext)
                T_w_new = self.get_wall_temperature(h_int, h_ext, eta_ext)
                dev_T_w = abs(T_w_new - T_w)
                if dev_T_w <= tol_T_w:
                    self._T_wall = T_w_new
                    break
                T_w = T_w_new
                i += 1
            else:
                raise ValueError('no solution found for wall temperature')
        return self._T_wall

    @property
    def int(self) -> '_InternalSinglePhaseHeatTransferSurface':
        return self.internal_surface

    @property
    def ext(self) -> '_ExternalSinglePhaseHeatTransferSurface':
        return self.external_surface


TTube = (
    internal_flow.CircularTube |
    flow_condensation.HorizontalTube |
    flow_boiling.Tube
)


class _InternalHeatTransferSurface(ABC):

    def __init__(self, parent: PlainFinTubeHeatExchangerCore):
        self.parent = parent
        self.geometry: geometry.TubeBankInside | None = None
        self.tube: TTube | None = None
        self._fluid_mean: FluidState | None = None
        self._m_dot: Quantity | None = None
        self._Q: Quantity | None = None

    @property
    def geo(self) -> geometry.TubeBankInside:
        return self.geometry

    def _get_fluid_mean(self) -> FluidState:
        return self._fluid_mean

    def _set_fluid_mean(self, val: FluidState) -> None:
        self._fluid_mean = val
        self.parent._T_wall = None

    fluid_mean = property(_get_fluid_mean, _set_fluid_mean)

    def _get_m_dot(self) -> Quantity:
        return self._m_dot

    def _set_m_dot(self, val: Quantity) -> Quantity:
        self._m_dot = val.to('kg / s')

    m_dot = property(_get_m_dot, _set_m_dot)

    def _get_Q(self) -> Quantity:
        return self._Q

    def _set_Q(self, val: Quantity) -> None:
        self._Q = val

    Q = property(_get_Q, _set_Q)

    @property
    @abstractmethod
    def h(self) -> Quantity:
        """Returns the average convection heat transfer coefficient on the
        internal side of the heat exchanger.
        """
        ...

    @abstractmethod
    def get_heat_trf_coeff(self, *args, **kwargs) -> Quantity:
        """Returns the convection heat transfer coefficient on the refrigerant
        side.
        """
        ...

    @property
    def R(self) -> Quantity:
        """Returns thermal resistance at internal surface."""
        h_int = self.h.to('W / (m ** 2 * K)')
        A_int = self.geometry.A.to('m ** 2')
        R_int = 1 / (h_int * A_int)  # K / W
        return R_int

    def _get_m_dot_tube(self) -> Quantity:
        """Returns the mass flow rate of internal fluid in a single tube
        assuming that the fluid is fed into the first row of the heat exchanger
        and leaves the heat exchanger through the last row.
        """
        G_rfg = self.m_dot / (self.geometry.A_o / self.geometry.N_r)
        m_dot_tube = G_rfg * (PI * self.geometry.D_i ** 2 / 4)
        return m_dot_tube


class _InternalSinglePhaseHeatTransferSurface(_InternalHeatTransferSurface):
    
    def __init__(
        self,
        parent: PlainFinTubeHeatExchangerCore,
        geometry: geometry.TubeBankInside,
        tube_length: Quantity | None = None
    ) -> None:
        super().__init__(parent)
        self.geometry = geometry
        self.tube = internal_flow.CircularTube(  # assumption: smooth tube (e = 0.0 mm)
            Di=self.geometry.D_i,
            L=tube_length or self.geometry.L2 or Q_(float('nan'), 'm'),
            fluid=None
        )

    @property
    def h(self) -> Quantity:
        """Returns the average convection heat transfer coefficient on the
        internal side of the heat exchanger.
        """
        return self.get_heat_trf_coeff()

    def get_heat_trf_coeff(self, *args, **kwargs) -> Quantity:
        """Returns the convection heat transfer coefficient on the refrigerant
        side.
        """
        self.tube.fluid = self._fluid_mean
        self.tube.m_dot = self._get_m_dot_tube()
        h_int = self.tube.avg_heat_transfer_coefficient()
        if isinstance(h_int, tuple):
            h_int = h_int[0]  # h_int for laminar flow with constant heat flux
        return h_int.to('W / (m ** 2 * K)')


class _InternalCondensingHeatTransferSurface(_InternalHeatTransferSurface):

    def __init__(
        self,
        parent: PlainFinTubeHeatExchangerCore,
        geometry: geometry.TubeBankInside
    ) -> None:
        super().__init__(parent)
        self.geometry = geometry
        self.tube = flow_condensation.HorizontalTube(
            D=self.geometry.D_h,
            fluid=Fluid('R134a')
        )

    def _set_fluid_mean(self, val: FluidState) -> None:
        super().fluid_mean = val
        self.tube.fluid = val.fluid

    @property
    def h(self) -> Quantity:
        """Returns the average convection heat transfer coefficient on the
        internal side of the heat exchanger.
        """
        h_int = self.get_heat_trf_coeff(self.parent.T_wall)
        return h_int.to('W / (m ** 2 * K)')

    def get_heat_trf_coeff(self, T_wall: Quantity) -> Quantity:
        """Returns the convection heat transfer coefficient on the refrigerant
        side based on the given wall temperature.
        """
        self.tube.m_dot = self._get_m_dot_tube()
        h_int = self.tube.heat_trf_coeff(
            x=self._fluid_mean.x,
            T_sat=self.fluid_mean.T,
            T_surf=T_wall
        )
        return h_int.to('W / (m ** 2 * K)')


class _InternalBoilingHeatTransferSurface(_InternalHeatTransferSurface):

    def __init__(
        self,
        parent: PlainFinTubeHeatExchangerCore,
        geometry: geometry.TubeBankInside
    ) -> None:
        super().__init__(parent)
        self.geometry = geometry
        self.tube = flow_boiling.Tube(
            Dh=self.geometry.D_h,
            A=self.geometry.D_i ** 2 / 4,
            fluid=Fluid('R134a'),
            tube_orient='horizontal'
        )

    def _set_fluid_mean(self, val: FluidState) -> None:
        super().fluid_mean = val
        self.tube.fluid = val.fluid

    @property
    def h(self) -> Quantity:
        """Returns the average convection heat transfer coefficient on the
        internal side of the heat exchanger.
        """
        return self.get_heat_trf_coeff()

    def get_heat_trf_coeff(self) -> Quantity:
        """Returns the convection heat transfer coefficient on the refrigerant
        side based on the given wall temperature.
        """
        self.tube.m_dot = self._get_m_dot_tube()
        h_int = self.tube.heat_trf_coeff(
            x=self._fluid_mean.x,
            T_sat=self.fluid_mean.T,
            q_s=self._Q / self.geometry.A
        )
        return h_int.to('W / (m ** 2 * K)')


class _ExternalHeatTransferSurface(ABC):

    def __init__(self, parent: PlainFinTubeHeatExchangerCore) -> None:
        self.parent = parent
        self.geometry: geometry.PlainFinStaggeredTBO | None = None
        self.k_f: Quantity | None = None
        self._fluid_in: FluidState | HumidAir | None = None
        self._fluid_out: FluidState | HumidAir | None = None
        self._fluid_mean: FluidState | HumidAir | None = None
        self._m_dot: Quantity | None = None

    @property
    def geo(self) -> geometry.PlainFinStaggeredTBO:
        return self.geometry

    def _get_fluid_in(self) -> FluidState | HumidAir:
        return self._fluid_in

    def _set_fluid_in(self, val: FluidState | HumidAir) -> None:
        self._fluid_in = val
        self.parent._T_wall = None

    fluid_in = property(_get_fluid_in, _set_fluid_in)

    def _get_fluid_out(self) -> FluidState | HumidAir:
        return self._fluid_out

    def _set_fluid_out(self, val: FluidState | HumidAir) -> None:
        self._fluid_out = val

    fluid_out = property(_get_fluid_out, _set_fluid_out)

    def _get_fluid_mean(self) -> FluidState | HumidAir:
        return self._fluid_mean

    def _set_fluid_mean(self, val: FluidState | HumidAir) -> None:
        self._fluid_mean = val
        self.parent._T_wall = None

    fluid_mean = property(_get_fluid_mean, _set_fluid_mean)

    def _get_m_dot(self) -> Quantity:
        return self._m_dot

    def _set_m_dot(self, val: Quantity) -> Quantity:
        self._m_dot = val.to('kg / s')

    m_dot = property(_get_m_dot, _set_m_dot)

    @property
    def h(self) -> Quantity:
        """Returns the average convection heat transfer coefficient on the
        external side of the heat exchanger without correction for variable
        fluid properties.
        """
        h = self.get_heat_trf_coeff(self.parent.T_wall)
        return h.to('W / (m ** 2 * K)')

    @property
    def R(self) -> Quantity:
        """Returns thermal resistance at the external surface with correction
        for variable fluid properties."""
        h_ext = self.h.to('W / (m ** 2 * K)')
        A_ext = self.geometry.A.to('m ** 2')
        eta_ext = self.get_eta(h_ext)
        R_ext = 1 / (eta_ext * h_ext * A_ext)  # K / W
        return R_ext

    def get_eta(self, h: Quantity) -> float:
        """Returns overall heat transfer efficiency of the finned, external
        surface.
        """
        fin = self._get_fin()
        fin.h_avg = h
        eta_fin = fin.efficiency
        A_f = self.geometry.A_f  # m ** 2
        A_ext = self.geometry.A  # m ** 2
        r = (A_f / A_ext).to('m ** 2 / m ** 2').m
        eta = 1 - r * (1 - eta_fin)
        return eta

    @property
    def eta(self) -> float:
        return self.get_eta(self.h)

    def _get_mass_velocity(self) -> Quantity:
        """Returns the mass velocity of external fluid."""
        A_o_ext = self.geometry.A_o.to('m ** 2')
        m_dot_air = self.m_dot.to('kg / s')
        G_air = m_dot_air / A_o_ext
        return G_air

    def _get_reynolds_number(self, D: Quantity) -> float:
        """Returns the Reynolds number of external fluid based on diameter `D`
        (either the hydraulic diameter, or the outside diameter of the tubes).
        """
        G_ext = self._get_mass_velocity()
        mu = self._fluid_mean.mu.to('Pa * s')
        D.ito('m')
        Re_ext = (G_ext * D / mu).m
        return Re_ext

    def _get_prandtl_number(self) -> float:
        """Returns the Prandtl number of the external fluid."""
        rho = self._fluid_mean.rho.to('kg / m ** 3').m
        mu = self._fluid_mean.mu.to('Pa * s').m
        k = self._fluid_mean.k.to('W / (m * K)').m
        cp = self._fluid_mean.cp.to('J / (kg * K)').m
        nu = mu / rho
        alpha = k / (rho * cp)
        Pr = nu / alpha
        return Pr

    @staticmethod
    def _get_nusselt_number(j: float, Re_D_h: float, Pr: float) -> float:
        """Returns the Nusselt number of the external fluid."""
        Nu = j * Re_D_h * Pr ** (1 / 3)
        return Nu

    def _get_fin(self) -> Fin:
        """Configures and returns the fin on the external, finned surface."""
        M = L = self.geometry.S_t / 2
        fin = PlainContinuousFin(
            r_i=self.geometry.D_r / 2,
            t=self.geometry.t_f,
            k=self.k_f,
            M=M,
            L=L
        )
        return fin

    def _get_flow_regime(self) -> str:
        """Returns flow regime of external fluid."""
        Re_D_h = self._get_reynolds_number(D=self.geometry.D_h)
        if Re_D_h > 2300:
            return 'turbulent'
        else:
            return 'laminar'

    def _correct_ext_heat_trf_coeff(
        self,
        Nu: float,
        T_w: Quantity
    ) -> Quantity:
        """Correct external heat transfer coefficient for variable fluid
        properties effects.
        """
        D_h = self.geometry.D_h
        Nu = correct_nusselt_number(
            Nu_cp=Nu,
            T_w=T_w.to('K'),
            flow_regime=self._get_flow_regime(),
            thermal_regime='cooling',
            fluid=self._fluid_mean
        )
        h_corr = Nu * self._fluid_mean.k / D_h
        return h_corr.to('W / (m ** 2 * K)')

    @abstractmethod
    def get_heat_trf_coeff(self, *args, **kwargs) -> Quantity:
        """Returns convection heat transfer coefficient."""
        ...

    @abstractmethod
    def get_pressure_drop(self, *args, **kwargs) -> Quantity:
        """Returns the pressure drop of the external fluid across the
        surface.
        """
        ...

    @property
    def dP(self) -> Quantity:
        return self.get_pressure_drop(self._fluid_in, self._fluid_out)


class _ExternalSinglePhaseHeatTransferSurface(_ExternalHeatTransferSurface):
    
    def __init__(
        self,
        parent: PlainFinTubeHeatExchangerCore,
        geometry: geometry.PlainFinStaggeredTBO,
        k_f: Quantity = Q_(237.0, 'W / (m * K)')
    ) -> None:
        super().__init__(parent)
        self.geometry = geometry
        self.k_f = k_f

    def get_heat_trf_coeff(self, T_wall: Quantity) -> Quantity:
        """Returns convection heat transfer coefficient based on given
        wall temperature.
        """
        Re_D_o = self._get_reynolds_number(D=self.geometry.D_o)
        Pr = self._get_prandtl_number()
        j = ht_correlation(
            Re_D_r=Re_D_o,
            D_r=self.geometry.D_r.to('m').m,
            D_h=self.geometry.D_h.to('m').m,
            S_t=self.geometry.S_t.to('m').m,
            S_l=self.geometry.S_l.to('m').m,
            N_f=self.geometry.N_f.to('1 / m').m,
            N_r=self.geometry.N_r
        )
        Re_D_h = self._get_reynolds_number(D=self.geometry.D_h)
        Nu = self._get_nusselt_number(j, Re_D_h, Pr)
        h_ext = self._correct_ext_heat_trf_coeff(Nu, T_wall)
        return h_ext.to('W / (m ** 2 * K)')

    def _get_fanning_friction_factor(self) -> float:
        """Returns the Fanning friction factor."""
        f = ff_correlation(
            Re_D_r=self._get_reynolds_number(D=self.geo.D_o),
            D_r=self.geo.D_r.to('m').m,
            S_t=self.geo.S_t.to('m').m,
            S_l=self.geo.S_l.to('m').m,
            N_f=self.geo.N_f.to('1 / m').m,
            N_r=self.geo.N_r
        )
        f = correct_friction_factor(
            f_cp=f,
            T_w=self.parent.T_wall,
            flow_regime=self._get_flow_regime(),
            thermal_regime='cooling',
            fluid=self._fluid_mean
        )
        return f

    def get_pressure_drop(
        self,
        fluid_in: FluidState | HumidAir,
        fluid_out: FluidState | HumidAir
    ) -> Quantity:
        G_ext = self._get_mass_velocity()  # kg / (s * m ** 2)
        D_h = self.geo.D_h.to('m')
        L = self.geo.L2.to('m')
        sigma = self.geo.sigma
        rho_a_in = fluid_in.rho
        rho_a_out = fluid_out.rho
        rho_a_avg = 2 * rho_a_in * rho_a_out / (rho_a_in + rho_a_out)
        f = self._get_fanning_friction_factor()
        k = (G_ext ** 2) / (2 * rho_a_in)
        a = f * (4 * L / D_h) * (rho_a_in / rho_a_avg)
        b = (1 + sigma ** 2) * (rho_a_in / rho_a_out - 1)
        dP_a = k * (a + b)
        return dP_a
