from math import pi as PI
from abc import ABC, abstractmethod
from scipy import optimize
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

# heat transfer correlation for plain fin-tube heat exchanger
ht_correlation = correlations.plain_flat_fin.heat_transfer.j_Wang_and_Chi

# Fanning friction factor correlation for plain fin-tube heat exchanger.
ff_correlation = correlations.plain_flat_fin.friction_factor.f_Wang_and_Chi


class PlainFinTubeHeatExchangerCore:
    """Model of the core of a plain fin-tube heat exchanger, used to calculate
    the overall heat transfer conductance of the heat exchanger core.

    Attributes and properties
    -------------------------
    internal_surface or int:
        Reference to the internal heat transfer surface of the heat exchanger
        core.
    external_surface or ext:
        Reference to the external heat transfer surface of the heat exchanger
        core.
    L2:
        External flow length of the heat exchanger core; depth of heat exchanger
        core.
    m_dot_int:
        Mass flow rate of internal fluid (refrigerant).
    m_dot_ext:
        Mass flow rate of external fluid (air).
    UA:
        Overall heat transfer conductance of heat exchanger core.
    T_wall:
        Average wall temperature between internal and external side of the
        heat exchanger core.
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
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
        condensing: bool = False,
        boiling: bool = False,
        N_r: int | None = None,
        **kwargs
    ) -> None:
        """Creates an instance of `PlainFinTubeHeatExchangerCore`.

        Parameters
        ----------
        L1:
            Width of heat exchanger core; first dimension of frontal area;
            tube length.
        L3:
            Height of heat exchanger core; second dimension of frontal area.
        S_t:
            Lateral or transverse pitch, i.e., the spacing between tubes of the
            same row.
        S_l:
            Longitudinal pitch, i.e., the spacing between tubes of two adjacent
            tube rows.
        D_i:
            Inner diameter of tubes.
        D_o:
            Outer diameter of tubes.
        t_f:
            Thickness of plain fins.
        N_f:
            Fin density, i.e., number of fins per unit length of tube.
        k_fin: optional
            Thermal conductivity of fin material; default value is for aluminum.
        condensing: optional
            If True, the internal fluid (refrigerant) inside the tubes is
            condensing, i.e., the case for the condensing part of an air
            condenser.
        boiling: optional
            If True, the internal fluid (refrigerant) inside the tubes is
            boiling, i.e., the case for the boiling part of air evaporator.
        N_r: optional
            Number of rows.
        """
        # Create the internal heat transfer surface.
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
        # Create the external heat transfer surface.
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
        # Set L2 if the number of rows is specified:
        if N_r is not None:
            self.L2 = N_r * S_l
        # Initialize internal variables.
        self._T_wall: Quantity | None = None

    @property
    def L2(self) -> Quantity:
        """Gets external flow length of the heat exchanger core."""
        return self.internal_surface.geometry.L2

    @L2.setter
    def L2(self, val: Quantity) -> None:
        """Sets external flow length of the heat exchanger core."""
        self.internal_surface.geometry.L2 = val
        self.external_surface.geometry.L2 = val
        self.internal_surface.tube.L = val

    @property
    def m_dot_int(self) -> Quantity:
        """Gets mass flow rate of the internal fluid (refrigerant)."""
        return self.internal_surface.m_dot

    @m_dot_int.setter
    def m_dot_int(self, val: Quantity) -> Quantity:
        """Sets mass flow rate of the internal fluid (refrigerant)."""
        self.internal_surface.m_dot = val

    @property
    def m_dot_ext(self) -> Quantity:
        """Gets the mass flow rate of the external fluid (air)."""
        return self.external_surface.m_dot

    @m_dot_ext.setter
    def m_dot_ext(self, val: Quantity) -> Quantity:
        """Sets the mass flow rate of the external fluid (air)."""
        self.external_surface.m_dot = val

    @property
    def UA(self) -> Quantity:
        """Gets overall heat transfer conductance of heat exchanger."""
        R_int = self.internal_surface.R
        R_foul_int = Q_(0.0, 'K / W')  # we ignore fouling on the internal side
        R_cond = Q_(0.0, 'K / W')  # we ignore resistance to conduction
        R_ext = self.external_surface.R
        R_foul_ext = Q_(0.0, 'K / W')  # we ignore fouling on the external side
        R_tot = R_int + R_foul_int + R_cond + R_ext + R_foul_ext
        UA = 1 / R_tot
        return UA

    def _get_wall_temperature(
        self,
        h_int: Quantity,
        h_ext: Quantity,
        eta_ext: float
    ) -> Quantity:
        """Returns the average wall temperature of the tubes.

        Parameters
        ----------
        h_int:
            Internal heat transfer coefficient.
        h_ext:
            External heat transfer coefficient.
        eta_ext:
            Overall efficiency of finned surface.
        """
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
        """Gets average wall temperature of tubes."""
        if self._T_wall is None:
            def _eq(T_w: float) -> float:
                T_w = Q_(T_w, 'K')
                if isinstance(self.int, _InternalBoilingHeatTransferSurface):
                    h_int = self.int.h
                else:
                    # noinspection PyProtectedMember
                    h_int = self.int._get_heat_trf_coeff(T_w)
                h_ext = self.ext.get_heat_trf_coeff(T_w)
                eta_ext = self.ext.get_eta(h_ext)
                T_w_new = self._get_wall_temperature(h_int, h_ext, eta_ext)
                dev = T_w_new.to('K') - T_w.to('K')
                return dev.m

            if isinstance(self.ext.fluid_mean, HumidAir):
                T_fluid_ext = self.ext.fluid_mean.Tdb.to('K')
            else:
                T_fluid_ext = self.ext.fluid_mean.T.to('K')
            T_fluid_min = min(
                self.int.fluid_mean.T.to('K'),
                T_fluid_ext
            ).m + 1.e-3
            T_fluid_max = max(
                self.int.fluid_mean.T.to('K'),
                T_fluid_ext
            ).m - 1.e-3
            try:
                sol = optimize.root_scalar(_eq, bracket=[T_fluid_min, T_fluid_max])
                self._T_wall = Q_(sol.root, 'K')
            except ValueError:
                self._T_wall = Q_((T_fluid_min + T_fluid_max) / 2, 'K')
        return self._T_wall

    @property
    def int(self) -> '_InternalSinglePhaseHeatTransferSurface':
        """Gets internal surface of heat exchanger core."""
        return self.internal_surface

    @property
    def ext(self) -> '_ExternalSinglePhaseHeatTransferSurface':
        """Gets external surface of heat exchanger core."""
        return self.external_surface


TTube = (
    internal_flow.CircularTube |
    flow_condensation.HorizontalTube |
    flow_boiling.Tube
)


class _InternalHeatTransferSurface(ABC):
    """General model for the internal heat transfer surface of the heat exchanger
    core.

    Attributes and properties
    -------------------------
    geometry or geo:
        Geometry of the internal heat transfer surface (instance of
        `TubeBankInside`).
    fluid_mean:
        Mean or average state of the internal fluid.
    tube:
        Reference to instance of annotation type TTube (refers to different
        classes that implement the calculation of the convection heat transfer
        coefficient depending on the internal fluid phase: single-phase, boiling
        or condensing).
    m_dot:
        Mass flow rate of internal fluid.
    Q:
        Heat transfer rate to or from internal fluid.
    h:
        Internal convection heat transfer coefficient.
    R:
        Thermal resistance at internal side of heat exchanger core.
    """
    def __init__(self, parent: PlainFinTubeHeatExchangerCore):
        """Creates the common part for all types of internal heat transfer
        surface.

        Parameters
        ----------
        parent:
            Reference to the `PlainFinTubeHeatExchangerCore` object that created
            this `_InternalHeatTransferSurface` object.
        """
        self.parent = parent
        self.geometry: geometry.TubeBankInside | None = None
        self.tube: TTube | None = None
        self._fluid_mean: FluidState | None = None
        self._m_dot: Quantity | None = None
        self._Q: Quantity | None = None

    @property
    def geo(self) -> geometry.TubeBankInside:
        """Gets the geometry model of the internal heat transfer surface."""
        return self.geometry

    def _get_fluid_mean(self) -> FluidState:
        """Gets the mean fluid state of the internal fluid (refrigerant)."""
        return self._fluid_mean

    def _set_fluid_mean(self, val: FluidState) -> None:
        """Sets the mean fluid state of the internal fluid (refrigerant)."""
        self._fluid_mean = val
        self.parent._T_wall = None

    fluid_mean = property(_get_fluid_mean, _set_fluid_mean)

    def _get_m_dot(self) -> Quantity:
        """Gets the mass flow rate of internal fluid (refrigerant)."""
        return self._m_dot

    def _set_m_dot(self, val: Quantity) -> Quantity:
        """Sets the mass flow rate of internal fluid (refrigerant)."""
        self._m_dot = val.to('kg / s')

    m_dot = property(_get_m_dot, _set_m_dot)

    def _get_Q(self) -> Quantity:
        """Gets the heat transfer rate to or from the internal fluid
        (refrigerant).
        """
        return self._Q

    def _set_Q(self, val: Quantity) -> None:
        """Sets the heat transfer rate to or from the internal fluid
        (refrigerant).
        """
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
    def _get_heat_trf_coeff(self, *args, **kwargs) -> Quantity:
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
        # Mass velocity [kg / (s.m^2)] = mass flow rate of refrigerant divided
        # by the minimum flow area of a single row:
        G_rfg = self.m_dot / (self.geometry.A_o / self.geometry.N_r)
        m_dot_tube = G_rfg * (PI * self.geometry.D_i ** 2 / 4)
        return m_dot_tube


class _InternalSinglePhaseHeatTransferSurface(_InternalHeatTransferSurface):
    """Model of internal heat transfer surface in case of single-phase internal
    fluid (superheated or subcooled refrigerant). Derived from general class
    `_InternalHeatTransferSurface`.

    Attributes and properties
    -------------------------
    See base class `_InternalHeatTransferSurface`.
    """
    def __init__(
        self,
        parent: PlainFinTubeHeatExchangerCore,
        geometry: geometry.TubeBankInside,
        tube_length: Quantity | None = None
    ) -> None:
        """
        Creates a `_InternalSinglePhaseHeatTransferSurface` instance.

        Parameters
        ----------
        parent:
            Reference to the `PlainFinTubeHeatExchangerCore` object that created
            this `_InternalHeatTransferSurface` object.
        geometry:
            Geometry of the internal heat transfer surface (instance of
            `TubeBankInside`).
        tube_length: optional
            The length of a single tube circuit, i.e., the flow length of
            internal fluid (if not specified, it is assumed that it equals the
            external flow length).
        """
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
        return self._get_heat_trf_coeff()

    def _get_heat_trf_coeff(self, *args, **kwargs) -> Quantity:
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
    """Model of internal heat transfer surface in case of condensing internal
    fluid (two-phase refrigerant). Derived from general class
    `_InternalHeatTransferSurface`.

    Attributes and properties
    -------------------------
    See base class `_InternalHeatTransferSurface`.
    """
    def __init__(
        self,
        parent: PlainFinTubeHeatExchangerCore,
        geometry: geometry.TubeBankInside
    ) -> None:
        """
        Creates `_InternalCondensingHeatTransferSurface` instance.

        Parameters
        ----------
        parent:
            Reference to the `PlainFinTubeHeatExchangerCore` object that created
            this `_InternalHeatTransferSurface` object.
        geometry:
            Geometry of the internal heat transfer surface (instance of
            `TubeBankInside`).
        """
        super().__init__(parent)
        self.geometry = geometry
        self.tube = flow_condensation.HorizontalTube(
            D=self.geometry.D_h,
            fluid=Fluid('R134a')
        )

    def _set_fluid_mean(self, val: FluidState) -> None:
        """Sets the mean fluid state of the internal fluid (refrigerant)."""
        super().fluid_mean = val
        self.tube.fluid = val.fluid

    @property
    def h(self) -> Quantity:
        """Returns the average convection heat transfer coefficient on the
        internal side of the heat exchanger.
        """
        h_int = self._get_heat_trf_coeff(self.parent.T_wall)
        return h_int.to('W / (m ** 2 * K)')

    def _get_heat_trf_coeff(self, T_wall: Quantity) -> Quantity:
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
    """Model of internal heat transfer surface in case of boiling internal
    fluid (two-phase refrigerant). Derived from general class
    `_InternalHeatTransferSurface`.

    Attributes and properties
    -------------------------
    See base class `_InternalHeatTransferSurface`.
    """
    def __init__(
        self,
        parent: PlainFinTubeHeatExchangerCore,
        geometry: geometry.TubeBankInside
    ) -> None:
        """Creates `_InternalBoilingHeatTransferSurface` instance.

        Parameters
        ----------
        parent:
            Reference to the `PlainFinTubeHeatExchangerCore` object that created
            this `_InternalHeatTransferSurface` object.
        geometry:
            Geometry of the internal heat transfer surface (instance of
            `TubeBankInside`).
        """
        super().__init__(parent)
        self.geometry = geometry
        self.tube = flow_boiling.Tube(
            Dh=self.geometry.D_h,
            A=self.geometry.D_i ** 2 / 4,
            fluid=Fluid('R134a'),
            tube_orient='horizontal'
        )

    def _set_fluid_mean(self, val: FluidState) -> None:
        """Sets the mean fluid state of the internal fluid (refrigerant)."""
        super().fluid_mean = val
        self.tube.fluid = val.fluid

    @property
    def h(self) -> Quantity:
        """Returns the average convection heat transfer coefficient on the
        internal side of the heat exchanger.
        """
        return self._get_heat_trf_coeff()

    def _get_heat_trf_coeff(self) -> Quantity:
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
    """General model for the external heat transfer surface of the heat exchanger
    core.

    Attributes and properties
    -------------------------
    geometry or geo:
        Geometry of the external heat transfer surface (instance of
        `PlainFinStaggeredTBO`).
    fluid_in:
        External fluid state at heat exchanger's inlet.
    fluid_out:
        External fluid state at heat exchanger's outlet.
    fluid_mean:
        Mean or average state of the external fluid through the heat exchanger
        core.
    m_dot:
        Mass flow rate of external fluid.
    h:
        External convection heat transfer coefficient.
    R:
        Thermal resistance at external side of heat exchanger core.
    eta:
        Overall efficiency of the external, finned surface.
    dP:
        Pressure drop of external fluid.
    """
    def __init__(self, parent: PlainFinTubeHeatExchangerCore) -> None:
        """Creates the common part for all types of external heat transfer
        surface.

        Parameters
        ----------
        parent:
            Reference to the `PlainFinTubeHeatExchangerCore` object that created
            this `_ExternalHeatTransferSurface` object.
        """
        self.parent = parent
        self.geometry: geometry.PlainFinStaggeredTBO | None = None
        self.k_f: Quantity | None = None
        self._fluid_in: FluidState | HumidAir | None = None
        self._fluid_out: FluidState | HumidAir | None = None
        self._fluid_mean: FluidState | HumidAir | None = None
        self._m_dot: Quantity | None = None

    @property
    def geo(self) -> geometry.PlainFinStaggeredTBO:
        """Gets the geometry model of the internal heat transfer surface."""
        return self.geometry

    def _get_fluid_in(self) -> FluidState | HumidAir:
        """Gets the external fluid state at the inlet of the heat exchanger."""
        return self._fluid_in

    def _set_fluid_in(self, val: FluidState | HumidAir) -> None:
        """Sets the external fluid state at the inlet of the heat exchanger."""
        self._fluid_in = val
        self.parent._T_wall = None

    fluid_in = property(_get_fluid_in, _set_fluid_in)

    def _get_fluid_out(self) -> FluidState | HumidAir:
        """Gets the external fluid state at the outlet of the heat exchanger."""
        return self._fluid_out

    def _set_fluid_out(self, val: FluidState | HumidAir) -> None:
        """Sets the external fluid state at the outlet of the heat exchanger."""
        self._fluid_out = val

    fluid_out = property(_get_fluid_out, _set_fluid_out)

    def _get_fluid_mean(self) -> FluidState | HumidAir:
        """Gets the mean fluid state through the heat exchanger."""
        return self._fluid_mean

    def _set_fluid_mean(self, val: FluidState | HumidAir) -> None:
        """Sets the mean fluid state through the heat exchanger."""
        self._fluid_mean = val
        self.parent._T_wall = None

    fluid_mean = property(_get_fluid_mean, _set_fluid_mean)

    def _get_m_dot(self) -> Quantity:
        """Gets the mass flow rate of external fluid."""
        return self._m_dot

    def _set_m_dot(self, val: Quantity) -> Quantity:
        """Sets the mass flow rate of external fluid."""
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
        """Gets overall heat transfer efficiency of the finned, external
        surface.
        """
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
        """Returns the pressure drop of the external fluid across the
        surface.
        """
        return self.get_pressure_drop(self._fluid_in, self._fluid_out)


class _ExternalSinglePhaseHeatTransferSurface(_ExternalHeatTransferSurface):
    """Model of external heat transfer surface in case of single-phase external
    fluid. Derived from general class `_ExternalHeatTransferSurface`.

    Attributes and properties
    -------------------------
    See base class `_ExternalHeatTransferSurface`.
    """
    def __init__(
        self,
        parent: PlainFinTubeHeatExchangerCore,
        geometry: geometry.PlainFinStaggeredTBO,
        k_f: Quantity = Q_(237.0, 'W / (m * K)')
    ) -> None:
        """Creates `_ExternalSinglePhaseHeatTransferSurface` instance.

        Parameters
        ----------
        parent:
            Reference to the `PlainFinTubeHeatExchangerCore` object that created
            this `_ExternalHeatTransferSurface` object.
        geometry:
            Geometry of the external heat transfer surface (instance of
            `PlainFinStaggeredTBO`).
        """
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
        """Returns the pressure drop of the external fluid across the
        surface.
        """
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
