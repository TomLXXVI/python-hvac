from math import pi as PI
from abc import ABC, abstractmethod
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import FluidState
from . import laminar_flow
from . import turbulent_flow

lam_fri = laminar_flow.friction_factor
tur_fri = turbulent_flow.friction_factor
lam_nus = laminar_flow.nusselt_number
tur_nus = turbulent_flow.nusselt_number

Q_ = Quantity

Re_crit = 2300


def hydraulic_diameter(Ac: Quantity, P: Quantity) -> Quantity:
    """Calculates the hydraulic diameter of a tube or duct.

    Parameters
    ----------
    Ac:
        Cross-sectional area of tube or duct.
    P:
        Wetted perimeter of tube or duct.
    """
    Dh = 4 * Ac / P
    return Dh.to('m')


def conv_heat_trf_coeff(Nu: float, k: Quantity, Dh: Quantity) -> Quantity:
    """Returns the convection heat transfer coefficient.

    Parameters
    ----------
    Nu:
        Nusselt number.
    k:
        Heat conduction coefficient of fluid.
    Dh:
        Hydraulic diameter of tube or duct.
    """
    k = k.to('W / (m * K)')
    Dh = Dh.to('m')
    h = Nu * k / Dh
    return h


class Tube(ABC):
    """Common abstract base class from which classes of tubes or ducts with
    different shapes are derived.
    """

    def __init__(
        self,
        Ac: Quantity,
        Dh: Quantity,
        L: Quantity,
        fluid: FluidState | None
    ):
        """Creates the part that is common to all classes of tubes or ducts.

        Parameters
        ----------
        Ac:
            Cross-sectional area of tube.
        Dh:
            Hydraulic diameter of tube.
        L:
            Length of the tube.
        fluid:
            State of the fluid; object that contains the values of the
            thermophysical properties of the fluid that flows through the tube.
            It is assumed that these properties don't change while flowing
            through the tube.
        """
        self.Ac = Ac.to('m ** 2')
        self.Dh = Dh.to('m')
        self.L = L.to('m')
        self.fluid = fluid
        self._m_dot: Quantity | None = None

    @property
    def m_dot(self) -> Quantity:
        """Returns the mass flow rate through the tube."""
        return self._m_dot

    @m_dot.setter
    def m_dot(self, v: Quantity) -> None:
        """Sets the mass flow rate through the tube."""
        self._m_dot = v.to('kg / s')

    @property
    def V_dot(self) -> Quantity:
        """Returns the volume flow rate through the tube."""
        V_dot = self._m_dot / self.fluid.rho
        return V_dot.to('m ** 3 / s')

    @V_dot.setter
    def V_dot(self, v: Quantity) -> None:
        """Sets the volume flow rate through the tube."""
        self._m_dot = self.fluid.rho * v.to('m ** 3 / s')

    def mean_velocity(self) -> Quantity:
        """Returns the mean flow velocity of the fluid through the tube."""
        if self.m_dot is not None:
            u_m = self.m_dot / (self.fluid.rho * self.Ac)
            return u_m.to('m / s')
        else:
            raise ValueError('Flow rate unknown. Set `m_dot` or `V_dot`.')

    def reynolds_number(self) -> float:
        """Returns the Reynolds number of the fluid flow through the tube."""
        u_m = self.mean_velocity()
        Re = u_m * self.Dh * self.fluid.rho / self.fluid.mu
        return Re.to('frac').m

    def laminar_hydrodynamic_entry_length(self) -> Quantity:
        """Returns the distance from the inlet where the momentum boundary layer
        converges in case the flow is laminar. Beyond this distance the flow is
        hydrodynamically fully developed, i.e. the axial flow velocity profile
        remains constant. The region close to the inlet is referred to as the
        hydrodynamically developing region.
        """
        Re = self.reynolds_number()
        x_fd_h_lam = 0.063 * Re * self.Dh
        return x_fd_h_lam.to('m')

    def turbulent_hydrodynamic_entry_length(self, u_inf: Quantity) -> Quantity:
        """Returns the distance from the inlet where the flow transitions from
        laminar to turbulent flow before becoming fully developed.
        """
        Re_x_crit = 5.e5
        u_inf = u_inf.to('m / s').m
        rho = self.fluid.rho.to('kg / m ** 3').m
        mu = self.fluid.mu.to('Pa * s').m
        x_crit = Re_x_crit * mu / (rho * u_inf)
        return Q_(x_crit, 'm')

    def prandtl_number(self) -> float:
        """Returns the Prandtl number of the fluid in the tube."""
        nu = self.fluid.mu / self.fluid.rho
        alpha = self.fluid.k / (self.fluid.rho * self.fluid.cp)
        Pr = nu / alpha
        return Pr.to('frac').m

    def laminar_thermal_entry_length(self) -> Quantity:
        """Returns the distance from the inlet where laminar flow becomes
        thermally fully developed, i.e. the convection heat transfer coefficient
        becomes constant. (In the thermally developing region the boundary layer
        is still growing and because its thickness is smaller, its thermal
        resistance is also smaller, which means that in the thermally developing
        region the heat transfer coefficient is larger, but decreasing towards
        its fully developed value.)

        Notes
        -----
        The ratio of the laminar thermal entry length to the laminar
        hydrodynamic length is approximately equal to the Prandtl number.
        """
        Re = self.reynolds_number()
        Pr = self.prandtl_number()
        x_fd_t_lam = 0.063 * Re * Pr * self.Dh
        return x_fd_t_lam.to('m')

    def get_flow_condition(self) -> str:
        """Returns the type of flow of the fluid in the tube, either 'laminar'
        or 'turbulent' depending on the Reynolds number.
        """
        Re = self.reynolds_number()
        if Re < Re_crit:
            return 'laminar'
        else:
            return 'turbulent'

    def dimensionless_tube_length(self, x: Quantity | None = None) -> float:
        """Returns the ratio of the tube length to the product of the Reynolds
        number and hydraulic diameter of the tube. This ratio is used to calculate
        the Graetz number and the friction factor.
        """
        Re = self.reynolds_number()
        if x is None:
            L_tilde = self.L / (Re * self.Dh)
        else:
            L_tilde = x.to('m') / (Re * self.Dh)
        return L_tilde.to('frac').m

    @abstractmethod
    def friction_factor(self) -> float:
        """Returns the average friction factor along the length of the tube.

        In case of laminar flow the friction factor depends on the shape of
        the tube. In case of turbulent flow the friction factor depends on
        the roughness of the internal tube wall surface.
        """
        ...

    @abstractmethod
    def fully_dev_friction_factor(self) -> float:
        """Returns the fully developed friction factor."""
        ...

    def _avg_turbulent_friction_factor(self, e: Quantity) -> float:
        """Returns the average friction factor in case of turbulent flow.

        Parameters
        ----------
        e:
            Absolute roughness of the internal tube wall surface.
        """
        e = e.to('m')
        Dh = self.Dh.to('m')
        e_r = (e / Dh).m
        if 0.0 <= e_r <= 1e-6:
            f_avg = tur_fri.SmoothTube.average_friction_factor(
                Re=self.reynolds_number(),
                Dh=self.Dh,
                L=self.L
            )
        else:
            f_avg = tur_fri.RoughTube.average_friction_factor(
                Re=self.reynolds_number(),
                Dh=self.Dh,
                L=self.L,
                e=e
            )
        return f_avg

    def _fully_dev_turbulent_friction_factor(self, e: Quantity) -> float:
        """Returns the (local) friction factor in case of fully developed
        turbulent flow.

        Parameters
        ----------
        e:
            Absolute roughness of the internal tube wall surface.
        """
        e = e.to('m')
        Dh = self.Dh.to('m')
        e_r = (e / Dh).m
        if 0.0 <= e_r <= 1e-6:
            f_fd = tur_fri.SmoothTube.fully_developed_local_friction_factor(
                Re=self.reynolds_number(),
            )
        else:
            f_fd = tur_fri.RoughTube.fully_developed_local_friction_factor(
                Re=self.reynolds_number(),
                e_r=e_r
            )
        return f_fd

    def pressure_drop(self) -> Quantity:
        """Returns the pressure drop across the length of the tube."""
        f_avg = self.friction_factor()
        u_m = self.mean_velocity()
        rho = self.fluid.rho
        dp = f_avg * (self.L / self.Dh) * (rho * u_m ** 2) / 2
        return dp.to('Pa')

    def graetz_number(self, x: Quantity | None = None):
        """Returns the Graetz number."""
        Pr = self.prandtl_number()
        L_tilde = self.dimensionless_tube_length(x)
        Gz = Pr / L_tilde
        return Gz

    @abstractmethod
    def avg_heat_transfer_coefficient(self, x: Quantity | None = None) -> Quantity | tuple[Quantity, Quantity]:
        """Returns the average convection heat transfer coefficient between
        the inlet of the tube and position `x` along the tube.

        In case the flow is laminar, the value of the heat transfer coefficient
        will depend on the thermal boundary conditions at the wall surface of the
        tube. Two conditions are considered here for circular and rectangular
        tubes: (1) uniform heat flux along the surface of the tube wall, and (2)
        uniform wall surface temperature. Therefore, two values are returned
        from the method. The first value applies to uniform heat flux, and the
        second value applies to uniform wall temperature.
        For annular ducts only one thermal boundary condition is considered here:
        the outer wall of the annulus is (perfectly) insulated, while the inner
        wall has a uniform temperature. So, only one value is returned.

        In the flow is turbulent, there is also only one value to be returned.
        """
        ...

    def avg_nusselt_number(self, x: Quantity | None = None) -> float:
        """Returns the average Nusselt number between the inlet of the tube
        and position `x` along the tube.

        In case of laminar flow the value of the Nusselt number will
        depend on the thermal boundary conditions at the wall surface of the
        tube. Two conditions are considered here for circular and rectangular
        tubes: (1) uniform heat flux along the surface of the tube wall, and (2)
        uniform wall surface temperature. Therefore, two values are returned
        from the method. The first value applies to uniform heat flux, and the
        second value applies to uniform wall temperature.
        For annular ducts only one thermal boundary condition is considered here:
        the outer wall of the annulus is (perfectly) insulated, while the inner
        wall has a uniform temperature. So, only one value is returned.

        In case of turbulent flow there is also only one value to be returned.
        """
        h_avg = self.avg_heat_transfer_coefficient(x)
        if isinstance(h_avg, tuple):
            Nu_avg = []
            for elem in h_avg:
                elem = elem.to('W / (m ** 2 * K)').m
                Dh = self.Dh.to('m').m
                k = self.fluid.k.to('W / (m * K)').m
                Nu_avg.append(elem * Dh / k)
            Nu_avg = tuple(Nu_avg)
        else:
            h_avg = h_avg.to('W / (m ** 2 * K)').m
            Dh = self.Dh.to('m').m
            k = self.fluid.k.to('W / (m * K)').m
            Nu_avg = h_avg * Dh / k
        return Nu_avg

    @abstractmethod
    def fully_dev_heat_transfer_coefficient(self) -> Quantity | tuple[Quantity, Quantity]:
        """Returns the local fully developed convection heat transfer coefficient.

        In case of laminar flow the value of the heat transfer coefficient will
        depend on the thermal boundary conditions at the wall surface of the
        tube. Two conditions are considered here for circular and rectangular
        tubes: (1) uniform heat flux along the surface of the tube wall, and (2)
        uniform wall surface temperature. Therefore, two values are returned
        from the method. The first value applies to uniform heat flux, and the
        second value applies to uniform wall temperature.
        For annular ducts only one thermal boundary condition is considered here:
        the outer wall of the annulus is (perfectly) insulated, while the inner
        wall has a uniform temperature. So, only one value is returned.

        In case of turbulent flow there is also only one value to be returned.
        """
        ...

    def fully_dev_nusselt_number(self) -> float | tuple[float, float]:
        """Returns the fully developed Nusselt number.

        In case of laminar flow the value of the Nusselt number will
        depend on the thermal boundary conditions at the wall surface of the
        tube. Two conditions are considered here for circular and rectangular
        tubes: (1) uniform heat flux along the surface of the tube wall, and (2)
        uniform wall surface temperature. Therefore, two values are returned
        from the method. The first value applies to uniform heat flux, and the
        second value applies to uniform wall temperature.
        For annular ducts only one thermal boundary condition is considered here:
        the outer wall of the annulus is (perfectly) insulated, while the inner
        wall has a uniform temperature. So, only one value is returned.

        In case of turbulent flow there is also only one value to be returned.
        """
        h = self.fully_dev_heat_transfer_coefficient()
        if isinstance(h, tuple):
            Nu = []
            for elem in h:
                elem = elem.to('W / (m ** 2 * K)').m
                Dh = self.Dh.to('m').m
                k = self.fluid.k.to('W / (m * K)').m
                Nu.append(elem * Dh / k)
            Nu = tuple(Nu)
        else:
            h = h.to('W / (m ** 2 * K)').m
            Dh = self.Dh.to('m').m
            k = self.fluid.k.to('W / (m * K)').m
            Nu = h * Dh / k
        return Nu


class CircularTube(Tube):
    """Derived class of the common `Tube` class that represents tubes with a
    circular cross-section.
    """

    def __init__(
        self,
        Di: Quantity,
        L: Quantity,
        fluid: FluidState | None,
        e: Quantity = Q_(0.0, 'mm')
    ) -> None:
        """
        Creates a `CircularTube` instance.

        Parameters
        ----------
        Di:
            Internal diameter of the tube.
        L:
            Length of the tube.
        fluid:
            State of the fluid; object that contains the values of the
            thermophysical properties of the fluid that flows through the tube.
            It is assumed that these properties don't change while flowing
            through the tube.
        e:
            Absolute roughness of the internal tube wall surface.
        """
        Ac = PI * Di ** 2 / 4
        super().__init__(Ac, Di, L, fluid)
        self.e = e.to('m')
        self.e_r = (self.e / self.Dh).m

    def friction_factor(self) -> float:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            f_avg = lam_fri.CircularTube.average_friction_factor(
                L_tilde=self.dimensionless_tube_length(),
                Re=self.reynolds_number()
            )
        else:  # if flow_condition == 'turbulent':
            f_avg = self._avg_turbulent_friction_factor(self.e)
        return f_avg

    def fully_dev_friction_factor(self) -> float:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            f_fd = lam_fri.CircularTube.fully_developed_local_friction_factor(
                Re=self.reynolds_number()
            )
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(self.e)
        return f_fd

    def avg_heat_transfer_coefficient(self, x: Quantity | None = None) -> Quantity | tuple[Quantity, Quantity]:
        """See docstring of base class `Tube`."""
        if x is None: x = self.L
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            UHF = lam_nus.CircularTube.UniformHeatFlux
            UWT = lam_nus.CircularTube.UniformWallTemperature
            Nu_avg_uhf = UHF.average_nusselt_number(
                Pr=self.prandtl_number(),
                Gz=self.graetz_number(x)
            )
            Nu_avg_uwt = UWT.average_nusselt_number(
                Pr=self.prandtl_number(),
                Gz=self.graetz_number(x)
            )
            h_avg_uhf = conv_heat_trf_coeff(Nu_avg_uhf, self.fluid.k, self.Dh)
            h_avg_uwt = conv_heat_trf_coeff(Nu_avg_uwt, self.fluid.k, self.Dh)
            return h_avg_uhf, h_avg_uwt
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(e=self.e)
            Nu_avg = tur_nus.average_nusselt_number(
                f_fd=f_fd,
                Re=self.reynolds_number(),
                Pr=self.prandtl_number(),
                Dh=self.Dh,
                L=x
            )
            h_avg = conv_heat_trf_coeff(Nu_avg, self.fluid.k, self.Dh)
            return h_avg
    
    def fully_dev_heat_transfer_coefficient(self) -> Quantity | tuple[Quantity, Quantity]:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            UHF = lam_nus.CircularTube.UniformHeatFlux
            UWT = lam_nus.CircularTube.UniformWallTemperature
            Nu_fd_uhf = UHF.fully_developed_local_nusselt_number()
            Nu_fd_uwt = UWT.fully_developed_local_nusselt_number()
            h_fd_uhf = conv_heat_trf_coeff(Nu_fd_uhf, self.fluid.k, self.Dh)
            h_fd_uwt = conv_heat_trf_coeff(Nu_fd_uwt, self.fluid.k, self.Dh)
            return h_fd_uhf, h_fd_uwt
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(e=self.e)
            Nu_fd = tur_nus.fully_developed_local_nusselt_number(
                f_fd=f_fd,
                Re=self.reynolds_number(),
                Pr=self.prandtl_number()
            )
            h_fd = conv_heat_trf_coeff(Nu_fd, self.fluid.k, self.Dh)
            return h_fd


class RectangularTube(Tube):
    """Derived class of the common `Tube` class that represents tubes with a
    rectangular cross-section.
    """
    def __init__(
        self,
        a: Quantity,
        b: Quantity,
        L: Quantity,
        fluid: FluidState | None,
        e: Quantity = Q_(0.0, 'mm')
    ) -> None:
        """
        Creates a `RectangularTube` instance.

        Parameters
        ----------
        a:
            First internal dimension of the rectangular cross-section.
        b:
            Second internal dimension of the rectangular cross-section.
        L:
            Length of the tube.
        fluid:
            State of the fluid; object that contains the values of the
            thermophysical properties of the fluid that flows through the tube.
            It is assumed that these properties don't change while flowing
            through the tube.
        e:
            Absolute roughness of the internal tube wall surface.
        """
        Ac = a * b
        Dh = hydraulic_diameter(Ac, P=2 * (a + b))
        super().__init__(Ac, Dh, L, fluid)
        self.a = a
        self.b = b
        self.e = e.to('m')
        self.e_r = (self.e / self.Dh).m

    def friction_factor(self) -> float:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            f_avg = lam_fri.RectangularTube.average_friction_factor(
                L_tilde=self.dimensionless_tube_length(),
                Re=self.reynolds_number(),
                a=self.a,
                b=self.b
            )
            return f_avg
        else:  # if flow_condition == 'turbulent':
            f_avg = self._avg_turbulent_friction_factor(self.e)
        return f_avg

    def fully_dev_friction_factor(self) -> float:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            f_fd = lam_fri.RectangularTube.fully_developed_local_friction_factor(
                Re=self.reynolds_number(),
                a=self.a,
                b=self.b
            )
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(self.e)
        return f_fd

    def avg_heat_transfer_coefficient(self, x: Quantity | None = None) -> Quantity | tuple[Quantity, Quantity]:
        """See docstring of base class `Tube`."""
        if x is None: x = self.L
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            UHF = lam_nus.RectangularTube.UniformHeatFlux
            UWT = lam_nus.RectangularTube.UniformWallTemperature
            Nu_avg_uhf = UHF.average_nusselt_number(
                a=self.a,
                b=self.b,
                Pr=self.prandtl_number(),
                Gz=self.graetz_number(x)
            )
            Nu_avg_uwt = UWT.average_nusselt_number(
                a=self.a,
                b=self.b,
                Pr=self.prandtl_number(),
                Gz=self.graetz_number(x)
            )
            h_avg_uhf = conv_heat_trf_coeff(Nu_avg_uhf, self.fluid.k, self.Dh)
            h_avg_uwt = conv_heat_trf_coeff(Nu_avg_uwt, self.fluid.k, self.Dh)
            return h_avg_uhf, h_avg_uwt
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(e=self.e)
            Nu_avg = tur_nus.average_nusselt_number(
                f_fd=f_fd,
                Re=self.reynolds_number(),
                Pr=self.prandtl_number(),
                Dh=self.Dh,
                L=x
            )
            h_avg = conv_heat_trf_coeff(Nu_avg, self.fluid.k, self.Dh)
            return h_avg
    
    def fully_dev_heat_transfer_coefficient(self) -> Quantity | tuple[Quantity, Quantity]:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            UHF = lam_nus.RectangularTube.UniformHeatFlux
            UWT = lam_nus.RectangularTube.UniformWallTemperature
            Nu_fd_uhf = UHF.fully_developed_local_nusselt_number(self.a, self.b)
            Nu_fd_uwt = UWT.fully_developed_local_nusselt_number(self.a, self.b)
            h_fd_uhf = conv_heat_trf_coeff(Nu_fd_uhf, self.fluid.k, self.Dh)
            h_fd_uwt = conv_heat_trf_coeff(Nu_fd_uwt, self.fluid.k, self.Dh)
            return h_fd_uhf, h_fd_uwt
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(e=self.e)
            Nu_fd = tur_nus.fully_developed_local_nusselt_number(
                f_fd=f_fd,
                Re=self.reynolds_number(),
                Pr=self.prandtl_number()
            )
            h_fd = conv_heat_trf_coeff(Nu_fd, self.fluid.k, self.Dh)
            return h_fd


class AnnularDuct(Tube):
    """Derived class of the common `Tube` class that represents tubes with a
    concentric, annular cross-section.
    """
    def __init__(
        self,
        ri: Quantity,
        ro: Quantity,
        L: Quantity,
        fluid: FluidState | None,
        e: Quantity = Q_(0.0, 'mm')
    ) -> None:
        """
        Creates a `AnnularDuct` instance.

        Parameters
        ----------
        ri:
            Radius of inner tube.
        ro:
            Radius of outer tube.
        L:
            Length of the duct.
        fluid:
            State of the fluid; object that contains the values of the
            thermophysical properties of the fluid that flows through the tube.
            It is assumed that these properties don't change while flowing
            through the tube.
        e:
            Absolute roughness of the internal tube wall surface.
        """
        Ac = PI * (ro ** 2 - ri ** 2)
        Dh = hydraulic_diameter(Ac, P=2 * PI * (ri + ro))
        super().__init__(Ac, Dh, L, fluid)
        self.ri = ri
        self.ro = ro
        self.e = e.to('m')
        self.e_r = (self.e / self.Dh).m

    def friction_factor(self) -> float:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            f_avg = lam_fri.AnnularDuct.average_friction_factor(
                L_tilde=self.dimensionless_tube_length(),
                Re=self.reynolds_number(),
                ri=self.ri,
                ro=self.ro
            )
            return f_avg
        else:  # if flow_condition == 'turbulent':
            f_avg = self._avg_turbulent_friction_factor(self.e)
        return f_avg

    def fully_dev_friction_factor(self) -> float:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            f_fd = lam_fri.AnnularDuct.fully_developed_local_friction_factor(
                Re=self.reynolds_number(),
                ri=self.ri,
                ro=self.ro
            )
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(self.e)
        return f_fd

    def avg_heat_transfer_coefficient(self, x: Quantity | None = None) -> Quantity:
        """See docstring of base class `Tube`."""
        if x is None: x = self.L
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            # UHF = lam_nus.AnnularDuct.UniformHeatFlux --> not implemented
            UWT = lam_nus.AnnularDuct.UniformWallTemperature
            # Nu_avg_uhf = UHF.average_nusselt_number(
            #     a=self.a,
            #     b=self.b,
            #     Pr=self.prandtl_number(),
            #     Gz=self.graetz_number(x)
            # )
            Nu_avg_uwt = UWT.average_nusselt_number(
                ri=self.ri,
                ro=self.ro,
                Pr=self.prandtl_number(),
                Gz=self.graetz_number(x)
            )
            # h_avg_uhf = conv_heat_trf_coeff(Nu_avg_uhf, self.fluid.k, self.Dh)
            h_avg_uwt = conv_heat_trf_coeff(Nu_avg_uwt, self.fluid.k, self.Dh)
            return h_avg_uwt
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(e=self.e)
            Nu_avg = tur_nus.average_nusselt_number(
                f_fd=f_fd,
                Re=self.reynolds_number(),
                Pr=self.prandtl_number(),
                Dh=self.Dh,
                L=x
            )
            h_avg = conv_heat_trf_coeff(Nu_avg, self.fluid.k, self.Dh)
            return h_avg
    
    def fully_dev_heat_transfer_coefficient(self) -> Quantity:
        """See docstring of base class `Tube`."""
        flow_condition = self.get_flow_condition()
        if flow_condition == 'laminar':
            UWT = lam_nus.AnnularDuct.UniformWallTemperature
            Nu_fd_uwt = UWT.fully_developed_local_nusselt_number(self.ri, self.ro)
            h_fd_uwt = conv_heat_trf_coeff(Nu_fd_uwt, self.fluid.k, self.Dh)
            return h_fd_uwt
        else:  # if flow_condition == 'turbulent':
            f_fd = self._fully_dev_turbulent_friction_factor(e=self.e)
            Nu_fd = tur_nus.fully_developed_local_nusselt_number(
                f_fd=f_fd,
                Re=self.reynolds_number(),
                Pr=self.prandtl_number()
            )
            h_fd = conv_heat_trf_coeff(Nu_fd, self.fluid.k, self.Dh)
            return h_fd


def find_volume_flow_rate(
    tube: Tube,
    dp: Quantity,
    V_dot_max: Quantity
) -> Quantity:
    """
    Finds the volume flow rate through a tube when the pressure drop across the
    tube is given.

    Parameters
    ----------
    tube:
        Instance of a subclass of `Tube`.
    dp:
        The pressure drop across the tube.
    V_dot_max:
        Maximum possible volume flow rate through the tube that will be used by
        the root finding algorithm for closing the search area. This value may
        be guessed.

    If no solution was found below `V_dot_max, the searching algorithm will try
    to find a solution by shifting the search area. The search area is then
    bounded by the previous upper limit and a new upper limit which is 1/10
    of `V_dot_max` above the previous upper limit. If no solution is found after
    10 trials, a `ValueError` is raised.
    """
    dp_target = dp.to('Pa').m
    V_dot_i = Q_(1.e-9, 'm ** 3 / s').m
    V_dot_f = V_dot_max.to('m ** 3 / s').m
    V_dot_incr = V_dot_f / 10
    max_iter = 10

    def _eq(V_dot: float) -> float:
        tube.V_dot = Q_(V_dot, 'm ** 3 / s')
        dp_ = tube.pressure_drop().to('Pa').m
        return dp_target - dp_

    i = 0
    while i <= max_iter:
        try:
            V = root_scalar(_eq, bracket=[V_dot_i, V_dot_f]).root
        except ValueError:
            if i < max_iter:
                V_dot_i = V_dot_f
                V_dot_f += V_dot_incr
                i += 1
                continue
            else:
                raise ValueError(f'no solution found below {V_dot_f:~P}') from None
        else:
            return Q_(V, 'm ** 3 / s')
