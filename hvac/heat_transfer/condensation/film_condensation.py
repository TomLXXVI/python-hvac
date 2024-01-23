"""
In this module, correlations are implemented for film condensation on a
few different geometries.

The correlations were taken from:
Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.
"""
import warnings
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid
from ..finned_surface import CircularRectangularFin

Q_ = Quantity


class VerticalPlate:
    """Class for calculating the condensate mass or volume flow rate and
    average heat transfer coefficient for film condensation of a fluid on a
    cool vertical flat surface (plate).

    To use this class, instantiate the class with the fixed parameters of the
    problem. Next, set the surface temperature of the cooled plate and the
    saturation temperature of the condensing fluid through the appropriate
    setters `T_surf` and `T_sat`. Then, set an initial guess for the condensate
    mass flow rate through setter `m_dot`. The actual condensate mass flow rate
    is solved by iteration.
    To get the actual condensate mass flow rate, use the `m_dot` getter.
    To get the actual condensate volume flow rate, use the `V_dot` getter.
    To get the average heat transfer coefficient, use the `h_avg` getter.
    """
    g: Quantity = Q_(9.81, 'm / s ** 2')

    @staticmethod
    def _reynolds_number(m_dot: float, W: float, mu_f: float) -> float:
        Re_c = 4 * m_dot / (W * mu_f)
        return Re_c

    @staticmethod
    def _prandtl_number(rho_f: float, mu_f: float, k_f: float, cp_f: float) -> float:
        nu_f = mu_f / rho_f
        alpha_f = k_f / (rho_f * cp_f)
        Pr = nu_f / alpha_f
        return Pr

    @staticmethod
    def _heat_trf_coeff_Re_low(
        Re_c: float,
        k_f: float,
        mu_f: float,
        rho_f: float,
        rho_g: float,
    ) -> float:
        g = VerticalPlate.g.to('m / s ** 2').m
        n = mu_f ** 2
        d = rho_f * (rho_f - rho_g) * g
        z = (n / d) ** (-1 / 3) * k_f
        h_avg = 1.47 * (Re_c ** (-1 / 3)) * z
        return h_avg

    @staticmethod
    def _heat_trf_coeff_Re_medium(
        Re_c: float,
        k_f: float,
        mu_f: float,
        rho_f: float,
        rho_g: float,
    ) -> float:
        g = VerticalPlate.g.to('m / s ** 2').m
        n = mu_f ** 2
        d = rho_f * (rho_f - rho_g) * g
        z = (n / d) ** (-1 / 3) * k_f
        h_avg = (Re_c / (1.08 * (Re_c ** 1.22) - 5.2)) * z
        return h_avg

    @staticmethod
    def _heat_trf_coeff_Re_high(
        Re_c: float,
        k_f: float,
        mu_f: float,
        rho_f: float,
        Pr_f: float,
        rho_g: float,
    ) -> float:
        g = VerticalPlate.g.to('m / s ** 2').m
        n = mu_f ** 2
        d = rho_f * (rho_f - rho_g) * g
        z = (n / d) ** (-1 / 3) * k_f
        h_avg = (Re_c / (8750 + (58 / (Pr_f ** 0.5)) * ((Re_c ** 0.75) - 253))) * z
        return h_avg

    @staticmethod
    def _condensate_mass_flow_rate(
        h_avg: float,
        h_fg: float,
        cp_f: float,
        T_sat: float,
        T_surf: float,
        A: float
    ) -> float:
        h_fg_star = h_fg + 0.68 * cp_f * (T_sat - T_surf)
        # latent heat of vaporization corrected to account for the sensible
        # cooling of the condensate
        m_dot = h_avg * A * (T_sat - T_surf) / h_fg_star
        return m_dot

    def __init__(
        self,
        W: Quantity,
        L: Quantity,
        fluid: Fluid,
        i_max: int = 20,
        tol: Quantity = Q_(1.e-4, 'kg / s')
    ) -> None:
        """Creates a `VerticalPlate` instance.

        Parameters
        ----------
        W: PlainQuantity
            The horizontal dimension (width) of the cooled vertical plate.
        L: PlainQuantity
            The vertical dimension (height) of the cooled vertical plate.
        fluid: Fluid
            The condensing fluid.
        i_max: int, default 20
            Maximum number of iterations to find the condensate mass flow rate.
        tol: float, default 0.001 kg/s
            The tolerated deviation between the current and previous value of
            the condensate mass flow rate at which the iteration is stopped.
        """
        self._W = W.to('m').m
        self._L = L.to('m').m
        self._A = self._W * self._L
        self.fluid = fluid
        self._i_max = i_max
        self._tol = tol.to('kg / s').m
        self._T_surf: float | None = None
        self._T_sat: float | None = None
        self._m_dot: float | None = None
        self._V_dot: float | None = None
        self._h_avg: float | None = None

    @property
    def T_surf(self) -> Quantity:
        """Returns the surface temperature of the cooled plate."""
        return Q_(self._T_surf, 'K')

    @T_surf.setter
    def T_surf(self, v: Quantity) -> None:
        """Set the surface temperature of the cooled plate."""
        self._T_surf = v.to('K').m
        self._h_avg = None

    @property
    def T_sat(self) -> Quantity:
        """Returns the (saturation) temperature of the fluid in contact with
        the cooled plate.
        """
        return Q_(self._T_sat, 'K')

    @T_sat.setter
    def T_sat(self, v: Quantity) -> None:
        """Set the (saturation) temperature of the fluid in contact with
        the cooled plate.
        """
        self._T_sat = v.to('K').m
        self._h_avg = None

    def _solve_by_iteration(self):
        # Find the condensate mass flow rate by iteration.
        # Once a solution is found for the mass flow rate within the range of
        # tolerance, the average heat transfer coefficient is also found.
        T_flm = (self.T_sat + self.T_surf) / 2
        liq_sat = self.fluid(T=T_flm, x=Q_(0.0, 'frac'))
        vap_sat = self.fluid(T=T_flm, x=Q_(1.0, 'frac'))
        mu_f = liq_sat.mu.to('Pa * s').m
        k_f = liq_sat.k.to('W / (m * K)').m
        rho_f = liq_sat.rho.to('kg / m ** 3').m
        cp_f = liq_sat.cp.to('J / (kg * K)').m
        h_f = liq_sat.h.to('J / kg').m
        rho_g = vap_sat.rho.to('kg / m ** 3').m
        h_g = vap_sat.h.to('J / kg').m
        h_fg = h_g - h_f
        i = 0
        self._m_dot = 0.001  # kg/s - initial guess
        while i <= self._i_max:
            Re_c = self._reynolds_number(self._m_dot, self._W, mu_f)
            if Re_c < 30.0:
                h_avg = self._heat_trf_coeff_Re_low(
                    Re_c, k_f, mu_f,
                    rho_f, rho_g
                )
            elif 30 <= Re_c <= 1600.0:
                h_avg = self._heat_trf_coeff_Re_medium(
                    Re_c, k_f, mu_f,
                    rho_f, rho_g
                )
            else:  # Re_c > 1600.0
                Pr_f = self._prandtl_number(rho_f, mu_f, k_f, cp_f)
                h_avg = self._heat_trf_coeff_Re_high(
                    Re_c, k_f, mu_f,
                    rho_f, Pr_f, rho_g
                )
            m_dot_new = self._condensate_mass_flow_rate(
                h_avg, h_fg, cp_f,
                self._T_sat, self._T_surf, self._A
            )
            if abs(self._m_dot - m_dot_new) <= self._tol:
                self._m_dot = m_dot_new
                self._V_dot = self._m_dot / rho_f
                self._h_avg = h_avg
                break
            self._m_dot = m_dot_new
            i += 1
        else:
            warnings.warn(
                f'reached maximum number of iterations {self._i_max}'
            )

    @property
    def m_dot(self) -> Quantity:
        """Returns the mass flow rate of condensate being produced."""
        if self._h_avg is None: self._solve_by_iteration()
        return Q_(self._m_dot, 'kg / s')

    # @m_dot.setter
    # def m_dot(self, v: Quantity) -> None:
    #     """Set initial guess of the mass flow rate of condensate being
    #     produced.
    #     """
    #     self._m_dot = v.to('kg / s').m
    #     self._h_avg = None

    @property
    def V_dot(self) -> Quantity:
        """Returns the volume flow rate of condensate being produced."""
        if self._h_avg is None: self._solve_by_iteration()
        return Q_(self._V_dot, 'm ** 3 / s')

    @property
    def h_avg(self) -> Quantity:
        """Returns the average heat transfer coefficient for film
        condensation.
        """
        if self._h_avg is None: self._solve_by_iteration()
        return Q_(self._h_avg, 'W / (m ** 2 * K)')


class HorizontalPlate:
    """Class for calculating the average heat transfer coefficient for film
    condensation of a fluid on a cool horizontal flat surface (plate).

    A distinction must be made between a downward and upward facing plate:
    In the case of a downward facing plate, instantiate from
    `HorizontalPlate.Downward`.
    In the case of an upward facing plate, instantiate from
    `HorizontalPlate.Upward`.

    After instantiation, set the surface temperature of the cool plate and the
    saturation temperature of the condensing fluid through the appropriate
    setters `T_surf` and `T_sat`. To get the average heat transfer coefficient,
    use the `h_avg` getter.
    """
    g: Quantity = Q_(9.81, 'm / s ** 2')

    class Downward:
        """Correlation of Gerstmann and Griffith (1967) for condensation on
        a horizontal, downward facing plate where condensate is removed in
        the form of droplets that form, grow, and detach.
        """
        @staticmethod
        def _rayleigh_number(
            T_surf: float,
            T_sat: float,
            rho_f: float,
            mu_f: float,
            k_f: float,
            cp_f: float,
            sigma: float,
            rho_g: float,
            h_fg: float,
            theta: float
        ) -> float:
            g = HorizontalPlate.g.to('m / s ** 2').m
            g = g * np.cos(theta)
            h_fg_star = h_fg + 0.68 * cp_f * (T_sat - T_surf)
            n1 = g * rho_f * (rho_f - rho_g) * h_fg_star
            d1 = mu_f * (T_sat - T_surf) * k_f
            n2 = sigma
            d2 = (rho_f - rho_g) * g
            Ra = (n1 / d1) * (n2 / d2) ** (3 / 2)
            return Ra

        @staticmethod
        def _heat_trf_coeff_Ra_low(
            Ra: float,
            rho_f: float,
            k_f: float,
            sigma: float,
            rho_g: float,
            theta: float
        ) -> float:
            g = HorizontalPlate.g.to('m / s ** 2').m
            g = g * np.cos(theta)
            n = sigma
            d = (rho_f - rho_g) * g
            z = (n / d) ** (-1 / 2) * k_f
            h_avg = 0.69 * (Ra ** 0.2) * z
            return h_avg

        @staticmethod
        def _heat_trf_coeff_Ra_high(
            Ra: float,
            rho_f: float,
            k_f: float,
            sigma: float,
            rho_g: float,
            theta: float
        ) -> float:
            g = HorizontalPlate.g.to('m / s ** 2').m
            g = g * np.cos(theta)
            n = sigma
            d = (rho_f - rho_g) * g
            z = (n / d) ** (-1 / 2) * k_f
            h_avg = 0.81 * (Ra ** 0.193) * z
            return h_avg

        def __init__(
            self,
            fluid: Fluid,
            theta: Quantity = Q_(0.0, 'deg')
        ) -> None:
            """Creates a `HorizontalPlate.Downward` instance.

            Parameters
            ----------
            fluid: Fluid
                The condensing fluid.
            theta: PlainQuantity
                Inclination angle of the plate. It should be less than 20°. If
                a value greater than 20° is given, a warning will be sent and
                `theta` will be limited to 20°.
            """
            self.fluid = fluid
            if theta.to('deg').m > 20.0:
                warnings.warn('theta too high and limited to 20°')
                theta = Q_(20, 'deg')
            self._theta = theta.to('rad').m
            self._T_surf: float | None = None
            self._T_sat: float | None = None

        @property
        def T_surf(self) -> Quantity:
            """Returns the surface temperature of the cooled plate."""
            return Q_(self._T_surf, 'K')

        @T_surf.setter
        def T_surf(self, v: Quantity) -> None:
            """Set the surface temperature of the cooled plate."""
            self._T_surf = v.to('K').m

        @property
        def T_sat(self) -> Quantity:
            """Returns the (saturation) temperature of the fluid in contact with
            the cooled plate.
            """
            return Q_(self._T_sat, 'K')

        @T_sat.setter
        def T_sat(self, v: Quantity) -> None:
            """Set the (saturation) temperature of the fluid in contact with
            the cooled plate.
            """
            self._T_sat = v.to('K').m

        @property
        def h_avg(self) -> Quantity:
            """Returns the average heat transfer coefficient for film
            condensation.
            """
            T_flm = (self.T_sat + self.T_surf) / 2
            liq_sat = self.fluid(T=T_flm, x=Q_(0.0, 'frac'))
            vap_sat = self.fluid(T=T_flm, x=Q_(1.0, 'frac'))
            mu_f = liq_sat.mu.to('Pa * s').m
            k_f = liq_sat.k.to('W / (m * K)').m
            rho_f = liq_sat.rho.to('kg / m ** 3').m
            cp_f = liq_sat.cp.to('J / (kg * K)').m
            h_f = liq_sat.h.to('J / kg').m
            sigma = liq_sat.sigma('N / m').m
            rho_g = vap_sat.rho.to('kg / m ** 3').m
            h_g = vap_sat.h.to('J / kg').m
            h_fg = h_g - h_f
            Ra = self._rayleigh_number(
                self._T_surf, self._T_sat, rho_f,
                mu_f, k_f, cp_f, sigma, rho_g, h_fg,
                self._theta
            )
            if Ra < 1.e8:
                if Ra < 1.e6:
                    warnings.warn('Ra smaller than lower limit')
                h_avg = self._heat_trf_coeff_Ra_low(
                    Ra, rho_f, k_f,
                    sigma, rho_g, self._theta
                )
                return Q_(h_avg, 'W / (m ** 2 * K)')
            if Ra >= 1.e8:
                if Ra > 1.e10:
                    warnings.warn('Ra greater than upper limit')
                h_avg = self._heat_trf_coeff_Ra_high(
                    Ra, rho_f, k_f,
                    sigma, rho_g, self._theta
                )
                return Q_(h_avg, 'W / (m ** 2 * K)')

    class Upward:
        """Correlation of Nimmo and Leppert (1970) for condensation on
        a horizontal, upward facing plate that is infinite in one direction
        and has length L in the other.
        """
        def __init__(
            self,
            fluid: Fluid,
            L: Quantity
        ) -> None:
            """Creates a `HorizontalPlate.Upward` instance.

            Parameters
            ----------
            fluid: Fluid
                The condensing fluid.
            L: PlainQuantity
                Length of the plate, assuming the width W is infinite (W >> L).
            """
            self.fluid = fluid
            self._L = L.to('m').m
            self._T_surf: float | None = None
            self._T_sat: float | None = None

        @property
        def T_surf(self) -> Quantity:
            """Returns the surface temperature of the cooled plate."""
            return Q_(self._T_surf, 'K')

        @T_surf.setter
        def T_surf(self, v: Quantity) -> None:
            """Set the surface temperature of the cooled plate."""
            self._T_surf = v.to('K').m

        @property
        def T_sat(self) -> Quantity:
            """Returns the (saturation) temperature of the fluid in contact with
            the cooled plate.
            """
            return Q_(self._T_sat, 'K')

        @T_sat.setter
        def T_sat(self, v: Quantity) -> None:
            """Set the (saturation) temperature of the fluid in contact with
            the cooled plate.
            """
            self._T_sat = v.to('K').m

        @property
        def h_avg(self) -> Quantity:
            """Returns the average heat transfer coefficient for film
            condensation.
            """
            T_flm = (self.T_sat + self.T_surf) / 2
            liq_sat = self.fluid(T=T_flm, x=Q_(0.0, 'frac'))
            vap_sat = self.fluid(T=T_flm, x=Q_(1.0, 'frac'))
            mu_f = liq_sat.mu.to('Pa * s').m
            k_f = liq_sat.k.to('W / (m * K)').m
            rho_f = liq_sat.rho.to('kg / m ** 3').m
            cp_f = liq_sat.cp.to('J / (kg * K)').m
            h_f = liq_sat.h.to('J / kg').m
            h_g = vap_sat.h.to('J / kg').m
            h_fg = h_g - h_f
            h_fg_star = h_fg + 0.68 * cp_f * (self._T_sat - self._T_surf)
            g = HorizontalPlate.g.to('m / s ** 2').m
            n = (rho_f ** 2) * g * h_fg_star * (self._L ** 3)
            d = mu_f * (self._T_sat - self._T_surf) * k_f
            h_avg = 0.82 * ((n / d) ** (1 / 5)) * k_f / self._L
            return Q_(h_avg, 'W / (m ** 2 * K)')


class HorizontalCylinder:
    """Class for calculating the average heat transfer coefficient for film
    condensation of a fluid on a cool horizontal cylinder.

    After instantiation, set the outer surface temperature of the cool cylinder
    and the saturation temperature of the condensing fluid through the appropriate
    setters `T_surf` and `T_sat`. To get the average heat transfer coefficient,
    use the `h_avg` getter.

    The average heat transfer coefficient is calculated with the correlation
    of Marto (1998).
    """
    g: Quantity = Q_(9.81, 'm / s ** 2')

    def __init__(
        self,
        fluid: Fluid,
        D: Quantity
    ) -> None:
        """Creates a `HorizontalCylinder` instance.

        Parameters
        ----------
        fluid: Fluid
            The condensing fluid.
        D: PlainQuantity
            Diameter of the cylinder.
        """
        self.fluid = fluid
        self._D = D.to('m').m
        self._T_surf: float | None = None
        self._T_sat: float | None = None

    @property
    def T_surf(self) -> Quantity:
        """Returns the outer surface temperature of the cooled cylinder."""
        return Q_(self._T_surf, 'K')

    @T_surf.setter
    def T_surf(self, v: Quantity) -> None:
        """Set the outer surface temperature of the cooled cylinder."""
        self._T_surf = v.to('K').m

    @property
    def T_sat(self) -> Quantity:
        """Returns the (saturation) temperature of the fluid in contact with
        the outer cooled cylinder surface.
        """
        return Q_(self._T_sat, 'K')

    @T_sat.setter
    def T_sat(self, v: Quantity) -> None:
        """Set the (saturation) temperature of the fluid in contact with
        the outer cooled cylinder surface.
        """
        self._T_sat = v.to('K').m

    @property
    def h_avg(self) -> Quantity:
        """Returns the average heat transfer coefficient for laminar film
        condensation.
        """
        T_flm = (self.T_sat + self.T_surf) / 2
        liq_sat = self.fluid(T=T_flm, x=Q_(0.0, 'frac'))
        vap_sat = self.fluid(T=T_flm, x=Q_(1.0, 'frac'))
        mu_f = liq_sat.mu.to('Pa * s').m
        k_f = liq_sat.k.to('W / (m * K)').m
        rho_f = liq_sat.rho.to('kg / m ** 3').m
        cp_f = liq_sat.cp.to('J / (kg * K)').m
        h_f = liq_sat.h.to('J / kg').m
        rho_g = vap_sat.rho.to('kg / m ** 3').m
        h_g = vap_sat.h.to('J / kg').m
        h_fg = h_g - h_f
        h_fg_star = h_fg + 0.68 * cp_f * (self._T_sat - self._T_surf)
        g = HorizontalCylinder.g.to('m / s ** 2').m
        n = rho_f * (rho_f - rho_g) * g * h_fg_star * (self._D ** 3)
        d = mu_f * (self._T_sat - self._T_surf) * k_f
        h_avg = 0.728 * ((n / d) ** (1 / 4)) * k_f / self._D
        return Q_(h_avg, 'W / (m ** 2 * K)')


class HorizontalCylinderBank(HorizontalCylinder):
    """Class for calculating the average heat transfer coefficient for film
    condensation of a fluid on cool bundles of tubes.

    After instantiation, set the surface temperature of the cool plate and the
    saturation temperature of the condensing fluid through the appropriate
    setters `T_surf` and `T_sat`. To get the average heat transfer coefficient,
    use the `h_avg` getter.

    The average heat transfer coefficient is calculated with the correlation
    of Kern (1958) as presented by Kakaç and Liu (1998).
    """

    def __init__(
        self,
        fluid: Fluid,
        D: Quantity,
        N_vert: int
    ) -> None:
        """Creates a `HorizontalCylinderBank` instance.

        Parameters
        ----------
        fluid: Fluid
            The condensing fluid.
        D: PlainQuantity
            Diameter of the cylinder.
        N_vert:
            The number of rows of tubes in the vertical direction.
        """
        super().__init__(fluid, D)
        self._N_vert = N_vert

    @property
    def h_avg(self) -> Quantity:
        h_avg_single = super().h_avg.magnitude
        h_avg = h_avg_single * (self._N_vert ** (-1 / 6))
        return Q_(h_avg, 'W / (m ** 2 * K)')


class HorizontalFinnedTube:
    g: Quantity = Q_(9.81, 'm / s ** 2')

    @staticmethod
    def _effective_diameter(
        eta_fin: float,
        A_f: float,
        A_eff: float,
        A_uf: float,
        L_tilde: float,
        D_r: float
    ) -> float:
        n1 = 1.30 * eta_fin * A_f
        d1 = A_eff * L_tilde ** (1 / 4)
        n2 = A_uf
        d2 = A_eff * D_r ** (1 / 4)
        D_eff = (n1 / d1 + n2 / d2) ** -4
        return D_eff

    @staticmethod
    def _side_surface_single_fin(D_o: float, D_r: float) -> float:
        A_f = 2 * (np.pi / 4) * (D_o ** 2 - D_r ** 2)
        return A_f

    @staticmethod
    def _exposed_tube_surface(D_o: float, p: float, th: float) -> float:
        A_uf = np.pi * D_o * (p - th)
        return A_uf

    @staticmethod
    def _effective_fin_area(eta_fin: float, A_f: float, A_uf: float) -> float:
        A_eff = eta_fin * A_f + A_uf
        return A_eff

    @staticmethod
    def _dimensionless_length(D_o: float, D_r: float) -> float:
        L_tilde = np.pi * (D_o ** 2 - D_r ** 2) / (4 * D_o)
        return L_tilde

    def __init__(
        self,
        D_r: Quantity,
        D_o: Quantity,
        th: Quantity,
        p: Quantity,
        k: Quantity,
        fluid: Fluid,
        i_max: int = 20,
        tol: Quantity = Q_(0.001, 'frac')
    ) -> None:
        """Creates a `HorizontalFinnedTube` instance.

        Parameters
        ----------
        D_r: PlainQuantity
            Root diameter of the tube.
        D_o: PlainQuantity
            Outer diameter of the finned tube.
        th: PlainQuantity
            Thickness of the fins.
        p: PlainQuantity
            Fin pitch.
        k: PlainQuantity
            Conductivity of fin material.
        fluid: Fluid
            The fluid condensing on the finned tube.
        i_max: int, default 20
            Maximum number of iterations to find the heat transfer coefficient.
        tol: float, default 0.01
            The tolerated deviation between the current and previous value of
            the fin efficiency at which the iteration is stopped.
        """
        self._D_r = D_r.to('m').m
        self._D_o = D_o.to('m').m
        self._th = th.to('m').m
        self._p = p.to('m').m
        self.fluid = fluid
        self._i_max = i_max
        self._tol = tol.to('frac').m
        self._T_surf: float | None = None
        self._T_sat: float | None = None
        self._fin = CircularRectangularFin(D_r / 2, D_o / 2, th, k)

    @property
    def T_surf(self) -> Quantity:
        """Returns the base temperature of the finned tube."""
        return Q_(self._T_surf, 'K')

    @T_surf.setter
    def T_surf(self, v: Quantity) -> None:
        """Set the base temperature of the finned tube."""
        self._T_surf = v.to('K').m
        self._fin.T_b = v

    @property
    def T_sat(self) -> Quantity:
        """Returns the (saturation) temperature of the fluid in contact with
        the outer surface of the finned tube.
        """
        return Q_(self._T_sat, 'K')

    @T_sat.setter
    def T_sat(self, v: Quantity) -> None:
        """Set the (saturation) temperature of the fluid in contact with
        the outer surface of the finned tube.
        """
        self._T_sat = v.to('K').m
        self._fin.T_fluid = v

    @property
    def h_avg(self) -> Quantity:
        """Returns the average heat transfer coefficient for laminar film
        condensation.
        """
        T_flm = (self.T_sat + self.T_surf) / 2
        liq_sat = self.fluid(T=T_flm, x=Q_(0.0, 'frac'))
        vap_sat = self.fluid(T=T_flm, x=Q_(1.0, 'frac'))
        mu_f = liq_sat.mu.to('Pa * s').m
        k_f = liq_sat.k.to('W / (m * K)').m
        rho_f = liq_sat.rho.to('kg / m ** 3').m
        cp_f = liq_sat.cp.to('J / (kg * K)').m
        h_f = liq_sat.h.to('J / kg').m
        h_g = vap_sat.h.to('J / kg').m
        h_fg = h_g - h_f
        h_fg_star = h_fg + 0.68 * cp_f * (self._T_sat - self._T_surf)
        g = HorizontalFinnedTube.g.to('m / s ** 2').m
        L_tilde = self._dimensionless_length(self._D_o, self._D_r)
        A_f = self._side_surface_single_fin(self._D_o, self._D_r)
        A_uf = self._exposed_tube_surface(self._D_o, self._p, self._th)
        eta_fin = 0.85  # initial guess
        h_avg = float('nan')
        i = 0
        while i <= self._i_max:
            A_eff = self._effective_fin_area(eta_fin, A_f, A_uf)
            D_eff = self._effective_diameter(eta_fin, A_f, A_eff, A_uf, L_tilde, self._D_r)
            n = (rho_f ** 2) * (k_f ** 3) * g * h_fg_star
            d = mu_f * (self._T_sat - self._T_surf) * D_eff
            h_avg = 0.689 * (n / d) ** (1 / 4)
            self._fin.h_avg = Q_(h_avg, 'W / (m ** 2 * K)')
            eta_fin_new = self._fin.efficiency
            if abs(eta_fin - eta_fin_new) < self._tol:
                break
            eta_fin = eta_fin_new
            i += 1
        else:
            warnings.warn(
                f'reached maximum number of iterations {self._i_max}'
            )
        return Q_(h_avg, 'W / (m ** 2 * K)')
