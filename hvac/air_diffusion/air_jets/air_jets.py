"""
HORIZONTAL AND INCLINED, ISOTHERMAL OR NON-ISOTHERMAL FREE AIR JETS

References
----------
[1] Hagström K., Sirén K., Zhivov A. M. (1999). Calculation Methods for Air
    Supply Design in Industrial Facilities. Helsinki University of Technology,
    Laboratory of Heating, Ventilating and Air Conditioning.
[2] Goodfellow, H. D., & Kosonen, R. (2020). Industrial Ventilation Design
    Guidebook: Volume 1: Fundamentals. Academic Press.
[3] Awbi, H. B. (2003). Ventilation of Buildings. Taylor & Francis.
"""
import warnings
from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import root_scalar
from scipy.special import erf
from hvac import Quantity
from hvac.fluids import FluidState
from .archimedes_number import archimedes_number

Q_ = Quantity


class Jet(ABC):

    def __init__(self):
        self.x_o = Q_(0.0, 'm')
        self.room_air = None

    def _find_supply_outlet_distance(self, U_o: Quantity) -> Quantity:
        """Finds the horizontal distance `x` for which the jet centerline
        velocity `U_o`, calculated with `centerline_velocity_zone3()`,
        equals the supply air outlet velocity `U_o`.
        Assuming that the jet air velocity in zone 1 (potential core region)
        and zone 2 (characteristic decay region) remains constant (which in fact
        is not true for zone 2), this function returns the approximate
        horizontal position of the supply outlet along the x-axis with respect
        to the origin (x = 0) in the "zone 3 equations" (see the functions that
        end with "_zone3" in their name).
        """
        def _fun(x: float) -> float:
            x = Q_(x, 'm')
            U_x = self.centerline_velocity_zone3(x)
            dev = (U_x - U_o).magnitude
            return dev

        x_min = Q_(0.01, 'm')
        x_max = self.throw_zone3(U_x=Q_(0.05, 'm / s')).to('m')
        sol = root_scalar(_fun, bracket=(x_min.m, x_max.m))
        x_o = Q_(sol.root, 'm')
        return x_o

    def _temperature_correction(self, x: Quantity) -> Quantity:
        """Calculates a correction of the dimensionless axial distance
        `sqrt(A_o)/x` as `sqrt(rho_x / rho_r)` = `sqrt(T_r / T_x)` due to
        mass density difference between supply air and return air (i.e.
        buoyancy). For isothermal jet flow, the temperature correction will be
        1.
        """
        T_x = self.centerline_temperature_zone3(x)
        tc = np.sqrt(self.room_air.T / T_x)
        return tc

    @abstractmethod
    def centerline_velocity_univ(self, x: Quantity) -> Quantity:
        """Returns the centerline velocity of the jet at a distance `x` from the
        supply opening. The equation used is "universal", as this equation is
        valid for zone 1, zone 2, and zone 3 of the air jet.
        """
        ...

    @abstractmethod
    def centerline_velocity_zone3(self, x: Quantity) -> Quantity:
        """Returns the centerline velocity of the jet at a distance `x` from the
        supply opening. Results are only valid in the fully developed turbulent
        region of the jet (zone 3).
        """
        ...

    @abstractmethod
    def velocity_profile_zone3(self, x: Quantity, y: Quantity) -> Quantity:
        """Returns the velocity in the cross-section of the jet at a distance
        `x` from the supply opening and a vertical/radial distance `y` from
        the centerline of the jet. Results are only valid in the fully developed
        turbulent region of the jet (zone 3).
        """
        ...

    def throw_zone3(self, U_x: Quantity) -> Quantity:
        """Returns the distance from the supply opening where the jet's
        centerline velocity equals the value `U_x`. Results are only valid for
        the fully developed turbulent region of an isothermal jet (zone 3).
        """
        U_x_target = U_x.to('m / s').magnitude

        def _fun(x: float) -> tuple[float, float]:
            U_x = self.centerline_velocity_zone3(Q_(x, 'm')).magnitude
            dev = U_x - U_x_target
            der = -U_x / x
            return dev, der

        sol = root_scalar(_fun, method='newton', fprime=True, x0=1.0)
        if sol.converged:
            throw = Q_(sol.root, 'm')
            return throw
        else:
            warnings.warn(
                f"throw length to {U_x:~P.1f} cannot be attained",
                category=RuntimeWarning
            )
            return Q_(float('nan'), 'm')

    @abstractmethod
    def entrainment_ratio_zone3(self, x: Quantity) -> Quantity:
        """Returns the fractional volume flow rate through a cross-section of
        the jet at a distance `x` from the supply opening with respect to the
        volume flow rate at the supply outlet. Results are only valid in the
        fully developed turbulent region of an isothermal jet (zone 3).
        """
        ...

    @abstractmethod
    def centerline_temperature_univ(self, x: Quantity) -> Quantity:
        """Returns the centerline temperature of the jet at a horizontal
        distance `x` from the supply opening. The equation used is "universal",
        as this equation is valid for zone 1, zone 2, and zone 3 of the air jet.
        """
        ...

    @abstractmethod
    def centerline_temperature_zone3(self, x: Quantity) -> Quantity:
        """Returns the centerline temperature of the jet at a distance `x` from
        the supply opening. Results are only valid in the fully developed
        turbulent region of the jet (zone 3).
        """
        ...

    @abstractmethod
    def temperature_profile_zone3(self, x: Quantity, y: Quantity) -> Quantity:
        """Returns the temperature in the cross-section of the jet at a distance
        `x` from the supply opening and a vertical or radial distance `y` from
        the centerline of the jet. Results are only valid in the fully developed
        turbulent region of the jet (zone 3).
        """
        ...

    @abstractmethod
    def trajectory(self, x: Quantity, alpha: Quantity = Q_(0, 'deg')) -> Quantity:
        """Returns the vertical distance `y` of the centerline of the jet at a
        horizontal distance `x` from the supply opening. Parameter `alpha` is
        the projection angle of the jet, measured with respect to a horizontal
        plane that goes through the center of the supply opening. The results
        are valid for projection angles between -45° (downward projection) and
        +45° (upward projection).
        """
        ...

    @abstractmethod
    def extremum(self, alpha: Quantity) -> tuple[Quantity, Quantity] | None:
        """Returns the coordinates of the highest or lowest point along the
        jet's trajectory. Parameter `alpha` is the projection angle of the
        jet. The function returns a result in case of a hot jet being projected
        downwards or a cold jet being projected upwards, otherwise `None` is
        returned.
        """
        ...

    @abstractmethod
    def local_archimedes_number_zone3(self, x: Quantity) -> float:
        """Returns the local Archimedes number of the jet at a distance `x`
        from the supply opening. As long as the local Archimedes number is
        smaller than 0.1 for a compact/radial jet or 0.15 for a linear jet, the
        jet can be considered as being unaffected by buoyancy forces (i.e.
        similar to an isothermal jet).
        """
        ...

    @abstractmethod
    def maximum_supply_air_temperature(self, B: Quantity, H: Quantity) -> Quantity:
        """Returns the supply air temperature for which the jet Archimedes
        number `Ar_x` at a distance `x = 0.22 * K1 * sqrt(BH)` is equal to 0.2.
        According to experimental studies by Grimitlyn (1978) this condition
        assures that a confined air jet still behaves like an isothermal jet
        in rooms with a height/width ratio between 0.3 and 1.0.
        """
        ...

    @abstractmethod
    def linear_jet_length(self) -> tuple[Quantity]:
        """Returns the length of the linear jet zone where buoyancy forces are
        negligibly small.
        """
        ...


class CompactJet(Jet):
    """Class that represents jets formed by cylindrical tubes, nozzles,
    square or rectangular openings with a small aspect ratio (< 40), unshaded or
    shaded by perforated plates, grilles, etc.
    """
    Pr = 0.75  # overall turbulent Prandtl number
    psi = 0.47

    def __init__(
        self,
        A_o: Quantity,
        U_o: Quantity,
        supply_air: FluidState,
        room_air: FluidState,
        K1: float = 7.0
    ) -> None:
        """Creates a `CompactJet` object.

        Parameters
        ----------
        A_o:
            Effective area of the diffuser opening.
        U_o:
            Average supply air velocity at the diffuser outlet.
        supply_air:
            State of the supply air.
        room_air:
            Mean state of the room air (return air).
        K1:
            Dynamic characteristic (throw constant) of the diffuser jet.
        """
        super().__init__()
        self.A_o = A_o.to('m**2')
        self.U_o = U_o.to('m / s')
        self.supply_air = supply_air
        self.room_air = room_air
        self.K1 = K1
        self.c = self._calculate_c()
        # Thermal characteristic of the diffuser jet:
        self.K2 = self.calculate_K2(K1)
        # Archimedes number at the supply outlet:
        d_o = np.sqrt(4 * self.A_o / np.pi).to('m')
        dT_o = (self.supply_air.T - self.room_air.T).to('K')
        self.Ar_o = archimedes_number(d_o, self.U_o, dT_o, self.room_air)
        # Offset distance from supply outlet:
        self.x_o = self._find_supply_outlet_distance(self.U_o)

    @classmethod
    def calculate_K2(cls, K1: float) -> float:
        K2 = np.sqrt((1 + cls.Pr) / 2) * K1
        return K2

    def _calculate_c(self) -> float:
        """The dynamic characteristic `K1` is defined as:
        ```
        K1 = theta * fi / sqrt(PI) * c
        ```
        Given `K1`, this function returns the constant `c` assuming `fi` = 1,
        i.e. assuming that the jet velocity distribution at the outlet is
        uniform.
        """
        theta = 1.0  # = sqrt (rho_o / rho_r), but K1 is specified for isothermal jet flow
        fi = 1.0  # uniform outlet velocity distribution
        c = theta * fi / (np.sqrt(np.pi) * self.K1)
        return c

    def centerline_velocity_univ(self, x: Quantity, AR: float = 1.0) -> Quantity:
        x = x.to('m').magnitude
        A_o = self.A_o.to('m**2').magnitude
        h_o = np.sqrt(A_o / AR)
        b_o = A_o / h_o
        k = erf(h_o / (2 * self.c * x)) * erf(b_o / (2 * self.c * x))
        U_x = np.sqrt(k) * self.U_o.to('m / s')
        return U_x

    def centerline_velocity_zone3(self, x: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        tc = self._temperature_correction(x)
        U_x = self.K1 * self.U_o * tc * np.sqrt(self.A_o) / x
        return U_x.to('m / s')

    def velocity_profile_zone3(self, x: Quantity, y: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        y = y.to('m')
        alpha = Q_(5.1, 'deg').to('rad')
        eta = y / (x * np.tan(alpha))
        k = np.exp(-0.7 * eta ** 2)
        U_x = self.centerline_velocity_zone3(x)
        U_y = k * U_x
        return U_y

    def entrainment_ratio_zone3(self, x: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        tc = self._temperature_correction(x)
        V_dot_frac = (2 / self.K1) * (1 / tc) * (x / np.sqrt(self.A_o))
        return V_dot_frac

    def centerline_temperature_univ(self, x: Quantity, AR: float = 1.0) -> Quantity:
        x = x.to('m').magnitude
        A_o = self.A_o.to('m**2').magnitude
        h_o = np.sqrt(A_o / AR)
        b_o = A_o / h_o
        p = np.sqrt((1 + self.Pr) / 2)
        b = b_o / (2 * self.c * x)
        h = h_o / (2 * self.c * x)
        n = erf(p * b) * erf(p * h)
        d = np.sqrt(erf(b) * erf(h))
        theta_x = n / d
        T_r = self.room_air.T.to('K')
        T_o = self.supply_air.T.to('K')
        T_x = T_r + theta_x * (T_o - T_r)
        return T_x

    def centerline_temperature_zone3(self, x: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        theta_x = self.K2 * np.sqrt(self.A_o) / x
        T_r = self.room_air.T.to('K')
        T_o = self.supply_air.T.to('K')
        T_x = T_r + theta_x * (T_o - T_r)
        return np.maximum(T_x, Q_(1.e-12, 'K'))

    def temperature_profile_zone3(self, x: Quantity, y: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        y = y.to('m')
        alpha_U05 = Q_(2 * 5.1, 'deg').to('rad')  # see ref. [1], table 4.2
        eta = y / (x * np.tan(alpha_U05) / np.sqrt(self.Pr))
        k = np.exp(-0.7 * eta ** 2)
        T_x = self.centerline_temperature_zone3(x)
        T_y = k * T_x
        return T_y

    def _sign(self, alpha: Quantity) -> int:
        if self.supply_air.T > self.room_air.T:
            # hot jet
            if alpha.m < 0:
                return -1
            else:
                return 1
        elif self.supply_air.T < self.room_air.T:
            # cold jet
            if alpha.m < 0:
                return 1
            else:
                return -1
        else:
            # isothermal jet
            return 0

    def trajectory(self, x: Quantity, alpha: Quantity = Q_(0, 'deg')) -> Quantity:
        x = x.to('m')
        alpha = alpha.to('rad')
        a = x * np.tan(alpha) / np.sqrt(self.A_o)
        b = (
            self.psi * (self.K2 / self.K1 ** 2) * self.Ar_o
            * (x / (np.cos(alpha) * np.sqrt(self.A_o))) ** 3
        )
        k = a + b
        y = k * np.sqrt(self.A_o)
        return y.to('m')

    def extremum(self, alpha: Quantity) -> tuple[Quantity, Quantity] | None:
        if self._sign(alpha) < 0:
            alpha = alpha.to('rad')
            n = np.cos(alpha) * np.sqrt(abs(np.sin(alpha))) * np.sqrt(self.A_o)
            d = np.sqrt(3 * self.psi * (self.K2 / self.K1 ** 2) * abs(self.Ar_o))
            x_v = n / d
            y_v = self.trajectory(x_v, alpha)
            return x_v, y_v
        return None

    def local_archimedes_number_zone3(self, x: Quantity) -> float:
        x = (self.x_o + x).to('m')
        # tc = self._temperature_correction(x)
        tc = 1.0
        Ar_x = (
            self.K2 / self.K1 ** 2
            * self.Ar_o * (1 / tc) * (x / np.sqrt(self.A_o)) ** 2
        )
        return Ar_x

    def maximum_supply_air_temperature(self, B: Quantity, H: Quantity) -> Quantity:
        B = B.to('m')
        H = H.to('m')
        g = Q_(9.81, 'm / s**2')
        Ar_x_max = 0.2
        T_r = self.room_air.T.to('K')
        k = Ar_x_max * T_r / (0.22 ** 2 * g)
        n = k * self.U_o ** 2 * np.sqrt(self.A_o)
        d = self.K2 * B * H
        dT_o = n / d
        if self.supply_air.T.to('K') < T_r:
            # cooling is expected
            dT_o = -dT_o
        T_o = T_r + dT_o
        return T_o

    def linear_jet_length(self) -> tuple[Quantity]:
        k = 0.5 * (
            (1 / abs(self.Ar_o)) ** (2 / 3)
            * (self.room_air.T / self.supply_air.T) ** (-1 / 3)
        )
        x_1 = k * np.sqrt(self.A_o)
        return x_1


class LinearJet(Jet):
    """Class that represents jets formed by slots or rectangular openings with
    a large aspect ratio (> 40).
    """
    Pr = 0.5  # overall turbulent Prandtl number
    psi = 0.4

    def __init__(
        self,
        h_o: Quantity,
        U_o: Quantity,
        supply_air: FluidState,
        room_air: FluidState,
        K1: float = 2.43,
        K2: float = 2.00
    ) -> None:
        """Creates a `LinearJet` object.

        Parameters
        ----------
        h_o:
            Height of the slot.
        U_o:
            Average air velocity at diffuser outlet.
        supply_air:
            State of supply air.
        room_air:
            Mean state of room air (return air).
        K1:
            Dynamic characteristic (throw constant) of the diffuser jet.
        K2:
            Thermal characteristic of the diffuser jet.
        """
        super().__init__()
        self.h_o = h_o.to('m')
        self.U_o = U_o.to('m / s')
        self.supply_air = supply_air
        self.room_air = room_air
        self.K1 = K1
        self.c = self._calculate_c()
        self.K2 = K2
        dT_o = (self.supply_air.T - self.room_air.T).to('K')
        self.Ar_o = archimedes_number(self.h_o, self.U_o, dT_o, self.room_air)
        # Offset distance from supply outlet:
        self.x_o = self._find_supply_outlet_distance(self.U_o)

    def _calculate_c(self) -> float:
        """The dynamic characteristic `K1` is calculated as:
        ```
        K1 = theta * fi / sqrt(PI) * c
        ```
        Given `K1`, this function returns the constant `c` assuming `fi` = 1,
        i.e. assuming that the jet velocity distribution at the outlet is
        uniform.
        """
        theta = 1.0  # np.sqrt(self.room_air.T.m / self.supply_air.T.m)
        fi = 1.0  # uniform outlet velocity distribution
        c = theta * fi / (np.sqrt(np.pi) * self.K1)
        return c

    def centerline_velocity_univ(self, x: Quantity) -> Quantity:
        x = x.to('m').magnitude
        h_o = self.h_o.to('m').magnitude
        k = erf(h_o / (2 * self.c * x))
        U_x = np.sqrt(k) * self.U_o.to('m / s')
        return U_x

    def centerline_velocity_zone3(self, x: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        tc = self._temperature_correction(x)
        U_x = self.K1 * self.U_o * tc * np.sqrt(self.h_o / x)
        return U_x.to('m / s')

    def velocity_profile_zone3(self, x: Quantity, y: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        y = y.to('m')
        alpha = Q_(5.95, 'deg').to('rad')
        eta = y / (x * np.tan(alpha))
        k = np.exp(-0.7 * eta ** 2)
        U_x = self.centerline_velocity_zone3(x)
        U_y = k * U_x
        return U_y

    def entrainment_ratio_zone3(self, x: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        tc = self._temperature_correction(x)
        V_dot_frac = np.sqrt(2 * (1 / tc) * x) / (self.K1 * np.sqrt(self.h_o))
        return V_dot_frac

    def centerline_temperature_univ(self, x: Quantity) -> Quantity:
        """Returns the centerline temperature of the jet at a horizontal
        distance `x` from the supply opening, assuming a slot with an infinite
        length. The equation used is "universal", as this equation is valid for
        zone 1, zone 2, and zone 3 of the air jet.
        """
        x = x.to('m').magnitude
        h_o = self.h_o.to('m').magnitude
        p = np.sqrt((1 + self.Pr) / 2)
        h = h_o / (2 * self.c * x)
        n = erf(p * h)
        d = np.sqrt(erf(h))
        theta_x = n / d
        T_r = self.room_air.T.to('K')
        T_o = self.supply_air.T.to('K')
        T_x = T_r + theta_x * (T_o - T_r)
        return T_x

    def centerline_temperature_zone3(self, x: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        theta_x = self.K2 * np.sqrt(self.h_o / x)
        T_r = self.room_air.T.to('K')
        T_o = self.supply_air.T.to('K')
        T_x = T_r + theta_x * (T_o - T_r)
        return T_x

    def temperature_profile_zone3(self, x: Quantity, y: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        y = y.to('m')
        alpha_U05 = Q_(2 * 5.9, 'deg').to('rad')  # see ref. [1], table 4.2
        eta = y / (x * np.tan(alpha_U05) / np.sqrt(self.Pr))
        k = np.exp(-0.7 * eta ** 2)
        T_x = self.centerline_temperature_zone3(x)
        T_y = k * T_x
        return T_y

    def _sign(self, alpha: Quantity) -> int:
        if self.supply_air.T > self.room_air.T:
            # hot jet
            if alpha.m < 0:
                return -1
            else:
                return 1
        elif self.supply_air.T < self.room_air.T:
            # cold jet
            if alpha.m < 0:
                return 1
            else:
                return -1
        else:
            # isothermal jet
            return 0

    def trajectory(self, x: Quantity, alpha: Quantity = Q_(0, 'deg')) -> Quantity:
        alpha = alpha.to('rad')
        a = x * np.tan(alpha) / self.h_o
        b = (
            self.psi * (self.K2 / self.K1 ** 2) * self.Ar_o
            * (x / self.h_o) ** (5 / 2)
        )
        k = a + b
        y = k * self.h_o
        return y.to('m')

    def extremum(self, alpha: Quantity) -> tuple[Quantity, Quantity] | None:
        if self._sign(alpha) < 0:
            alpha_ = self._sign(alpha) * alpha.to('rad')
            n = np.tan(alpha_)
            d = (
                self.psi * (self.K2 / self.K1 ** 2)
                * self.Ar_o * (5 / 2) * self.h_o ** (-3 / 2)
            )
            x_v = (n / d) ** (2 / 3)
            y_v = self.trajectory(x_v, alpha)
            return x_v, y_v
        return None

    def local_archimedes_number_zone3(self, x: Quantity) -> float:
        x = (self.x_o + x).to('m')
        Ar_x = (
            self.K2 / self.K1 ** 2
            * self.Ar_o * (x / self.h_o) ** (3 / 2)
        )
        return Ar_x

    def maximum_supply_air_temperature(self, B: Quantity, H: Quantity) -> Quantity:
        B = B.to('m')
        H = H.to('m')
        g = Q_(9.81, 'm / s**2')
        Ar_x_max = 0.2
        T_r = self.room_air.T.to('K')
        k = Ar_x_max * T_r / (0.22 ** 1.5 * g)
        n = np.sqrt(self.K1) * self.U_o ** 2 * np.sqrt(self.h_o)
        d = self.K2 * (B * H) ** 0.75
        dT_o = k * (n / d)
        if self.supply_air.T.to('K') < T_r:
            # cooling is expected
            dT_o = -dT_o
        T_o = T_r + dT_o
        return T_o

    def linear_jet_length(self) -> tuple[Quantity]:
        k = 0.5 * (
            (1 / abs(self.Ar_o)) ** (2 / 3)
            * (self.room_air.T / self.supply_air.T) ** (-1 / 3)
        )
        x_1 = k * self.h_o
        return x_1


class RadialJet(Jet):
    """Class that represents radial and also incomplete radial jets. Radial jets
    are formed by ceiling circular diffusers with flat discs or multidiffusers
    that direct air horizontally in all directions. Incomplete radial jets
    are supplied through outlets with grills having diverging vanes.
    """
    Pr = 0.5  # overall turbulent Prandtl number

    def __init__(
        self,
        r_o: Quantity,
        h_o: Quantity,
        U_o: Quantity,
        supply_air: FluidState,
        room_air: FluidState,
        alpha: Quantity = Q_(0, 'deg'),
        K1: float = 1.05
    ) -> None:
        """Creates a `RadialJet` object.

        Parameters
        ----------
        r_o:
            Effective radius of the diffuser.
        h_o:
            Effective slot height of the diffuser.
        U_o:
            Average air velocity at diffuser outlet.
        supply_air:
            State of supply air.
        room_air:
            Mean state of room air (return air).
        alpha:
            Diffuser outlet angle. This angle is zero for radial jet diffusers
            (i.e. the default).
        K1:
            Dynamic characteristic (throw constant) of the diffuser jet.
        """
        super().__init__()
        self.r_o = r_o.to('m')
        self.h_o = h_o.to('m')
        self.alpha = alpha.to('rad')
        self.U_o = U_o.to('m / s')
        self.supply_air = supply_air
        self.room_air = room_air
        self.K1 = K1
        self.c = self._calculate_c()
        # Thermal characteristic of the diffuser jet:
        self.K2 = self._calculate_K2()
        # Archimedes number at the supply outlet:
        self.A_o = 2 * np.pi * self.r_o * self.h_o
        dT_o = (self.supply_air.T - self.room_air.T).to('K')
        self.Ar_o = archimedes_number(np.sqrt(self.A_o), self.U_o, dT_o, self.room_air)
        # Offset distance from supply outlet:
        self.x_o = self._find_supply_outlet_distance(self.U_o)
        # A number of methods below are identical to the methods of `CompactJet`:
        self._compact_jet = CompactJet(
            A_o=self.A_o,
            U_o=self.U_o,
            supply_air=self.supply_air,
            room_air=self.room_air,
            K1=self.K1
        )
        self._compact_jet.K2 = self.K2
        self._compact_jet.x_o = self.x_o

    def _calculate_c(self) -> float:
        """Derives the diffuser constant `c` from the dynamic characteristic
        `K1`.
        """
        h_o = self.h_o.magnitude
        r_o = self.r_o.magnitude
        alpha = self.alpha.magnitude
        a = (h_o / r_o * np.cos(alpha)) ** 2
        b = np.sqrt(a)
        c = -(self.K1 ** 2)
        D = b ** 2 - 4 * a * c
        if D < 0:
            raise ValueError('no real solution for constant `c`')
        else:
            x1 = -b + np.sqrt(D) / (2 * a)
            x2 = -b - np.sqrt(D) / (2 * a)
            if x1 > 0:
                return x1
            elif x2 > 0:
                return x2
            else:
                raise ValueError('negative solution for constant `c`')

    def _calculate_K2(self) -> float:
        """Calculates the thermal characteristic `K2` of the radial diffuser."""
        beta = 2.0 * np.pi  # jet spreads air in all directions (360°)
        theta = np.sqrt(self.room_air.T.m / self.supply_air.T.m)
        fi = 1.0  # uniform outlet velocity distribution
        alpha = self.alpha.magnitude
        a = np.sqrt(
            (1 + self.Pr)
            / (4 * np.sqrt(np.pi) * self.c * np.sin(alpha))
        )
        b = theta / (fi * np.sqrt(beta))
        K2 = a * b
        return K2

    def centerline_velocity_univ(self, x: Quantity) -> Quantity:
        # uses the same calculation routine as `CompactJet`
        U_x = self._compact_jet.centerline_temperature_univ(x)
        return U_x

    def centerline_velocity_zone3(self, x: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        a = self.c * (self.h_o / self.r_o) * np.cos(self.alpha)
        n = np.sqrt(a * (a + 1))
        d = np.sqrt(x * (x - self.r_o)) / self.r_o
        U_x = (n / d) * self.U_o
        return U_x.to('m / s')

    def velocity_profile_zone3(self, x: Quantity, y: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        y = y.to('m')
        alpha = Q_(10.5, 'deg').to('rad')
        eta = y / (x * np.tan(alpha))
        k = np.exp(-0.7 * eta ** 2)
        U_x = self.centerline_velocity_zone3(x)
        U_y = k * U_x
        return U_y

    def entrainment_ratio_zone3(self, x: Quantity) -> float:
        x = (self.x_o + x).to('m')
        tc = self._temperature_correction(x)
        V_dot_frac = 0.069 * (1 / tc) * x / np.sqrt(self.r_o * self.h_o)
        return V_dot_frac

    def centerline_temperature_univ(self, x: Quantity) -> Quantity:
        # uses the same calculation routine as `CompactJet`
        T_x = self._compact_jet.centerline_temperature_univ(x)
        return T_x

    def centerline_temperature_zone3(self, x: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        theta_x = self.K2 * self.r_o / x
        T_r = self.room_air.T.to('K')
        T_o = self.supply_air.T.to('K')
        T_x = T_r + theta_x * (T_o - T_r)
        return T_x

    def temperature_profile_zone3(self, x: Quantity, y: Quantity) -> Quantity:
        x = (self.x_o + x).to('m')
        y = y.to('m')
        alpha_05 = Q_(2 * 10.5, 'deg').to('rad')  # see ref. [1], table 4.2
        eta = y / (x * np.tan(alpha_05) / np.sqrt(self.Pr))
        k = np.exp(-0.7 * eta ** 2)
        T_x = self.centerline_temperature_zone3(x)
        T_y = k * T_x
        return T_y

    def trajectory(self, x: Quantity, alpha: Quantity = Q_(0, 'deg')) -> Quantity:
        # uses the same calculation routine as `CompactJet`
        y = self._compact_jet.trajectory(x, alpha)
        return y

    def extremum(self, alpha: Quantity) -> tuple[Quantity, Quantity] | None:
        # uses the same calculation routine as `CompactJet`
        x_v, y_v = self._compact_jet.extremum(alpha)
        return x_v, y_v

    def local_archimedes_number_zone3(self, x: Quantity) -> float:
        # uses the same calculation routine as `CompactJet`
        Ar_x = self._compact_jet.local_archimedes_number_zone3(x)
        return Ar_x

    def maximum_supply_air_temperature(self, B: Quantity, H: Quantity) -> Quantity:
        # uses the same calculation routine as `CompactJet`
        T_o = self._compact_jet.maximum_supply_air_temperature(B, H)
        return T_o

    def linear_jet_length(self) -> tuple[Quantity]:
        # uses the same calculation routine as `CompactJet`
        x_1 = self._compact_jet.linear_jet_length()
        return x_1
