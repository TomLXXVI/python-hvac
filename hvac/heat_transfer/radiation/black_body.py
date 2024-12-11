import numpy as np
from scipy.integrate import quad
from hvac import Quantity


Q_ = Quantity


class BlackBody:
    """A blackbody:
    -   absorbs all the radiation that is incident upon it (regardless of the
        direction or wavelength of the incident radiation).
    -   emits radiation uniformly in all directions (i.e., it is a diffuse
        emitter).
    -   emits the maximum possible amount of radiation at a given temperature
        and wavelength.
    """
    C1 = Q_(3.742e8, 'W * µm**4 / m**2')
    C2 = Q_(14_388, 'µm * K')
    sigma = Q_(5.67e-8, 'W / (m**2 * K**4)')  # Stefan-Boltzmann constant

    def __init__(self, T: Quantity, A: Quantity = Q_(1.0, 'm**2')):
        """Creates a `BlackBody` instance.

        Parameters
        ----------
        T:
            Temperature of the blackbody.
        A:
            Surface area of the blackbody.
        """
        self.T = T.to('K')
        self.A = A.to('m**2')

    def spectral_emissive_flux(self, wavelength: Quantity) -> Quantity:
        """Returns the blackbody spectral emissive power per unit area at the
        given wavelength according to Planck's law.
        """
        wl = wavelength.to('µm')
        n = self.C1
        d = wl**5 * (np.exp(self.C2 / (wl * self.T)) - 1)
        q_dot_s = n / d
        return q_dot_s

    def total_emissive_flux(self) -> Quantity:
        """Returns the total blackbody emissive power per unit area."""
        q_dot = self.sigma * self.T ** 4
        return q_dot

    def peak_spectral_emissive_flux(self) -> tuple[Quantity, Quantity]:
        """Returns the peak emissive power per unit area in the blackbody
        spectrum and the wavelength at which the peak is situated according to
        Wien's law.
        """
        wl_max = Q_(2897.8, 'µm * K') / self.T
        q_dot_s_max = self.spectral_emissive_flux(wl_max)
        return q_dot_s_max, wl_max

    def energy_fraction(self, wl1: Quantity, wl2: Quantity):
        """Returns the fraction of the total blackbody emission that occurs
        between wavelength `wl1` and wavelength `wl2`.
        """
        C1 = self.C1.m
        sigma = self.sigma.m
        C2 = self.C2.m
        
        def __fun__(x: float) -> float:
            n = C1
            d = sigma * x**5 * (np.exp(C2 / x) - 1)
            r = n / d
            return r

        wl1 = wl1.to('µm').m
        wl2 = wl2.to('µm').m
        T = self.T.to('K').m
        F = quad(__fun__, wl1 * T, wl2 * T)
        return F[0]

    def total_emissive_power(self) -> Quantity:
        """Returns the total blackbody emissive power emitted from the surface."""
        q_dot = self.total_emissive_flux()
        Q_dot = self.A * q_dot
        return Q_dot


def net_radiation_exchange(
    bb1: BlackBody,
    bb2: BlackBody,
    vf12: float | None = None,
    vf21: float | None = None
) -> Quantity:
    """Returns the net rate of radiation exchange from blackbody surface 1 to
    blackbody surface 2.

    Parameters
    ----------
    bb1:
        Blackbody surface 1
    bb2:
        Blackbody surface 2
    vf12:
        View factor from surface 1 to surface 2, i.e. the ratio of radiation
        leaving surface 1 that goes directly to surface 2 to the total radiation
        leaving surface 1.
    vf21:
        View factor from surface 2 to surface 1.

    Notes
    -----
    Either `vf12` or `vf21` should be specified. If both are specified, `vf12`
    takes precedence.
    """
    q_dot1 = bb1.total_emissive_flux()
    q_dot2 = bb2.total_emissive_flux()
    if vf12 is not None:
        Q_dot12 = bb1.A * vf12 * (q_dot1 - q_dot2)
        return Q_dot12
    elif vf21 is not None:
        Q_dot12 = bb2.A * vf21 * (q_dot1 - q_dot2)
        return Q_dot12


def space_resistance(
    bb1: BlackBody,
    bb2: BlackBody,
    vf12: float | None = None,
    vf21: float | None = None
) -> Quantity:
    """Returns the resistance to radiation heat transfer from blackbody surface
    1 to blackbody surface 2.

    Parameters
    ----------
    bb1:
        Blackbody surface 1
    bb2:
        Blackbody surface 2
    vf12:
        View factor from surface 1 to surface 2, i.e. the ratio of radiation
        leaving surface 1 that goes directly to surface 2 to the total radiation
        leaving surface 1.
    vf21:
        View factor from surface 2 to surface 1.

    Notes
    -----
    Either `vf12` or `vf21` should be specified. If both are specified, `vf12`
    takes precedence.
    """
    if vf12 is not None:
        R12 = 1 / (bb1.A * vf12)
        return R12
    elif vf21 is not None:
        R12 = 1 / (bb2.A * vf21)
        return R12
