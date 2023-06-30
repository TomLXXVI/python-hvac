"""Heat exchanger surface geometrical characteristics and dimensions for
tube-fin heat exchangers.

References
----------
[1]     Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger
        Design. John Wiley & Sons. (§ 8.2)
"""
from typing import TypeVar
import numpy as np
from abc import ABC, abstractmethod
from hvac import Quantity

Q_ = Quantity


def hydraulic_diameter(
    A_o: Quantity | None = None,
    L: Quantity | None = None,
    A: Quantity | None = None,
    sigma: float | None = None,
    alpha: Quantity | None = None
) -> Quantity:
    """Generalized relation for hydraulic diameter.

    Parameters
    ----------
    A_o:
        minimum free flow area [m ** 2]
    L:
        flow length [m]
    A:
        total heat transfer area [m ** 2]
    sigma:
        ratio of free flow area to frontal area [-]
    alpha:
        ratio of total heat transfer area to total volume of heat exchanger [m]

    Notes
    -----
    Two groups of input parameters are possible:
    1. A_o, L, and A
    2. sigma and alpha

    Returns
    -------
    hydraulic diameter [m]
    """
    if all([A_o, L, A]):
        D_h = 4 * A_o * L / A
    elif all([sigma, alpha]):
        D_h = 4 * sigma / alpha
    else:
        D_h = None
    return D_h.to('m')


def total_number_of_tubes(
    L2: Quantity,
    L3: Quantity,
    S_t: Quantity,
    S_l: Quantity,
    arrangement: str = 'staggered'
) -> float:
    """Total number of tubes in the heat exchanger.

    Parameters
    ----------
    L2:
        Length of the tube bank parallel to the direction of external flow [m].
    L3:
        Length of the tube bank perpendicular to the direction of external flow
        [m].
    S_t:
        Lateral or transverse pitch, i.e. spacing between tubes of the same
        row [m].
    S_l:
        Longitudinal pitch, i.e. spacing between tubes of two adjacent tube
        rows [m].
    arrangement: {'inline', 'staggered' (default)}
        Indicates if the tube bank either has an inline or staggered arrangement
        of the tubes.

    Returns
    -------
    Total number of tubes in the heat exchanger [-].
    """
    if arrangement == 'inline':
        N_t = (L2 * L3) / (S_t * S_l)
    else:
        N_t = (L3 / S_t) * (L2 / S_l + 1) * 0.5
        N_t += (L3 / S_t - 1) * (L2 / S_l - 1) * 0.5
    return N_t.to('m / m').m


class TubeBankInside:
    """Defines the internal geometry of the tube bank.

    Parameters
    ----------
    S_t:
        Lateral or transverse pitch, i.e. spacing between tubes of the
        same row [m].
    S_l:
        Longitudinal pitch, i.e. spacing between tubes of two adjacent tube
        rows [m].
    D_i:
        Inside diameter of the tubes [m].
    alpha: optional
        Ratio of total heat transfer surface area to total volume of the heat
        exchanger [1 / m].
    L1: optional
        Tube length in the direction of inside flow, i.e. the length of the
        tube available for heat transfer with the outside flow [m].
    L2: optional
        Length of the tube bank parallel to the direction of external flow
        [m].
    L3: optional
        Length of the tube bank perpendicular to the direction of external
        flow [m].
    t_header:
        Header thickness; L1 + 2 * t_header is the tube length that needs
        to be taken into account for determining pressure drop along one
        tube [m].
    tube_arrangement: {'inline', 'staggered' (default)}
        Indicates if the tube bank either has an inline or staggered
        arrangement of the tubes.
    """
    def __init__(
        self,
        S_t: Quantity,
        S_l: Quantity,
        D_i: Quantity,
        alpha: Quantity | None = None,
        L1: Quantity | None = None,
        L2: Quantity | None = None,
        L3: Quantity | None = None,
        t_header: Quantity = Q_(0.0, 'm'),
        tube_arrangement: str = 'staggered'
    ) -> None:
        self.S_t = S_t.to('m')
        self.S_l = S_l.to('m')
        self.D_i = D_i.to('m')
        self._alpha = alpha.to('1 / m') if alpha is not None else None
        self.L1 = L1.to('m') if L1 is not None else None
        self.L2 = L2.to('m') if L2 is not None else None
        self.L3 = L3.to('m') if L3 is not None else None
        self.t_header = t_header.to('m')
        self.tube_arrangement = tube_arrangement

    @property
    def N_t(self) -> float | None:
        """Total number of tubes in the heat exchanger [-]."""
        if (self.L2 is not None) and (self.L3 is not None):
            return total_number_of_tubes(
                self.L2, self.L3, self.S_t,
                self.S_l, self.tube_arrangement
            )

    @property
    def A(self) -> Quantity | None:
        """Total heat transfer area [m ** 2]."""
        if (self.L1 is not None) and (self.N_t is not None):
            A = np.pi * self.D_i * self.L1 * self.N_t
            return A.to('m ** 2')

    @property
    def A_o(self) -> Quantity | None:
        """Minimum free flow area [m ** 2]."""
        if self.N_t is not None:
            A_o = np.pi * (self.D_i ** 2) / 4 * self.N_t
            return A_o.to('m ** 2')

    @property
    def N_r1(self) -> float | None:
        """Number of tubes of 1st row."""
        if self.L3 is not None:
            N_r1 = self.L3 / self.S_t
            return N_r1.to('m / m').m

    @property
    def N_r(self) -> float | None:
        """Number of rows."""
        if self.L2 is not None:
            N_r = self.L2 / self.S_l
            return N_r.to('m / m').m

    @property
    def A_fr(self) -> Quantity | None:
        """Core frontal area [m ** 2]."""
        if all([self.L2, self.L3]):
            A_fr = self.L2 * self.L3
            return A_fr.to('m ** 2')

    @property
    def D_h(self) -> Quantity:
        """Hydraulic diameter [m]."""
        return self.D_i.to('m')

    @property
    def V(self) -> Quantity | None:
        """Heat exchanger volume [m ** 3]."""
        if all([self.L1, self.L2, self.L3]):
            V = self.L1 * self.L2 * self.L3
            return V.to('m ** 3')

    @property
    def sigma(self) -> float:
        """Ratio of free flow to frontal area."""
        if all([self.L2, self.L3, self.N_t]):
            sigma = (np.pi / 4) * self.D_i ** 2 * self.N_t / (self.L2 * self.L3)
        else:
            sigma = (np.pi / 4) * self.D_i ** 2 / (self.S_t * self.S_l)
        return sigma.to('m ** 2 / m ** 2').m

    @property
    def alpha(self) -> Quantity | None:
        """Ratio of total heat transfer area to total volume of the heat
        exchanger core [1/m].
        """
        if (self.A is not None) and (self.V is not None):
            alpha = self.A / self.V
        else:
            alpha = self._alpha
        return alpha.to('1 / m')

    @property
    def heat_transfer_length(self) -> Quantity | None:
        """Tube length for heat transfer [m]."""
        return self.L1.to('m')

    @property
    def pressure_drop_length(self) -> Quantity | None:
        """Tube length for pressure drop [m]."""
        if self.L1 is not None:
            L = self.L1 + 2 * self.t_header
            return L.to('m')


class TubeBankOutside(ABC):
    """Defines the general external geometry of a tube bank.

    Parameters
    ----------
    S_t:
        Lateral or transverse pitch, i.e. the spacing between tubes of the
        same row [m].
    S_l:
        Longitudinal pitch, i.e. the spacing between tubes of two adjacent tube
        rows [m].
    alpha: optional
        Ratio of total heat transfer surface area to total volume of the heat
        exchanger [1 / m].
    sigma: optional
        Ratio of the free flow area to the frontal area of the heat exchanger [-].
    L1: optional
        Tube length in the direction of internal flow, i.e. the length of the
        tube available for heat transfer with external flow or the width of
        the heat exchanger core [m].
    L2: optional
        Length of the heat exchanger core parallel to the direction of external
        flow, i.e. the depth of the heat exchanger core [m].
    L3: optional
        Length of the heat exchanger core perpendicular to the direction of
        external flow, i.e. the height of the heat exchanger core [m].
    """
    def __init__(
        self,
        S_t: Quantity,
        S_l: Quantity,
        alpha: Quantity | None = None,
        sigma: float | None = None,
        L1: Quantity | None = None,
        L2: Quantity | None = None,
        L3: Quantity | None = None
    ) -> None:
        self.S_t = S_t.to('m')
        self.S_l = S_l.to('m')
        self._alpha = alpha.to('1 / m') if alpha is not None else None
        self._sigma = sigma
        self.L1 = L1.to('m') if L1 is not None else None
        self.L2 = L2.to('m') if L2 is not None else None
        self.L3 = L3.to('m') if L3 is not None else None

    @property
    @abstractmethod
    def N_t(self) -> float | None:
        ...

    @property
    @abstractmethod
    def A(self) -> Quantity | None:
        """Total heat transfer area [m ** 2]."""
        ...

    @property
    @abstractmethod
    def A_o(self) -> Quantity | None:
        """Minimum free flow area [m ** 2]."""
        ...

    @property
    def N_r1(self) -> float | None:
        """Number of tubes of 1st row."""
        if self.L3 is not None:
            N_r1 = self.L3 / self.S_t
            return N_r1.to('m / m').m

    @property
    def N_r(self) -> float | None:
        """Number of rows."""
        if self.L2 is not None:
            N_r = self.L2 / self.S_l
            return N_r.to('m / m').m

    @property
    def A_fr(self) -> Quantity | None:
        """Frontal area [m ** 2]."""
        if (self.L1 is not None) and (self.L3 is not None):
            A_fr = self.L1 * self.L3
            return A_fr.to('m ** 2')

    @property
    def V(self) -> Quantity | None:
        """Heat exchanger volume [m ** 3]."""
        if all([self.L1, self.L2, self.L3]):
            V = self.L1 * self.L2 * self.L3
            return V.to('m ** 3')

    @property
    def D_h(self) -> Quantity | None:
        """Hydraulic diameter [m]."""
        if all([self.A_o, self.L2, self.A]):
            D_h = hydraulic_diameter(self.A_o, self.L2, self.A)
        elif (self._alpha is not None) and (self._sigma is not None):
            D_h = hydraulic_diameter(sigma=self._sigma, alpha=self._alpha)
        else:
            D_h = None
        return D_h

    @property
    def sigma(self) -> float | None:
        """Ratio of free flow area to frontal area [-]."""
        if (self.A_o is not None) and (self.A_fr is not None):
            sigma = (self.A_o / self.A_fr).to('m ** 2 / m ** 2').m
        else:
            sigma = self._sigma
        return sigma

    @property
    def alpha(self) -> Quantity | None:
        """Ratio of total heat transfer area to total volume of the heat exchanger
        core [1/m].
        """
        if (self.A is not None) and (self.V is not None):
            alpha = self.A / self.V
        else:
            alpha = self._alpha
        return alpha.to('1 / m')

    @property
    def pressure_drop_length(self) -> Quantity | None:
        """Flow length for pressure drop [m]."""
        return self.L2.to('m')


TTubeBankOutside = TypeVar('TTubeBankOutside', bound=TubeBankOutside)


class BareInlineTBO(TubeBankOutside):
    """Defines the external geometry of a bare inline tube bank.

    Parameters
    ----------
    S_t:
        Lateral or transverse pitch, i.e. spacing between tubes of the
        same row [m].
    S_l:
        Longitudinal pitch, i.e. spacing between tubes of two adjacent tube
        rows [m].
    D_o:
        Outside diameter of the tubes [m].
    alpha: optional
        Ratio of total heat transfer surface area to total volume of the heat
        exchanger [1 / m].
    sigma: optional
        Ratio of the free flow area to the frontal area of the heat exchanger [-].
    L1: optional
        Tube length in the direction of inside flow, i.e. the length of the
        tube available for heat transfer with the outside flow [m].
    L2: optional
        Length of the tube bank parallel to the direction of external flow
        [m].
    L3: optional
        Length of the tube bank perpendicular to the direction of external
        flow [m].
    """
    def __init__(
        self,
        S_t: Quantity,
        S_l: Quantity,
        D_o: Quantity,
        alpha: Quantity | None = None,
        sigma: float | None = None,
        L1: Quantity | None = None,
        L2: Quantity | None = None,
        L3: Quantity | None = None
    ) -> None:
        super().__init__(S_t, S_l, alpha, sigma, L1, L2, L3)
        self.D_o = D_o.to('m')

    @property
    def N_t(self) -> float | None:
        if all([self.L2, self.L3]):
            return total_number_of_tubes(
                self.L2,
                self.L3,
                self.S_t,
                self.S_l,
                arrangement='inline'
            )

    @property
    def A(self) -> Quantity | None:
        """Total heat transfer area [m ** 2]."""
        if all([self.L1, self.N_t, self.L2, self.L3]):
            A_tb = np.pi * self.D_o * self.L1 * self.N_t
            A_hd = self.L2 * self.L3 - np.pi * self.D_o ** 2 / 4 * self.N_t
            A = A_tb + 2 * A_hd
            return A.to('m ** 2')

    @property
    def A_o(self) -> Quantity | None:
        """Minimum free flow area [m ** 2]."""
        if (self.L1 is not None) and (self.N_r1 is not None):
            A_o = (self.S_t - self.D_o) * self.L1 * self.N_r1
            return A_o.to('m ** 2')


class BareStaggeredTBO(BareInlineTBO):
    """Defines the outside geometry of the bare staggered tube bank.

    Parameters
    ----------
    See class `BareInlineTBO`.
    """
    @property
    def N_t(self) -> float | None:
        if all([self.L2, self.L3]):
            return total_number_of_tubes(
                self.L2,
                self.L3,
                self.S_t,
                self.S_l,
                arrangement='staggered'
            )

    @property
    def A_o(self) -> Quantity | None:
        """Minimum free flow area [m ** 2]."""
        if all([self.L1, self.L3]):
            x = (self.S_t - self.D_o) / 2
            y = np.sqrt((self.S_t / 2) ** 2 + self.S_l ** 2) - self.D_o
            z = min(2 * x, 2 * y)
            A_o = (
                (self.L3 / self.S_t - 1) * z
                + (self.S_t - self.D_o)
            ) * self.L1
            return A_o.to('m ** 2')


class CircularFinInlineTBO(TubeBankOutside):
    """Defines the external geometry of a circular finned inline tube bank.

    Parameters
    ----------
    See class `TubeBankOutside`

    Other Parameters
    ----------------
    D_o:
        Outer tube diameter [m]
    D_r:
        Effective diameter of the root of a circular fin. Depending upon the
        manufacturing techniques, it may be the tube outside diameter or tube
        outside diameter plus two times the fin collar thickness [m].
    D_f:
        Fin tip diameter [m].
    t_f:
        Thickness of a circular fin [m].
    N_f:
        The number of fins per unit length [1/m]
    """
    def __init__(
        self,
        S_t: Quantity,
        S_l: Quantity,
        D_o: Quantity,
        D_r: Quantity,
        D_f: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        alpha: Quantity | None = None,
        sigma: float | None = None,
        L1: Quantity | None = None,
        L2: Quantity | None = None,
        L3: Quantity | None = None
    ) -> None:
        super().__init__(S_t, S_l, alpha, sigma, L1, L2, L3)
        self.D_o = D_o.to('m')
        self.D_r = D_r.to('m')
        self.D_f = D_f.to('m')
        self.t_f = t_f.to('m')
        self.N_f = N_f.to('1 / m')

    @property
    def N_t(self) -> float | None:
        if all([self.L2, self.L3]):
            return total_number_of_tubes(
                self.L2,
                self.L3,
                self.S_t,
                self.S_l,
                arrangement='inline'
            )

    @property
    def A_p(self) -> Quantity | None:
        """Prime area, i.e. the tube surface area minus the area blocked by the
        fins [m ** 2].
        """
        if all([self.L1, self.L2, self.L3, self.N_t]):
            A_p = np.pi * self.D_r * self.L1 * (1 - self.t_f * self.N_f) * self.N_t
            A_p += 2 * (self.L2 * self.L3 - np.pi * self.D_o ** 2 / 4 * self.N_t)
            return A_p.to('m ** 2')

    @property
    def A_f(self) -> Quantity | None:
        """Fin surface area [m ** 2]."""
        if all([self.L1, self.N_t]):
            A_f = (
                2 * np.pi * (self.D_f ** 2 - self.D_r ** 2) / 4
                + np.pi * self.D_f * self.t_f
            )
            A_f *= self.N_f * self.L1 * self.N_t
            return A_f.to('m ** 2')

    @property
    def A(self) -> Quantity | None:
        """Total heat transfer area [m ** 2]."""
        if all([self.A_p, self.A_f]):
            A = self.A_p + self.A_f
            return A.to('m ** 2')

    @property
    def A_f_to_A(self) -> float:
        """Ratio of fin area to total heat transfer area."""
        # prime area per 1 m (L1 = 1 m) per tube (N_t = 1)
        A_p = np.pi * self.D_r * Q_(1, 'm') * (1 - self.t_f * self.N_f) * 1
        A_p += 2 * (self.S_l * self.S_t - np.pi * self.D_o ** 2 / 4 * 1)
        # fin area per 1 m (L1 = 1 m) per tube (N_t = 1)
        A_f = (
            2 * np.pi * (self.D_f ** 2 - self.D_r ** 2) / 4
            + np.pi * self.D_f * self.t_f
        )
        A_f *= self.N_f * Q_(1, 'm') * 1
        A_to_A_f = 1 + A_p / A_f
        A_f_to_A = 1 / A_to_A_f
        return A_f_to_A.to('m ** 2 / m ** 2').m

    @property
    def A_o(self) -> Quantity | None:
        """Minimum free flow area [m ** 2]."""
        if all([self.L1, self.L3]):
            A_o = (
                (self.S_t - self.D_r) * self.L1
                - (self.D_f - self.D_r) * self.t_f * self.N_f * self.L1
            ) * (self.L3 / self.S_t)
            return A_o.to('m ** 2')


class CircularFinStaggeredTBO(CircularFinInlineTBO):
    """Defines the external geometry of a circular finned staggered tube bank.

    Parameters
    ----------
    See class `CircularFinInlineTBO`.
    """
    @property
    def N_t(self) -> float | None:
        if all([self.L2, self.L3]):
            return total_number_of_tubes(
                self.L2, self.L3, self.S_t,
                self.S_l, arrangement='staggered'
            )

    @property
    def A_o(self) -> Quantity | None:
        """Minimum free flow area [m ** 2]."""
        if all([self.L1, self.L3]):
            x = (
                (self.S_t - self.D_r)
                - (self.D_f - self.D_r) * self.t_f * self.N_f
            ) / 2
            S_d = np.sqrt((self.S_t / 2) ** 2 + self.S_l ** 2)
            y = S_d - self.D_r - (self.D_f - self.D_r) * self.t_f * self.N_f
            z = min(2 * x, 2 * y)
            A_o = (
                (self.L3 / self.S_t - 1) * z
                + (self.S_t - self.D_r)
                - (self.D_f - self.D_r) * self.t_f * self.N_f
            ) * self.L1
            return A_o.to('m ** 2')


class PlainFinInlineTBO(TubeBankOutside):
    """Defines the external geometry of a continuous, plain finned, inline tube
    bank.

    Parameters
    ----------
    See class `TubeBankOutside`.

    Other Parameters
    ----------------
    D_o:
        Outer tube diameter [m]
    D_r:
        Effective diameter of the root of a circular fin. Depending upon the
        manufacturing techniques, it may be the tube outside diameter or tube
        outside diameter plus two times the fin collar thickness [m].
    t_f:
        Thickness of a circular fin [m].
    N_f:
        The number of fins per unit length [1/m]
    """
    def __init__(
        self,
        S_t: Quantity,
        S_l: Quantity,
        D_o: Quantity,
        D_r: Quantity,
        t_f: Quantity,
        N_f: Quantity,
        alpha: Quantity | None = None,
        sigma: float | None = None,
        L1: Quantity | None = None,
        L2: Quantity | None = None,
        L3: Quantity | None = None
    ) -> None:
        super().__init__(S_t, S_l, alpha, sigma, L1, L2, L3)
        self.D_o = D_o.to('m')
        self.D_r = D_r.to('m')
        self.t_f = t_f.to('m')
        self.N_f = N_f.to('1 / m')

    @property
    def N_t(self) -> float | None:
        if all([self.L2, self.L3]):
            return total_number_of_tubes(
                self.L2,
                self.L3,
                self.S_t,
                self.S_l,
                arrangement='inline'
            )

    @property
    def A_p(self) -> Quantity | None:
        """Prime area, i.e. the tube surface area minus the area blocked by the
        fins [m ** 2].
        """
        if all([self.L1, self.N_t, self.L2, self.L3]):
            A_p = np.pi * self.D_r * self.L1 * (1 - self.t_f * self.N_f) * self.N_t
            A_p += 2 * (self.L2 * self.L3 - np.pi * self.D_o ** 2 / 4 * self.N_t)
            return A_p.to('m ** 2')

    @property
    def A_f(self) -> Quantity | None:
        """Fin surface area [m ** 2]."""
        if all([self.L1, self.N_t, self.L2, self.L3]):
            a = 2 * (self.L2 * self.L3 - (np.pi * self.D_r ** 2) / 4 * self.N_t)
            a *= self.N_f * self.L1
            b = 2 * self.L3 * self.t_f * self.N_f * self.L1  # leading and trailing edges area
            A_f = a + b
            return A_f.to('m ** 2')

    @property
    def A(self) -> Quantity | None:
        """Total heat transfer area [m ** 2]."""
        if all([self.A_p, self.A_f]):
            A = self.A_p + self.A_f
            return A.to('m ** 2')

    @property
    def A_f_to_A(self) -> float:
        """Ratio of fin area to total heat transfer area."""
        # prime area per 1 m (L1 = 1 m) per tube (N_t = 1)
        A_p = np.pi * self.D_r * Q_(1, 'm') * (1 - self.t_f * self.N_f) * 1
        A_p += 2 * (self.S_l * self.S_t - np.pi * self.D_o ** 2 / 4 * 1)
        # fin area per 1 m (L1 = 1 m) per tube (N_t = 1)
        a = 2 * (self.S_l * self.S_t - (np.pi * self.D_r ** 2) / 4 * 1)
        a *= self.N_f * Q_(1, 'm')
        b = 2 * self.S_t * self.t_f * self.N_f * Q_(1, 'm')  # leading and trailing edges area
        A_f = a + b
        A_to_A_f = 1 + A_p / A_f
        A_f_to_A = 1 / A_to_A_f
        return A_f_to_A.to('m ** 2 / m ** 2').m

    @property
    def A_o(self) -> Quantity | None:
        """Minimum free flow area [m ** 2]."""
        if all([self.L1, self.L3]):
            a = (self.S_t - self.D_r) * self.L1
            b = (self.S_t - self.D_r) * self.t_f * self.N_f * self.L1
            A_o = (a - b) * (self.L3 / self.S_t)
            return A_o.to('m ** 2')


class PlainFinStaggeredTBO(PlainFinInlineTBO):
    """Defines the external geometry of a continuous, plain finned, staggered
    tube bank.

    Parameters
    ----------
    See class `PlainFinInlineTBO`.
    """
    @property
    def N_t(self) -> float | None:
        if all([self.L2, self.L3]):
            return total_number_of_tubes(
                self.L2,
                self.L3,
                self.S_t,
                self.S_l,
                arrangement='staggered'
            )

    @property
    def A_o(self) -> Quantity | None:
        """Minimum free flow area [m ** 2]."""
        if all([self.L1, self.L3]):
            x = 0.5 * (self.S_t - self.D_r) * (1 - self.t_f * self.N_f)
            S_d = np.sqrt((self.S_t / 2) ** 2 + self.S_l ** 2)
            y = S_d - self.D_r - (self.S_t - self.D_r) * self.t_f * self.N_f
            z = min(2 * x, 2 * y)
            a = (self.L3 / self.S_t - 1) * z
            b = self.S_t - self.D_r
            c = (self.S_t - self.D_r) * self.t_f * self.N_f
            A_o = (a + b - c) * self.L1
            return A_o.to('m ** 2')

    @property
    def zeta(self) -> float:
        """Contraction and expansion area ratio, needed for calculating the
        pressure losses at the leading and trailing fin edges due to sudden
        contraction at the entrance and sudden expansion at the exit of the
        heat exchanger.
        """
        zeta = 1 - (self.t_f * self.N_f).to('m / m').m
        return zeta
