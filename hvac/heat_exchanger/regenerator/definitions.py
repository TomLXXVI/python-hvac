"""
Implementation of definitions used for the heat transfer analysis of rotary
regenerators.

References
----------
Shah, R. K., & Sekulic, D. P. (2003). Fundamentals of Heat Exchanger Design.
John Wiley & Sons. Chapter 5: Thermal design theory for regenerators, § 5.1.2.
"""
import numpy as np
from hvac import Quantity

Q_ = Quantity


class Fluid:

    @staticmethod
    def heat_capacity_rate(m_dot: Quantity, cp: Quantity) -> Quantity:
        """Returns the heat capacity rate of the fluid.

        Parameters
        ----------
        m_dot:
            Mass flow rate of fluid.
        cp:
            Specific heat of fluid.
        """
        C_dot = m_dot * cp
        return C_dot.to('J / (s * K)')

    @staticmethod
    def _heat_capacitance_01(m: Quantity, cp: Quantity) -> Quantity:
        """Returns the heat capacitance of the fluid.

        Parameters
        ----------
        m:
            Mass of the fluid contained in the regenerator matrix at any instant
            of time.
        cp:
            Specific heat of the fluid.
        """
        C = m * cp
        return C.to('J / K')

    @staticmethod
    def _heat_capacitance_02(C_dot: Quantity, t_d: Quantity) -> Quantity:
        """Returns the heat capacitance of the fluid.

        Parameters
        ----------
        C_dot:
            Heat capacity rate of the fluid.
        t_d:
            Fluid dwell time (or residence time).
        """
        C = C_dot * t_d
        return C.to('J / K')

    @staticmethod
    def _heat_capacitance_03(
        C_dot: Quantity,
        L: Quantity,
        u_m: Quantity
    ) -> Quantity:
        """Returns the heat capacitance of the fluid.

        Parameters
        ----------
        C_dot:
            Heat capacity rate of the fluid.
        L:
            Length of the regenerator matrix.
        u_m:
            Mean fluid axial velocity.
        """
        C = C_dot * L / u_m
        return C.to('J / K')

    @staticmethod
    def heat_capacitance(
        m: Quantity | None = None,
        cp: Quantity | None = None,
        C_dot: Quantity | None = None,
        t_d: Quantity | None = None,
        L: Quantity | None = None,
        u_m: Quantity | None = None
    ) -> Quantity | None:
        """Returns the heat capacitance of the fluid.

        Parameters
        ----------
        m:
            Mass of the fluid contained in the regenerator matrix at any
            instant of time.
        cp:
            Specific heat of the fluid.
        C_dot:
            Heat capacity rate of the fluid.
        t_d:
            Fluid dwell time (or residence time).
        L:
            Length of the regenerator matrix.
        u_m:
            Mean fluid axial velocity.

        Valid groups of keyword arguments are:
        * `m` and `cp`, or
        * `C_dot` and `t_d`, or
        * `C_dot`, `L`, and `u_m`.
        """
        if all((m, cp)) is not False:
            C = Fluid._heat_capacitance_01(m, cp)
            return C
        elif all((C_dot, t_d)) is not False:
            C = Fluid._heat_capacitance_02(C_dot, t_d)
            return C
        elif all((C_dot, L, u_m)) is not False:
            C = Fluid._heat_capacitance_03(C_dot, L, u_m)
            return C
        else:
            return None


class Matrix:

    @staticmethod
    def heat_capacitance(m_w: Quantity, c_w: Quantity) -> Quantity:
        """Returns the heat capacitance of the regenerator matrix wall.

        Parameters
        ----------
        m_w:
            Mass of all matrices (disks).
        c_w:
            Specific heat of the matrix material.
        """
        C_r = m_w * c_w
        return C_r.to('J / K')

    @staticmethod
    def _heat_capacity_rate_01(m_w: Quantity, c_w: Quantity, N: Quantity) -> Quantity:
        """Returns the matrix wall heat capacity rate for a rotary regenerator.

        Parameters
        ----------
        m_w:
            Mass of all matrices (disks).
        c_w:
            Specific heat of the matrix material.
        N:
            Rotational speed of the rotary regenerator.
        """
        C_r = Matrix.heat_capacitance(m_w, c_w)
        C_dot_r = C_r * N
        return C_dot_r.to('J / (K * s)')

    @staticmethod
    def heat_capacity_rate_02(C_r: Quantity, N: Quantity) -> Quantity:
        """Returns the matrix wall heat capacity rate during the cold/hot
        period.

        Parameters
        ----------
        C_r:
            Total matrix heat capacitance.
        N:
            Rotational speed of the rotary regenerator.
        """
        C_rj_dot = C_r * N
        return C_rj_dot.to('J / (K * s)')

    @staticmethod
    def heat_capacity_rate(
        m_w: Quantity | None = None,
        c_w: Quantity | None = None,
        N: Quantity | None = None,
        C_r: Quantity | None = None
    ) -> Quantity | None:
        """Returns the matrix wall heat capacity rate during the cold/hot
        period.

        Parameters
        ----------
        m_w:
            Mass of all matrices (disks).
        c_w:
            Specific heat of the matrix material.
        N:
            Rotational speed of the rotary regenerator.
        C_r:
            Total matrix heat capacitance.

        Valid groups of keyword arguments are:
        * `m_w`, `c_w`, and `N`, or
        * `C_r` and `N`.
        """
        if all((m_w, c_w, N)) is not False:
            C_r_dot = Matrix._heat_capacity_rate_01(m_w, c_w, N)
            return C_r_dot
        elif all((C_r, N)) is not False:
            C_r_dot = Matrix.heat_capacity_rate_02(C_r, N)
            return C_r_dot
        else:
            return None

    @staticmethod
    def sector_heat_capacitance(C_r: Quantity, theta_j: Quantity) -> Quantity:
        """Returns the heat capacitance of the cold/hot sector of the
        regenerator matrix wall.

        Parameters
        ----------
        C_r:
            Total matrix heat capacitance.
        theta_j:
            Disk sector angle through which cold/hot gas flows.
        """
        theta_t = Q_(2 * np.pi, 'rad')
        C_rj = C_r * (theta_j.to('rad') / theta_t)
        return C_rj.to('J / K')

    @staticmethod
    def _sector_heat_transfer_area_01(
        A: Quantity,
        theta_j: Quantity
    ) -> Quantity:
        """Returns the heat transfer area of the cold/hot sector of the
        regenerator matrix wall.

        Parameters
        ----------
        A:
            Total heat transfer area of rotary regenerator.
        theta_j:
            Disk sector angle through which cold/hot gas flows.
        """
        if theta_j.units in ('deg', 'rad'):
            theta_j = theta_j.to('rad')
            theta_t = Q_(2 * np.pi, 'rad')
            r_j = theta_j.to('rad') / theta_t
        else:
            r_j = theta_j.to('frac')
        A_j = A * r_j
        return A_j.to('m**2')

    @staticmethod
    def _sector_heat_transfer_area_02(
        beta: Quantity,
        V: Quantity,
        theta_j: Quantity
    ) -> Quantity:
        """Returns the heat transfer area of the cold/hot sector of the
        regenerator matrix wall.

        Parameters
        ----------
        beta:
            Heat transfer surface area density or packing density
            (i.e. heat transfer surface area per unit of volume).
        V:
            Total volume of all matrices.
        theta_j:
            Disk sector angle through which cold/hot gas flows.
        """
        if theta_j.units in ('deg', 'rad'):
            theta_j = theta_j.to('rad')
            theta_t = Q_(2 * np.pi, 'rad')
            r_j = theta_j.to('rad') / theta_t
        else:
            r_j = theta_j.to('frac')
        A_j = beta * V * r_j
        return A_j

    @staticmethod
    def sector_heat_transfer_area(
        A: Quantity | None = None,
        theta_j: Quantity | None = None,
        beta: Quantity | None = None,
        V: Quantity | None = None
    ) -> Quantity | None:
        """Returns the heat transfer area of the cold/hot sector of the
        regenerator matrix wall.

        Parameters
        ----------
        A:
            Total heat transfer area of rotary regenerator.
        theta_j:
            Disk sector angle through which cold/hot gas flows.
        beta:
            Heat transfer surface area density or packing density
            (i.e. heat transfer surface area per unit of volume).
        V:
            Total volume of all matrices.

        Valid groups of keyword arguments are:
        * `A` and `theta_j`, or
        * `beta`, `V`, and `theta_j`.
        """
        if all((A, theta_j)) is not False:
            A_j = Matrix._sector_heat_transfer_area_01(A, theta_j)
            return A_j
        elif all((beta, V, theta_j)) is not False:
            A_j = Matrix._sector_heat_transfer_area_02(beta, V, theta_j)
            return A_j
        else:
            return None

    @staticmethod
    def total_heat_transfer_area(
        A_fr: Quantity,
        L_w: Quantity,
        beta: Quantity,
        seal_fraction: Quantity
    ) -> Quantity:
        """Returns the total matrix heat transfer surface area.

        Parameters
        ----------
        A_fr:
            Frontal area of matrix
        L_w:
            Rotor length.
        beta:
            Packing density (area per unit of volume).
        seal_fraction:
            Fraction of rotor face area covered by radial seals.
        """
        A_fr = A_fr.to('m**2')
        L_w = L_w.to('m')
        beta = beta.to('1 / m')
        V_w = A_fr * L_w
        A_w = beta * V_w * (1.0 - seal_fraction.to('frac'))
        return A_w.to('m**2')

    @staticmethod
    def _porosity_01(A_o: Quantity, A_fr: Quantity) -> Quantity:
        """The core or matrix porosity is the ratio of the flow area to the
        frontal area of the core in case the heat transfer surface is made of
        continuous flow passages.

        Parameters
        ----------
        A_o:
            Total flow area of the matrix.
        A_fr:
            Frontal area of the matrix.
        """
        sigma = A_o.to('m**2') / A_fr.to('m**2')
        return sigma

    @staticmethod
    def _porosity_02(V_void: Quantity, A_fr: Quantity, L: Quantity) -> Quantity:
        """The core or matrix porosity is a ratio of the void volume to the
        total core or matrix volume in case the heat transfer surface has
        interruptions (such as perforations) or is made up of porous materials.

        Parameters
        ----------
        V_void:
            Void volume.
        A_fr:
            Frontal area of the matrix.
        L:
            Length of the regenerator matrix.
        """
        sigma = V_void.to('m**3') / (A_fr * L).to('m**3')
        return sigma

    @staticmethod
    def porosity(
        A_o: Quantity | None = None,
        A_fr: Quantity | None = None,
        V_void: Quantity | None = None,
        L: Quantity | None = None
    ) -> Quantity | None:
        """Returns the core or matrix porosity.
        - The core or matrix porosity is the ratio of the flow area to the
        frontal area of the core in case the heat transfer surface is made of
        continuous flow passages.
        - The core or matrix porosity is a ratio of the void volume to the
        total core or matrix volume in case the heat transfer surface has
        interruptions (such as perforations) or is made up of porous materials.

        Parameters
        ----------
        A_o:
            Total flow area of the matrix.
        A_fr:
            Frontal area of the matrix.
        V_void:
            Void volume.
        L:
            Length of the regenerator matrix.

        Valid groups of keyword arguments are:
        * `A_o` and `A_fr`, or
        * `V_void`, `A_fr`, and `L`.
        """
        if all((A_o, A_fr)) is not False:
            sigma = Matrix._porosity_01(A_o, A_fr)
            return sigma
        elif all((V_void, A_fr, L)) is not False:
            sigma = Matrix._porosity_02(V_void, A_fr, L)
            return sigma
        else:
            return None

    @staticmethod
    def packing_density(
        A: Quantity | None = None,
        V: Quantity | None = None,
        A_fr: Quantity | None = None,
        L: Quantity | None = None
    ) -> Quantity:
        """Returns the heat transfer surface area density or packing density
        (i.e. the heat transfer surface area per unit of volume).

        Parameters
        ----------
        A:
            Total heat transfer area of rotary regenerator.
        V:
            Total volume of all matrices.
        A_fr:
            Frontal area of the matrix.
        L:
            Length of the regenerator matrix.

        Valid groups of keyword arguments are:
        * `A` and `V`, or
        * `A`, `A_fr`, and `L`.
        """
        if all((A, V)) is not False:
            beta = A / V
            return beta.to('1 / m')
        elif all((A, A_fr, L)) is not False:
            beta = A / (A_fr * L)
            return beta.to('1 / m')
        else:
            return None

    @staticmethod
    def hydraulic_diameter(sigma: Quantity, beta: Quantity) -> Quantity:
        """Returns the hydraulic diameter of the flow passages.

        Parameters
        ----------
        sigma:
            Core or matrix porosity.
        beta:
            Heat transfer surface area density or packing density.
        """
        D_h = 4 * sigma / beta
        return D_h.to('m')

    @staticmethod
    def mass(
        A_fr: Quantity,
        L_w: Quantity,
        rho_w: Quantity,
        sigma: Quantity
    ) -> Quantity:
        """Returns the mass of the matrix.

        Parameters
        ----------
        A_fr:
            Frontal area of matrix
        L_w:
            Rotor length.
        rho_w:
            Density of matrix material.
        sigma:
            Matrix porosity.
        """
        A_fr = A_fr.to('m**2')
        L_w = L_w.to('m')
        rho_w = rho_w.to('kg / m**3')
        sigma = sigma.to('frac')
        V_w = A_fr * L_w * (1.0 - sigma)
        m_w = rho_w * V_w
        return m_w.to('kg')
