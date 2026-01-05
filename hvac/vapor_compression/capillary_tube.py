"""
Calculating the length of a capillary tube.

References
----------
Stoecker, W. F., & Jones, J. W. (1982). Refrigeration and Air Conditioning.
McGraw-Hill Publishing Company.
"""
import math
from hvac import Quantity
from hvac.fluids import Fluid, FluidState

import warnings
from hvac.fluids import CoolPropWarning

warnings.filterwarnings("ignore", category=CoolPropWarning)

Q_ = Quantity
PI = math.pi


def reynolds(vel: Quantity, d_i: Quantity, fluid_state: FluidState) -> float:
    """
    Returns the Reynolds number of the fluid flow through the capillary tube.
    """
    vel = vel.to('m / s')
    d_i = d_i.to('m')
    rho = fluid_state.rho.to('kg / m ** 3')
    mu = fluid_state.mu.to('Pa * s')
    Re = vel * d_i * rho / mu
    return Re.to('frac').m


def friction_factor(Re: float) -> float:
    """
    Returns the friction factor of the fluid valid for Reynolds numbers in the
    lower range of the turbulent region.
    """
    return 0.33 / (Re ** 0.25)


def velocity(m_dot: Quantity, d_i: Quantity, fluid_state: FluidState) -> Quantity:
    """
    Returns the flowing velocity of the refrigerant.
    """
    m_dot = m_dot.to('kg / s')
    d_i = d_i.to('m')
    rho = fluid_state.rho.to('kg / m ** 3')
    A = PI * (d_i ** 2) / 4
    vel = m_dot / (A * rho)
    return vel


def m_dot_to_A(m_dot: Quantity, d_i: Quantity) -> Quantity:
    """
    Returns the ratio of the refrigerant mass flow rate to the cross-sectional
    area of the tube (= mass velocity).
    """
    m_dot = m_dot.to('kg / s')
    d_i = d_i.to('m')
    A = PI * (d_i ** 2) / 4
    _m_dot_to_A = m_dot / A
    return _m_dot_to_A.to('kg / (s * m ** 2')


class CapillaryTube:
    """
    Class to calculate the required length of a capillary tube given:
    -   the mass flow rate of refrigerant
    -   the condensing temperature
    -   the evaporating temperature
    -   the inside tube diameter
    """
    def __init__(
        self,
        d_i: Quantity,
        refrigerant: Fluid,
        m_dot: Quantity,
        T_cnd: Quantity,
        T_evp: Quantity,
        dT_step: Quantity = Q_(1, 'K')
    ) -> None:
        """
        Creates a `CapillaryTube` calculation object.

        Parameters
        ----------
        d_i: Quantity
            Internal diameter of the capillary tube.
        refrigerant: Fluid
            The refrigerant that flows through the capillary tube.
        m_dot: Quantity
            Mass flow rate of refrigerant.
        T_cnd: Quantity
            Condensing temperature.
        T_evp: Quantity
            Evaporating temperature.
        dT_step: Quantity, default 1 K
            The temperature step taken between two iterations of the calculation
            loop.
        """
        self.d_i = d_i.to('m')
        self.refrigerant = refrigerant
        self.m_dot = m_dot.to('kg / s')
        self.T_cnd = T_cnd.to('K')
        self.T_evp = T_evp.to('K')
        self.dT_step = dT_step.to('K')

        self._A = PI * (self.d_i ** 2) / 4
        self._m_dot_to_A = self.m_dot / self._A

    def get_length(self, show_steps: bool = True) -> Quantity:
        """Returns the required length of the capillary tube."""
        i = 0
        T_out, x_out = None, None
        delta_L_list = []
        while True:
            if i == 0:
                T_in = self.T_cnd
                x_in = Q_(0.0, 'frac')
            else:
                T_in = T_out
                x_in = x_out
            T_out = T_in - self.dT_step
            if T_out < self.T_evp:
                break
            x_out, delta_L = self._calculate_step(T_in, x_in, T_out)
            delta_L_list.append(delta_L)

            if show_steps:
                self.print_step(
                    i,
                    T_out,
                    self._rfg_state_out.P,
                    x_out,
                    self._rfg_state_out.rho,
                    self._rfg_state_out.h,
                    self._vel_out,
                    delta_L,
                    sum(delta_L_list)
                )

            i += 1
        L = sum(delta_L_list)
        return L

    def _calculate_step(
        self,
        T_in: Quantity,
        x_in: Quantity,
        T_out: Quantity
    ) -> tuple[Quantity, ...]:
        self._set_conditions_at_entrance(T_in, x_in)
        x_out = self._calculate_vapor_quality(T_out)
        self._set_conditions_at_exit(T_out, x_out)
        delta_L = self._calculate_increment()
        return x_out, delta_L

    def _set_conditions_at_entrance(
        self,
        T_in: Quantity,
        x_in: Quantity
    ) -> None:
        self._rfg_state_in = self.refrigerant(T=T_in.to('K'), x=x_in.to('frac'))
        self._vel_in = velocity(self.m_dot, self.d_i, self._rfg_state_in)
        Re_in = reynolds(self._vel_in, self.d_i, self._rfg_state_in)
        self._f_in = friction_factor(Re_in)

    def _calculate_vapor_quality(
        self,
        T_out: Quantity
    ) -> Quantity:
        rfg_sat_liq = self.refrigerant(T=T_out, x=Q_(0.0, 'frac'))
        rfg_sat_vap = self.refrigerant(T=T_out, x=Q_(1.0, 'frac'))

        rho_sat_vap = rfg_sat_vap.rho.to('kg / m**3').m
        rho_sat_liq = rfg_sat_liq.rho.to('kg / m**3').m
        _m_dot_to_A = self._m_dot_to_A.to('kg / (s * m**2)').m
        h_sat_vap = rfg_sat_vap.h.to('kJ / kg').m
        h_sat_liq = rfg_sat_liq.h.to('kJ / kg').m

        a = (1 / rho_sat_vap - 1 / rho_sat_liq) ** 2
        a *= _m_dot_to_A ** 2
        a *= (1 / 2)

        b1 = 1000 * (h_sat_vap - h_sat_liq)
        b2 = (1 / rho_sat_liq)
        b2 *= (1 / rho_sat_vap - 1 / rho_sat_liq)
        b2 *= _m_dot_to_A ** 2
        b = b1 + b2

        c1 = 1000 * (h_sat_liq - self._rfg_state_in.h.to('kJ / kg').m)
        c2 = (_m_dot_to_A ** 2) * (1 / 2) * ((1 / rho_sat_liq) ** 2)
        c3 = (self._vel_in.to('m / s').m ** 2) / 2
        c = c1 + c2 - c3

        d = math.sqrt(b ** 2 - 4 * a * c)
        x1 = (-b + d) / (2 * a)
        x2 = (-b - d) / (2 * a)
        x = max(x1, x2)
        return Q_(x, 'frac')

    def _set_conditions_at_exit(
        self,
        T_out: Quantity,
        x: Quantity
    ) -> None:
        self._rfg_state_out = self.refrigerant(T=T_out, x=x)
        self._vel_out = self._m_dot_to_A * (1 / self._rfg_state_out.rho)
        Re_out = reynolds(self._vel_out, self.d_i, self._rfg_state_out)
        self._f_out = friction_factor(Re_out)

    def _calculate_increment(self) -> Quantity:
        k1 = self._rfg_state_in.P - self._rfg_state_out.P
        k2 = self.m_dot * (self._vel_out - self._vel_in) / self._A
        k = k1 - k2
        f_mean = (self._f_in + self._f_out) / 2
        vel_mean = (self._vel_in + self._vel_out) / 2
        delta_L = k * (1 / f_mean)
        delta_L *= (2 / vel_mean)
        delta_L *= (1 / self._m_dot_to_A) * self.d_i
        return delta_L.to('m')

    @staticmethod
    def print_step(
        step: int,
        T: Quantity,
        P: Quantity,
        x: Quantity,
        rho: Quantity,
        h: Quantity,
        vel: Quantity,
        delta_L: Quantity | None = None,
        L: Quantity | None = None
    ) -> None:
        print(
            f"{step:<4}",
            f"{T.to('degC'):>15~P.3f}",
            f"{P.to('kPa'):>15~P.1f}",
            f"{x.to('frac'):>15~P.3f}",
            f"{(1 / rho).to('m**3 / kg'):>15~P.6f}",
            f"{h.to('kJ / kg'):>15~P.3f}",
            f"{vel.to('m / s'):>15~P.3f}",
            f"{delta_L.to('m'):>15~P.4f}",
            f"{L.to('m'):>15~P.4f}",
            sep="\t", end="\n"
        )


if __name__ == '__main__':

    def main():
        R22 = Fluid("R22")
        capillary_tube = CapillaryTube(
            d_i=Q_(1.63, 'mm'),
            refrigerant=R22,
            m_dot=Q_(0.01, 'kg / s'),
            T_cnd=Q_(40, 'degC'),
            T_evp=Q_(-3, 'degC')
        )
        L = capillary_tube.get_length()
        print(f"required capillary tube length = {L.to('cm'):~P.3f}")

    main()
