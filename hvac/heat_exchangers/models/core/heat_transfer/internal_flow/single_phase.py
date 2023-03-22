from typing import Optional
from abc import ABC, abstractmethod
import numpy as np
from hvac import Quantity
from hvac.fluids import FluidState, HumidAir
from hvac.heat_exchangers.models.core.geometry import CoilGeometry

Q_ = Quantity

# Friction factors are Moody or Darcy friction factors (= 4x Fanning friction factor).


class AbstractSinglePhaseFlow(ABC):

    def __init__(
        self,
        fluid: FluidState,
        v: Quantity,
        Dh: Quantity,
        L: Quantity,
        e: Optional[Quantity] = None,
    ) -> None:
        """
        Parameters
        ----------
        fluid: Fluid
            (Average) fluid condition.
        v: Quantity
            Flow velocity.
        Dh: Quantity
            Hydraulic diameter.
        L: Quantity
            Tube length.
        e: Quantity
            Tube wall roughness.
        """
        self.fluid = fluid
        self.v = v.to('m / s')
        self.Dh = Dh.to('m')
        self.e = e.to('m')
        self.L = L.to('m')

    @property
    def Re(self) -> Quantity:
        Re = self.fluid.rho * self.v * self.Dh / self.fluid.mu
        return Re.to('')

    @property
    def Pr(self) -> Quantity:
        Pr = self.fluid.cp * self.fluid.mu / self.fluid.k
        return Pr.to('')

    @property
    def f_fd(self) -> Quantity:
        """Friction factor for fully developed laminar/turbulent flow."""
        f_fd = self._f_fd()
        return f_fd

    @property
    def f_avg(self) -> Quantity:
        """Average friction factor for laminar/turbulent flow along tube length L."""
        f_avg = self._f_avg()
        return f_avg

    @property
    def dP(self) -> Quantity:
        """Pressure drop for tube length L."""
        dP = self.f_avg * (self.L / self.Dh) * (0.5 * self.fluid.rho * self.v ** 2)
        return dP

    @property
    def Nu_fd(self) -> Quantity:
        """Nusselt number for fully developed laminar/turbulent flow."""
        Nu_fd = self._Nu_fd()
        return Nu_fd

    @property
    def Nu_avg(self) -> Quantity:
        """Average Nusselt number for laminar/turbulent flow along tube length L."""
        Nu_avg = self._Nu_avg()
        return Nu_avg

    @property
    def h_fd(self) -> Quantity:
        """Heat transfer coefficient for fully developed laminar/turbulent flow."""
        h = self.Nu_fd * self.fluid.k / self.Dh
        return h.to('W / (m ** 2 * K)')

    @property
    def h_avg(self) -> Quantity:
        """Average heat transfer coefficient for laminar/turbulent flow along tube length L."""
        h_avg = self.Nu_avg * self.fluid.k / self.Dh
        return h_avg.to('W / (m ** 2 * K)')

    @abstractmethod
    def _f_fd(self) -> Quantity:
        ...

    @abstractmethod
    def _f_avg(self) -> Quantity:
        ...

    @abstractmethod
    def _Nu_fd(self) -> Quantity:
        ...

    @abstractmethod
    def _Nu_avg(self) -> Quantity:
        ...


class TurbulentFlow(AbstractSinglePhaseFlow):
    """Implements the abstract methods of class `AbstractSinglePhaseFlow` for turbulent flow."""

    def _f_fd_Offor_Alabi(self) -> Quantity:
        # friction factor for fully developed turbulent flow in a rough duct
        # with 4e3 < Re < 1e8 and 1e-6 < e/Dh < 5e-2
        # according to Offor and Alabi (2016)
        C1 = 3.71
        C2 = -1.975
        C3 = 3.93
        C4 = 1.092
        C5 = 7.627
        C6 = 395.9
        C7 = -2
        a = self.e / (C1 * self.Dh)
        b = C2 / self.Re
        c = (self.e / (C3 * self.Dh)) ** C4
        d = C5 / (self.Re + C6)
        f_fd = (-2.0 * np.log10(a + b * np.log(c + d))) ** C7
        return f_fd

    def _f_fd_Li(self) -> Quantity:
        # friction factor for fully developed turbulent flow in a smooth duct
        # with 4e3 < Re < 1e7 and e/Dh -> 0
        # according to Li et al. (2011)
        C1 = -0.0015702
        C2 = 0.3942031
        C3 = 2.5341533
        a = C1 / np.log(self.Re)
        b = C2 / np.log(self.Re) ** 2
        c = C3 / np.log(self.Re) ** 3
        f_fd = 4 * (a + b + c)
        return f_fd

    def _f_fd(self) -> Quantity:
        if 1e-6 < self.e / self.Dh < 5e-2:
            # rough duct
            f_fd = self._f_fd_Offor_Alabi()
        else:
            # smooth duct
            f_fd = self._f_fd_Li()
        return f_fd

    def _f_avg(self) -> Quantity:
        f_avg = self.f_fd * (1 + (self.Dh / self.L) ** 0.7)
        return f_avg

    def _Nu_fd_Gnielinski(self) -> Quantity:
        # Nusselt number for fully developed turbulent flow with
        # 0.5 < Pr < 2000 and 2300 < Re < 5e6
        # according to Gnielinski (1976)
        n = (self.f_fd / 8) * (self.Re - 1000) * self.Pr
        d = 1 + 12.7 * (self.Pr ** (2 / 3) - 1) * np.sqrt(self.f_fd / 8)
        Nu_fd = n / d
        return Nu_fd

    def _Nu_fd_Notter_Sleicher(self) -> Quantity:
        # Nusselt number for fully developed turbulent flow with
        # 0.004 < Pr < 0.1 and 1e4 < Re < 1e6
        # according to Notter and Sleicher (1972).
        Nu_fd_T = 4.8 + 0.0156 * (self.Re ** 0.85) * (self.Pr ** 0.93)
        Nu_fd_H = 6.3 + 0.0167 * (self.Re ** 0.85) * (self.Pr ** 0.93)
        Nu_fd = (Nu_fd_T + Nu_fd_H) / 2
        return Nu_fd

    def _Nu_fd(self) -> Quantity:
        if (0.004 < self.Pr < 0.1) and (1e4 < self.Re < 1e6):
            return self._Nu_fd_Notter_Sleicher()
        elif (0.5 < self.Pr < 2000) and (2300 < self.Re < 5e8):
            return self._Nu_fd_Gnielinski()

    def _Nu_avg(self) -> Quantity:
        # Average Nusselt number for turbulent flow along tube length L (according to Kakaç et al. 1987)
        C = 1.0
        m = 0.7
        Nu_avg = self.Nu_fd * (1 + C * (self.L / self.Dh) ** -m)
        return Nu_avg


class AbstractLaminarFlow(AbstractSinglePhaseFlow, ABC):
    """Implements common properties of laminar flow, regardless of conduit shape."""

    @property
    def L_tilde(self) -> Quantity:
        # dimensionless tube length
        # the local friction factor will approach the fully developed value when L_tilde > 0.063
        L_dl = self.L / (self.Dh * self.Re)
        return L_dl

    @property
    def Gz(self) -> Quantity:
        # Graetz number
        Gz = self.Pr / self.L_tilde
        return Gz


class LaminarFlowCircularConduit(AbstractLaminarFlow):
    """Implements abstract methods of class `AbstractSinglePhaseFlow` for laminar flow in circular conduits."""

    def _f_fd(self) -> Quantity:
        f_fd = 64.0 / self.Re
        return f_fd

    def _f_avg(self) -> Quantity:
        t = 0.215 * (self.L_tilde ** -0.5)
        n = 0.0196 * (self.L_tilde ** -1) + 1 - 0.215 * (self.L_tilde ** -0.5)
        d = 1 + 0.00021 * (self.L_tilde ** -2)
        f_avg_to_f_fd = t + n / d
        f_avg = f_avg_to_f_fd * self.f_fd
        return f_avg

    @staticmethod
    def _Nu_fd_h() -> Quantity:
        # local Nusselt number for laminar, hydrodynamically and thermally fully developed flow
        # in a circular tube and uniform heat flux
        return 4.36

    @staticmethod
    def _Nu_fd_t() -> Quantity:
        # local Nusselt number for laminar, hydrodynamically and thermally fully developed flow
        # in a circular tube and uniform wall temperature
        return 3.66

    def _Nu_fd(self) -> Quantity:
        Nu_fd_t = self._Nu_fd_t()
        Nu_fd_h = self._Nu_fd_h()
        Nu_fd = (Nu_fd_t + Nu_fd_h) / 2
        return Nu_fd

    def _Nu_avg_t(self) -> Quantity:
        # average Nusselt number for simultaneously developing flow with a constant wall temperature
        t = self._Nu_fd_t()
        n = (0.049 + 0.020 / self.Pr) * (self.Gz ** 1.12)
        d = 1 + 0.065 * (self.Gz ** 0.7)
        Nu_avg_t = t + n / d
        return Nu_avg_t

    def _Nu_avg_h(self) -> Quantity:
        # average Nusselt number for simultaneously developing flow with constant heat flux
        t = self._Nu_fd_h()
        n = (0.1156 + 0.08569 / (self.Pr ** 0.4)) * self.Gz
        d = 1 + 0.1158 * (self.Gz ** 0.6)
        Nu_avg_h = t + n / d
        return Nu_avg_h

    def _Nu_avg(self) -> Quantity:
        Nu_avg_t = self._Nu_avg_t()
        Nu_avg_h = self._Nu_avg_h()
        Nu_avg = (Nu_avg_t + Nu_avg_h) / 2
        return Nu_avg


class SinglePhaseFlow(AbstractSinglePhaseFlow):
    """
    Implements the abstract methods of class `AbstractSinglePhaseFlow` taking the kind of flow (turbulent or
    laminar) into account.
    """
    Re_crit: int = 2300

    def __init__(
        self,
        fluid: FluidState,
        v: Quantity,
        Dh: Quantity,
        L: Quantity,
        e: Quantity,
        shape_of_duct: Optional[str] = None
    ) -> None:
        super().__init__(fluid, v, Dh, L, e)
        self._turbulent_flow = TurbulentFlow(fluid, v, Dh, L, e)
        if (shape_of_duct is None) or (shape_of_duct == 'circular'):
            self._laminar_flow = LaminarFlowCircularConduit(fluid, v, Dh, L, e)
        # may be extended for other duct shapes (rectangular, annular)

    def _f_fd(self) -> float:
        if self.Re > self.Re_crit:
            return self._turbulent_flow.f_fd
        else:
            return self._laminar_flow.f_fd

    def _f_avg(self) -> float:
        if self.Re > self.Re_crit:
            return self._turbulent_flow.f_avg
        else:
            return self._laminar_flow.f_avg

    def _Nu_fd(self) -> float:
        if self.Re > self.Re_crit:
            return self._turbulent_flow.Nu_fd
        else:
            return self._laminar_flow.Nu_fd

    def _Nu_avg(self) -> float:
        if self.Re > self.Re_crit:
            return self._turbulent_flow.Nu_avg
        else:
            return self._laminar_flow.Nu_avg


class SinglePhaseHeatTransfer:

    def __init__(
        self,
        coolant: FluidState,
        coil_geometry: CoilGeometry,
        mc_tube: Quantity,
        ma: Quantity,
        vfa: Quantity,
        air: HumidAir,
        e: Optional[Quantity] = None
    ) -> None:
        """
        Parameters
        ----------
        coolant: Fluid
            Coolant condition the thermodynamic properties will be based upon.
        coil_geometry: CoilGeometry
            The geometrical data about the coil.
        mc_tube: Quantity
            Mass flow rate of coolant inside tubes.
        ma: Quantity
            Mass flow rate of air.
        vfa: Quantity
            Air face velocity.
        air: HumidAir
            Entering air condition.
        e: Quantity
            Tube wall roughness.
        """
        Afa = ma / (vfa * air.rho)              # face area
        Ao = coil_geometry.Ao_to_Afa * Afa      # external surface area of 1 row
        Ai = (1 / coil_geometry.Ao_to_Ai) * Ao  # internal surface area of 1 row
        Dh = coil_geometry.Di                   # inner tube diameter
        A_tube = np.pi * Dh ** 2 / 4            # internal cross-section of 1 tube
        vc = mc_tube / (coolant.rho * A_tube)   # coolant velocity in tubes
        L = (Ai / (np.pi * Dh))                 # length of all tubes in 1 row
        e = Q_(0, 'mm') if e is None else e
        self._single_phase_flow = SinglePhaseFlow(coolant, vc, Dh, L, e, None)

    @property
    def h(self) -> Quantity:
        h = self._single_phase_flow.h_avg
        return h
