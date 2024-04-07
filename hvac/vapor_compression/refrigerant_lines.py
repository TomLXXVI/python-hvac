from abc import ABC
from hvac import Quantity
from hvac.fluids import FluidState
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube


Q_ = Quantity
g = Q_(9.81, 'm / s ** 2')


class RefrigerantLineSection:

    def __init__(
        self,
        L_eq: Quantity,
        D_int: Quantity,
        elevation: Quantity = Q_(0, 'm'),
        e: Quantity = Q_(0.0015, 'mm')
    ) -> None:
        self.L_eq = L_eq.to('m')
        self.elevation = elevation.to('m')
        self.D_int = D_int.to('m')
        self.e = e.to('m')
        self._rfg_state: FluidState | None = None
        self._m_dot: Quantity | None = None

    @property
    def rfg_state(self) -> FluidState:
        return self._rfg_state

    @rfg_state.setter
    def rfg_state(self, v: FluidState) -> None:
        self._rfg_state = v

    @property
    def m_dot(self) -> Quantity:
        return self._m_dot

    @m_dot.setter
    def m_dot(self, v: Quantity) -> Quantity:
        self._m_dot = v.to('kg / s')

    @property
    def pressure_drop(self) -> Quantity:
        tube = CircularTube(
            Di=self.D_int,
            L=self.L_eq,
            fluid=self._rfg_state,
            e=Q_(0.0015, 'mm')
        )
        tube.m_dot = self._m_dot
        dP_friction = tube.pressure_drop()
        if self.elevation.magnitude > 0.0:
            dP_elevation = g * self._rfg_state.rho * self.elevation
        else:
            dP_elevation = Q_(0.0, 'Pa')
        dP = dP_friction + dP_elevation
        return dP.to('Pa')


class RefrigerantLine(ABC):

    def __init__(self):
        self.sections: list[RefrigerantLineSection] = []
        self._rfg_state: FluidState | None = None

    def add_section(
        self,
        L_eq: Quantity,
        D_int: Quantity,
        elevation: Quantity = Q_(0, 'm'),
        e: Quantity = Q_(0.0015, 'mm')
    ) -> None:
        """Adds a piping section to the refrigerant line.

        Parameters
        ----------
        L_eq:
            Total length of the section including the equivalent lengths of
            fittings and accessories in the section.
        D_int:
            Inside diameter of the piping section.
        elevation:
            Elevation of the section outlet with respect to the inlet.
        e:
            Copper tube wall roughness.

        Notes
        -----
        Pipe sections should be added in the order that they appear in the
        refrigerant line, following the direction of refrigerant flow.
        """
        section = RefrigerantLineSection(L_eq, D_int, elevation, e)
        self.sections.append(section)

    @property
    def rfg_state(self) -> FluidState:
        return self._rfg_state

    @rfg_state.setter
    def rfg_state(self, v: FluidState) -> None:
        """Sets the refrigerant state at the entrance of the refrigerant line.
        In the case of a liquid line, this is the state of the refrigerant at
        the condenser outlet. In the case of a discharge line, this is the state
        of the refrigerant at the compressor outlet. In the case of a suction
        line, this is the state of the refrigerant at the evaporator outlet.
        """
        Rfg = v.fluid
        n = len(self.sections)
        for i in range(n):
            if i == 0:
                self.sections[i].rfg_state = v
            if i < n - 1:
                dP = self.sections[i].pressure_drop
                P = self.sections[i].rfg_state.P - dP
                T = self.sections[i].rfg_state.T
                self.sections[i + 1].rfg_state = Rfg(T=T, P=P)

    @property
    def m_dot(self) -> Quantity:
        return self.sections[0].m_dot

    @m_dot.setter
    def m_dot(self, v: Quantity) -> Quantity:
        for section in self.sections:
            section.m_dot = v.to('kg / s')

    @property
    def pressure_drop(self) -> Quantity:
        dP = sum(section.pressure_drop for section in self.sections)
        return dP


class SuctionLine(RefrigerantLine):
    pass


class DischargeLine(RefrigerantLine):
    pass


class LiquidLine(RefrigerantLine):
    pass
