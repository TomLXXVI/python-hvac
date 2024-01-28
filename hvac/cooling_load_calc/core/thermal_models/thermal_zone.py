"""NODAL THERMAL ZONE MODEL
---------------------------
Lumped capacitance model of a thermal zone or space with a variable zone air
temperature. This model is incorporated in class `UnconditionedZone`: see
`hvac.building.unconditioned_zone.py`.

In this type of model the zone air temperature is considered to be variable.
Its value follows from an energy balance of the zone air node, i.e., the
sum of the heat gains minus the heat rate extracted by the cooling system.

Solving the model means to determine at each time moment of the selected day
the zone air temperature, the heat gains or losses, and the heat rate extracted
by the cooling system or supplied by the heating system.

The heat transfer through exterior building elements is modeled by linear
thermal networks. The interior thermal mass of the zone is represented in the
model by a single temperature node.
"""
from __future__ import annotations
from typing import Callable, TYPE_CHECKING
from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from hvac import Quantity

if TYPE_CHECKING:
    from hvac.cooling_load_calc.core.construction_assembly import (
        ConstructionAssembly,
        ConstructionLayer
    )
    from hvac.cooling_load_calc.core.fenestration import Window
    from hvac.cooling_load_calc.core.building_element import (
        ExteriorBuildingElement,
        InteriorBuildingElement
    )
    from hvac.cooling_load_calc.building.conditioned_zone import SpaceVentilation


Q_ = Quantity
dt_sec: int = 3600  # time step in seconds


def set_time_step(dt_hr: float) -> float:
    """Sets the time step for solving the system of node equations in a linear
    thermal network.

    Parameters
    ----------
    dt_hr:
        Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4 means
        the time step for the calculations is one quarter of an hour.

    Returns
    -------
    The time step converted to seconds.
    """
    global dt_sec
    dt_sec = dt_hr * 3600
    return dt_sec


class TemperatureNode(ABC):
    """Represents a general temperature node.

    A temperature node has a unit thermal capacity `C`, expressed per unit of
    area, i.e., in SI units J/(K.m²).
    The node can be connected on two sides (the left and the right side) to
    either another temperature node in the linear network, or to the adjacent
    environment through a thermal resistor. The thermal resistor on the left is
    designated `R1`, while the resistor on the right is designated `R2`. The
    thermal resistance is expressed on the basis of one unit area, i.e., in SI
    units K.m²/W.
    The surface area, `A`, associated with a temperature node is considered
    separately.
    """
    def __init__(self):
        self.ID: str = ''
        self._C: float = 0.0                             # J / K
        self._A: float | list[float] = 1.0               # m**2
        self._R1: float | list[float] = float('inf')     # K / W
        self._R2: float = float('inf')                   # K / W
        self._F_rad: float = 0.46                        # -
        self._T_sa: Callable[[float], Quantity] = None

    @classmethod
    @abstractmethod
    def create(cls, *args, **kwargs) -> TemperatureNode:
        """Creates a `TemperatureNode` object."""
        ...

    @abstractmethod
    def get_a_coefficients(self) -> list[float]:
        """Returns the coefficients of the node equation to put into the
        coefficient matrix A of the linear thermal network to be solved.
        """
        ...

    @abstractmethod
    def get_b_value(self, *args, **kwargs) -> float:
        """Returns the input value of the node equation to put in the input
        matrix B of the linear thermal network to be solved.
        """
        ...

    @property
    def A(self) -> Quantity:
        """Returns the surface area associated with the temperature node.
        In case of a zone air node or thermal storage node a `Quantity` array is
        returned with the surface areas of the exterior building elements
        surrounding the zone and the surface area associated with the interior
        thermal mass.
        """
        return Q_(self._A, 'm**2')

    @property
    def R1(self) -> Quantity:
        """Returns the unit thermal resistance R1.
        In case of a thermal storage node or a zone air node a `Quantity` array
        is returned with the unit thermal resistance between the thermal storage
        node/zone air node and each exterior building element surrounding the
        zone.
        """
        if isinstance(self._R1, list):
            if isinstance(self._A, list):
                return Q_([r1 * a for r1, a in zip(self._R1, self._A[:-1])], 'K * m**2 / W')
            else:
                return Q_([r1 for r1 in self._R1] * self._A, 'K * m**2 / W')
        else:
            return Q_(self._R1 * self._A, 'K * m**2 / W')

    @property
    def R2(self) -> Quantity:
        """Returns the unit thermal resistance R2."""
        if isinstance(self._A, list):
            # in case of a zone air node or thermal storage node: multiply with
            # the last element in `self._A` which is the surface area associated
            # with the thermal storage node which is connected to the zone air
            # node by R2.
            return Q_(self._R2 * self._A[-1], 'K * m**2 / W')
        else:
            return Q_(self._R2 * self._A, 'K * m**2 / W')

    @property
    def C(self) -> Quantity:
        """Returns the thermal capacity per unit area C of the node."""
        if isinstance(self._A, list):
            # in case of a thermal storage node or zone air node:
            return Q_(self._C / self._A[-1], 'J / (K * m**2)')
        else:
            return Q_(self._C / self._A, 'J / (K * m**2)')

    @property
    def F_rad(self) -> Quantity:
        """Returns the radiative fraction."""
        return Q_(self._F_rad, 'frac')


class ExteriorSurfaceNode(TemperatureNode):
    """Represents the first temperature node at the exterior side of an
    exterior building element. It is connected to the exterior temperature
    (sol-air temperature).
    """
    @classmethod
    def create(
        cls,
        ID: str,
        C: Quantity = Q_(0.0, 'J / (K * m**2)'),
        A: Quantity = Q_(1.0, 'm**2'),
        R1: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        R2: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        T_sa: Callable[[float], Quantity] | None = None
    ) -> TemperatureNode:
        """Creates an `ExteriorSurfaceNode` object.

        Parameters
        ----------
        ID:
            Name to identify the node in a linear thermal network.
        C:
            The thermal capacity of the node per unit area.
        A:
            The surface area associated with the node.
        R1:
            The unit thermal resistance between the preceding node and this node.
        R2:
            The unit thermal resistance between this node and the next node.
        T_sa:
            Function that returns the sol-air temperature at a given solar time
            moment in seconds from midnight (0 s).
        """
        node = cls()
        node.ID = ID
        node._A = A.to('m**2').magnitude
        node._C = C.to('J / (K * m**2)').magnitude * node._A
        node._R1 = R1.to('K * m**2 / W').magnitude / node._A
        node._R2 = R2.to('K * m**2 / W').magnitude / node._A
        node._T_sa = T_sa
        return node
    
    @property
    def _a1(self) -> float:
        """Coefficient in the node equation of the exterior surface node
        temperature.
        """
        a1 = -(3 + (2 * dt_sec / self._C) * (1 / self._R1 + 1 / self._R2))
        return a1

    @property
    def _a2(self) -> float:
        """Coefficient in the node equation of the adjacent building mass node
        temperature.
        """
        a2 = (2 * dt_sec / self._C) * (1 / self._R2)
        return a2

    def _b(self, T: list[float], t_sol_sec: float) -> float:
        """Calculates the input-side (rhs) of the node equation at the current
        time moment t or time index k (t = k * dt).

        Parameters
        ----------
        T:
            List with the previous node temperature at time index k-2 and
            time index k-1.
        t_sol_sec:
            The current solar time in seconds from midnight (0 s).
        """
        k = (2 * dt_sec / self._C) * (1 / self._R1)
        b = T[0] - 4 * T[1] - k * self._T_sa(t_sol_sec).to('K').m
        return b

    def get_a_coefficients(self) -> list[float]:
        return [self._a1, self._a2]

    def get_b_value(self, T: list[float], t_sol_sec: float) -> float:
        return self._b(T, t_sol_sec)


class BuildingMassNode(TemperatureNode):
    """Represents a temperature node inside the exterior building element.
    It is connected on the left side with an exterior surface node or another
    building mass node. On the right side it is connected with an interior
    surface node or another building mass node.
    """

    @classmethod
    def create(
        cls,
        ID: str,
        C: Quantity = Q_(0.0, 'J / (K * m**2)'),
        A: Quantity = Q_(1.0, 'm**2'),
        R1: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        R2: Quantity = Q_(float('inf'), 'K * m**2 / W')
    ) -> TemperatureNode:
        """Creates a `BuildingMassNode` object.

        Parameters
        ----------
        ID:
            Name to identify the node in a linear thermal network.
        C:
            The thermal capacity of the node per unit area.
        A:
            The surface area associated with the node.
        R1:
            The unit thermal resistance between the preceding node and this node.
        R2:
            The unit thermal resistance between this node and the next node.
        """
        node = cls()
        node.ID = ID
        node._A = A.to('m**2').magnitude
        node._C = C.to('J / (K * m**2)').magnitude * node._A
        node._R1 = R1.to('K * m**2 / W').magnitude / node._A
        node._R2 = R2.to('K * m**2 / W').magnitude / node._A
        return node
    
    @property
    def _a1(self) -> float:
        """Coefficient in the node equation of the preceding building mass node
        temperature or exterior surface node temperature.
        """
        a1 = 2 * dt_sec / (self._C * self._R1)
        return a1

    @property
    def _a2(self) -> float:
        """Coefficient in the node equation of this building mass node
        temperature.
        """
        a2 = -(3 + 2 * (dt_sec / self._C) * (1 / self._R1 + 1 / self._R2))
        return a2

    @property
    def _a3(self) -> float:
        """Coefficient in the node equation of the next building mass node
        temperature or interior surface node temperature.
        """
        a3 = 2 * dt_sec / (self._C * self._R2)
        return a3

    @staticmethod
    def _b(T: list[float]) -> float:
        """Calculates the input-side (rhs) of the node equation at the current
        time moment t or time index k (t = k * dt).

        Parameters
        ----------
        T:
            List with the previous node temperature at time index k-2 and
            time index k-1.
        """
        b = T[0] - 4 * T[1]
        return b

    def get_a_coefficients(self) -> list[float]:
        return [self._a1, self._a2, self._a3]

    def get_b_value(self, T: list[float]) -> float:
        return self._b(T)


class InteriorSurfaceNode(TemperatureNode):
    """Represents a temperature node on the interior surface of an exterior
    building element. It has no thermal capacity. On the right side the
    node is connected to the zone air temperature. It is also connected
    to the thermal storage node of the zone.
    As the interior surface node has no thermal capacity, the heat flow into
    the interior surface node equals the heat flow out to the zone air. The
    amount of convective heat transfer to the zone air node and the amount of
    radiative heat transfer to the thermal storage node depend on the value
    that was set for the radiative fraction `F_rad`.
    """
    @classmethod
    def create(
        cls,
        ID: str,
        A: Quantity = Q_(1.0, 'm**2'),
        R1: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        R2: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        F_rad: Quantity = Q_(0.46, 'frac')
    ) -> TemperatureNode:
        """Creates an `InteriorSurfaceNode` object.

        Parameters
        ----------
        ID:
            Name to identify the node in a linear thermal network.
        A:
            The surface area associated with the node.
        R1:
            The unit thermal resistance between the preceding node and this node.
        R2:
            The unit thermal resistance between this node and the next node.
        F_rad:
            The radiative fraction of the conductive heat gain through the
            exterior building elements.
        """
        node = cls()
        node.ID = ID
        node._A = A.to('m**2').magnitude
        node._C = 0.0
        node._R1 = R1.to('K * m**2 / W').magnitude / node._A
        node._R2 = R2.to('K * m**2 / W').magnitude / node._A
        node._F_rad = F_rad.to('frac').magnitude
        return node
    
    @property
    def _a1(self) -> float:
        """Coefficient in the node equation of the preceding building mass node
        temperature.
        """
        a1 = 1 / self._R1
        return a1

    @property
    def _a2(self) -> float:
        """Coefficient in the node equation of this internal surface node
        temperature.
        """
        a2 = -((1 / self._R1) + (1 / ((1 - self._F_rad) * self._R2)))
        return a2

    @property
    def _a3(self) -> float:
        """Coefficient in the node equation of the zone air node
        temperature.
        """
        a3 = 1 / ((1 - self._F_rad) * self._R2)
        return a3

    def get_a_coefficients(self) -> list[float]:
        return [self._a1, self._a2, self._a3]

    def get_b_value(self) -> float:
        return 0.0


class ThermalStorageNode(TemperatureNode):
    """Represents the interior thermal mass in a zone.
    Each exterior building element surrounding the zone is connected with its
    interior surface node to the thermal storage node of the zone.
    """
    def __init__(self):
        super().__init__()
        self.windows: list[Window] | None = None
        self.ext_doors: list[ExteriorBuildingElement] | None = None
        self.int_build_elems: list[InteriorBuildingElement] | None = None
        self.Q_dot_sol_rd: Callable[[float], Quantity] | None = None
        self.Q_dot_ihg_rd: Callable[[float], Quantity] | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        C: Quantity = Q_(0.0, 'J / (K * m**2)'),
        A: list[Quantity] = Q_(1.0, 'm**2'),
        R1: list[Quantity] = Q_(float('inf'), 'K * m**2 / W'),
        R2: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        F_rad: Quantity = Q_(0.46, 'frac'),
        windows: list[Window] | None = None,
        ext_doors: list[ExteriorBuildingElement] | None = None,
        int_build_elems: list[InteriorBuildingElement] | None = None,
        Q_dot_sol_rd: Callable[[float], Quantity] | None = None,
        Q_dot_ihg_rd: Callable[[float], Quantity] | None = None
    ) -> TemperatureNode:
        """Creates a `ThermalStorageNode` object.

        Parameters
        ----------
        ID:
            Name to identify the node in a linear thermal network.
        C:
            The thermal capacity of the node per unit area.
        A:
            List with the surface areas of the exterior building elements and
            the last element being the surface area of the interior thermal mass.
        R1:
            List with the unit thermal resistances between each interior surface
            node and this thermal storage node.
        R2:
            The unit thermal resistance between this thermal storage node and
            the zone air node.
        F_rad:
            The radiative fraction of the conductive heat gain through the
            exterior building elements.
        windows:
            List of the windows that are present in the exterior envelope of
            the zone.
        ext_doors:
            List of the exterior doors that are present in the exterior envelope
            of the zone.
        int_build_elems:
            List of the interior building elements, including interior doors,
            that surround the zone.
        Q_dot_sol_rd:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the radiative part of the solar heat
            gain.
        Q_dot_ihg_rd:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the radiative part of the internal heat
            gains in the zone.
        """
        node = cls()
        node.ID = ID
        node._A = [a.to('m**2').magnitude for a in A]
        node._C = C.to('J / (K * m**2)').magnitude * node._A[-1]
        node._R1 = [
            r1.to('K * m**2 / W').magnitude / a
            for r1, a in zip(R1, node._A[:-1])
        ]
        node._R2 = R2.to('K * m**2 / W').magnitude / node._A[-1]
        node._F_rad = F_rad.to('frac').magnitude
        node.windows = windows
        node.ext_doors = ext_doors
        node.int_build_elems = int_build_elems
        node.Q_dot_sol_rd = Q_dot_sol_rd
        node.Q_dot_ihg_rd = Q_dot_ihg_rd
        return node
    
    @property
    def _a1(self) -> list[float]:
        """Coefficients in the node equation of all the interior surface node
        temperatures connected to the thermal storage node temperature.
        """
        a1 = [
            (2 * dt_sec * self._F_rad) / ((1 - self._F_rad) * r1 * self._C)
            for r1 in self._R1
        ]
        return a1

    @property
    def _a2(self) -> float:
        """Coefficient in the node equation of this thermal storage node
        temperature.
        """
        a2 = -(((2 * dt_sec) / (self._R2 * self._C)) + 3)
        return a2

    @property
    def _a3(self) -> float:
        """Coefficient in the node equation of the zone air node temperature."""
        if self.windows:
            UA_wnd = sum(
                wnd.F_rad * wnd.UA.to('W / K').m
                for wnd in self.windows
            )
        else:
            UA_wnd = 0.0
        if self.ext_doors:
            UA_edr = sum(
                edr.F_rad * edr.UA.to('W / K').m
                for edr in self.ext_doors
            )
        else:
            UA_edr = 0.0
        if self.int_build_elems:
            UA_ibe = sum(
                ibe.F_rad * ibe.UA.to('W / K').m
                for ibe in self.int_build_elems
            )
        else:
            UA_ibe = 0.0
        UA = UA_wnd + UA_edr + UA_ibe
        a3 = sum(
            (2 * dt_sec * self._F_rad) / (self._C * (1 - self._F_rad) * r1)
            for r1 in self._R1
        )
        a3 += (2 * dt_sec / self._C) * UA
        a3 = -a3 + ((2 * dt_sec) / (self._C * self._R2))
        return a3

    def _b(self, T: list[float], t_sol_sec: float) -> float:
        """Calculates the input-side (rhs) of the node equation at the current
        solar time moment t_sol_sec in seconds from midnight (0 s).

        Parameters
        ----------
        T:
            List with the previous node temperature at time index k-2 and
            time index k-1.
        t_sol_sec:
            The current solar time in seconds from midnight (0 s).
        """
        Q_dot_wnd, Q_dot_edr, Q_dot_ibe = 0.0, 0.0, 0.0
        if self.windows:
            Q_dot_wnd = sum(
                (wnd.F_rad * wnd.UA * wnd._ext_surf.T_db(t_sol_sec).to('K')).to('W').m
                for wnd in self.windows
            )
        if self.ext_doors:
            Q_dot_edr = sum(
                (edr.F_rad * edr.UA * edr._ext_surf.T_sa(t_sol_sec).to('K')).to('W').m
                for edr in self.ext_doors
            )
        if self.int_build_elems:
            Q_dot_ibe = sum(
                (ibe.F_rad * ibe.UA * ibe.T_adj(t_sol_sec).to('K')).to('W').m
                for ibe in self.int_build_elems
            )
        Q_dot_rad = (
            Q_dot_wnd + Q_dot_edr + Q_dot_ibe
            + self.Q_dot_sol_rd(t_sol_sec).to('W').m
            + self.Q_dot_ihg_rd(t_sol_sec).to('W').m
        )
        b = T[0] - 4 * T[1] - (2 * dt_sec / self._C) * Q_dot_rad
        return b

    def get_a_coefficients(self) -> list[float]:
        a = self._a1
        a.append(self._a2)
        a.append(self._a3)
        return a

    def get_b_value(self, T: list[float], t_sol_sec: float) -> float:
        return self._b(T, t_sol_sec)


class ZoneAirNode(TemperatureNode):

    def __init__(self):
        super().__init__()
        self.windows: list[Window] | None = None
        self.ext_doors: list[ExteriorBuildingElement] | None = None
        self.int_build_elems: list[InteriorBuildingElement] | None = None
        self.ventilation: SpaceVentilation | None = None
        self.Q_dot_sol_cv: Callable[[float], Quantity] | None = None
        self.Q_dot_ihg_cv: Callable[[float], Quantity] | None = None
        self.Q_dot_sys: Callable[[float, ...], Quantity] | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        C: Quantity = Q_(0.0, 'J / (K * m**2)'),
        A: list[Quantity] = Q_(1.0, 'm**2'),
        R1: list[Quantity] = Q_(float('inf'), 'K * m**2 / W'),
        R2: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        F_rad: Quantity = Q_(0.46, 'frac'),
        windows: list[Window] | None = None,
        ext_doors: list[ExteriorBuildingElement] | None = None,
        int_build_elems: list[InteriorBuildingElement] | None = None,
        ventilation: SpaceVentilation | None = None,
        Q_dot_sol_cv: Callable[[float], Quantity] | None = None,
        Q_dot_ihg_cv: Callable[[float], Quantity] | None = None,
        Q_dot_sys: Callable[[float, ...], Quantity] | None = None
    ) -> TemperatureNode:
        """Creates a `ZoneAirNode` object.

        Parameters
        ----------
        ID:
            Name to identify the node in a linear thermal network.
        C:
            The thermal capacity of the node per unit area.
        A:
            List with the surface areas of the exterior building elements,
            the second-to-last element being the surface area associated with
            the interior thermal mass, and the last element the floor area of
            the zone (used to determine the total thermal capacity of the zone
            air).
        R1:
            List with the unit thermal resistances between each interior surface
            node and this zone air node.
        R2:
            The unit thermal resistance between this zone air node and the
            thermal storage node.
        F_rad:
            The radiative fraction of the conductive heat gain through the
            exterior building elements.
        windows:
            List of the windows that are present in the exterior envelope of
            the zone.
        ext_doors:
            List of the exterior doors that are present in the exterior envelope
            of the zone.
        int_build_elems:
            List of the interior building elements, including interior doors,
            that surround the zone.
        ventilation:
            `SpaceVentilation` object that encapsulates the
            ventilation/infiltration model of the zone.
        Q_dot_sol_cv:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the convective part of the solar heat
            gain.
        Q_dot_ihg_cv:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the convective part of the internal heat
            gains in the zone.
        Q_dot_sys:
            A function with signature
                f(t_sol_sec: float, T_zone: float) -> Quantity
            which takes the time `t_sol_sec` in seconds from midnight (0 s) of
            the considered day and the zone air temperature in Kelvins (K). It
            returns the cooling capacity (`Quantity` object) of the cooling
            system (i.e. the heat rate extracted from the zone air by the
            cooling system). Note that heat extraction by the cooling system
            is associated with a positive sign. Should instead the system supply
            heat to the zone, this must be associated with a negative sign.
            If this parameter is left to None, it is assumed that no (operating)
            cooling system is present in the zone, i.e. the zone is truly
            unconditioned.
        """
        node = cls()
        node.ID = ID
        node._A = [a.to('m**2').magnitude for a in A]
        node._C = C.to('J / (K * m**2)').magnitude * node._A[-1]
        node._R1 = [
            r1.to('K * m**2 / W').m / a
            for r1, a in zip(R1, node._A[:-2])
        ]
        node._R2 = R2.to('K * m**2 / W').magnitude / node._A[-2]
        node._F_rad = F_rad.to('frac').magnitude
        node.windows = windows
        node.ext_doors = ext_doors
        node.int_build_elems = int_build_elems
        node.ventilation = ventilation
        node.Q_dot_sol_cv = Q_dot_sol_cv
        node.Q_dot_ihg_cv = Q_dot_ihg_cv
        node.Q_dot_sys = Q_dot_sys
        return node
    
    @property
    def _a1(self) -> list[float]:
        """Coefficients in the node equation of the interior surface node
        temperatures connected to this zone air node temperature.
        """
        a1 = [2 * dt_sec / (r1 * self._C) for r1 in self._R1]
        return a1

    @property
    def _a2(self) -> float:
        """Coefficient in the node equation of the thermal storage node
        temperature.
        """
        a2 = 2 * dt_sec / (self._R2 * self._C)
        return a2

    @property
    def _a3(self) -> float:
        """Coefficient in the node equation of this zone air node temperature."""
        UA_wnd, UA_edr, UA_ibe, UA_vent = 0.0, 0.0, 0.0, 0.0
        if self.windows:
            UA_wnd = sum(
                (1 - wnd.F_rad) * wnd.UA.to('W / K').m
                for wnd in self.windows
            )
        if self.ext_doors:
            UA_edr = sum(
                (1 - edr.F_rad) * edr.UA.to('W / K').m
                for edr in self.ext_doors
            )
        if self.int_build_elems:
            UA_ibe = sum(
                (1 - ibe.F_rad) * ibe.UA.to('W / K').m
                for ibe in self.int_build_elems
            )
        if self.ventilation:
            V_dot = (
                self.ventilation.V_dot_ext
                + self.ventilation.V_dot_sup
                + self.ventilation.V_dot_trf
            )
            UA_vent = 0.34 * V_dot
        UA = UA_wnd + UA_edr + UA_ibe + UA_vent
        a3 = sum(2 * dt_sec / (r1 * self._C) for r1 in self._R1)
        a3 += 2 * dt_sec / (self._R2 * self._C)
        a3 += 2 * dt_sec * UA / self._C
        a3 += 3
        a3 = -a3
        return a3

    def _b(self, T: list[float], t_sol_sec: float) -> float:
        """Calculates the input-side (rhs) of the node equation at the current
        solar time moment `t_sol_sec` in seconds from midnight (0 s).
        """
        Q_dot_wnd, Q_dot_edr, Q_dot_ibe, Q_dot_vent = 0.0, 0.0, 0.0, 0.0
        if self.windows:
            Q_dot_wnd = sum(
                ((1 - wnd.F_rad) * wnd.UA * wnd._ext_surf.T_db(t_sol_sec).to('K')).to('W').m
                for wnd in self.windows
            )
        if self.ext_doors:
            Q_dot_edr = sum(
                ((1 - edr.F_rad) * edr.UA * edr._ext_surf.T_sa(t_sol_sec).to('K')).to('W').m
                for edr in self.ext_doors
            )
        if self.int_build_elems:
            Q_dot_ibe = sum(
                ((1 - ibe.F_rad) * ibe.UA * ibe.T_adj(t_sol_sec).to('K')).to('W').m
                for ibe in self.int_build_elems
            )
        if self.ventilation:
            Q_dot_vent_ext = (
                0.34 * self.ventilation.V_dot_ext
                * self.ventilation.T_db_ext(t_sol_sec).to('K').m
            )
            Q_dot_vent_sup = (
                0.34 * self.ventilation.V_dot_sup
                * self.ventilation.T_sup(t_sol_sec).to('K').m
            )
            if self.ventilation.T_trf:
                Q_dot_vent_trf = (
                    0.34 * self.ventilation.V_dot_trf
                    * self.ventilation.T_trf(t_sol_sec).to('K').m
                )
            else:
                Q_dot_vent_trf = 0.0
            Q_dot_vent = Q_dot_vent_ext + Q_dot_vent_sup + Q_dot_vent_trf
        Q_dot_conv = (
            Q_dot_wnd + Q_dot_edr + Q_dot_ibe + Q_dot_vent
            + self.Q_dot_sol_cv(t_sol_sec).to('W').m
            + self.Q_dot_ihg_cv(t_sol_sec).to('W').m
        )
        r = (
            T[0] - 4 * T[1]
            + (2 * dt_sec / self._C)
            * (self.Q_dot_sys(t_sol_sec, T[1]).to('W').m - Q_dot_conv)
        )
        return r

    def get_a_coefficients(self) -> list[float]:
        a = self._a1
        a.append(self._a2)
        a.append(self._a3)
        return a

    def get_b_value(self, T: list[float], t_sol_sec) -> float:
        return self._b(T, t_sol_sec)


class NodalThermalZoneModelBuilder:
    """Creates the nodes of a nodal thermal zone model given the exterior
    building elements surrounding the zone, the specs of the zone's interior
    thermal mass, and the time-functions of radiative and convective heat gain
    in the zone (solar heat gain, internal heat gain, ventilation/infiltration
    heat gain).
    """
    def __init__(
        self,
        ebe_lst: list[ExteriorBuildingElement],
        F_rad: Quantity = Q_(0.46, 'frac'),
        A_tsn: Quantity = Q_(1.0, 'm**2'),
        C_tsn: Quantity = Q_(0.0, 'J / (K * m**2)'),
        R_tsn: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        A_zan: Quantity = Q_(1.0, 'm**2'),
        C_zan: Quantity = Q_(0.0, 'J / (K * m**2)'),
        windows: list[Window] | None = None,
        ext_doors: list[ExteriorBuildingElement] | None = None,
        int_build_elems: list[InteriorBuildingElement] | None = None,
        ventilation: SpaceVentilation | None = None,
        Q_dot_sol_cv: Callable[[float], Quantity] | None = None,
        Q_dot_sol_rd: Callable[[float], Quantity] | None = None,
        Q_dot_ihg_cv: Callable[[float], Quantity] | None = None,
        Q_dot_ihg_rd: Callable[[float], Quantity] | None = None,
        Q_dot_sys: Callable[[float, ...], Quantity] | None = None
    ) -> None:
        """Creates a `NodalThermalZoneModelBuilder` object.

        Parameters
        ----------
        ebe_lst:
            List of the opaque exterior building elements that surround the
            zone.
        F_rad:
            The radiative factor of the conduction heat gain trough the opaque
            exterior building elements.
        A_tsn:
            The surface area associated with the interior thermal mass of the
            zone.
        C_tsn:
            The effective unit thermal capacity of the interior thermal mass.
        R_tsn:
            The unit thermal resistance between the interior thermal mass and
            the zone air.
        A_zan:
            The floor area of the zone.
        C_zan:
            The thermal capacity of the zone air per unit floor area.
        windows:
            List of the windows that are present in the exterior envelope of
            the zone.
        ext_doors:
            List of the exterior doors that are present in the exterior envelope
            of the zone.
        int_build_elems:
            List of the interior building elements, including interior doors,
            that surround the zone.
        ventilation:
            `SpaceVentilation` object that encapsulates the
            ventilation/infiltration model of the zone.
        Q_dot_sol_cv:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the convective part of the solar heat
            gain.
        Q_dot_sol_rd:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the radiative part of the solar heat
            gain.
        Q_dot_ihg_cv:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the convective part of the internal heat
            gains in the zone.
        Q_dot_ihg_rd:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the radiative part of the internal heat
            gains in the zone.
        Q_dot_sys:
            A function with signature
                f(t_sol_sec: float, T_zone: float) -> Quantity
            which takes the time `t_sol_sec` in seconds from midnight (0 s) of
            the considered day and the zone air temperature in Kelvins (K). It
            returns the cooling capacity (`Quantity` object) of the cooling
            system (i.e. the heat rate extracted from the zone air by the
            cooling system). Note that heat extraction by the cooling system
            is associated with a positive sign. Should instead the system supply
            heat to the zone, this must be associated with a negative sign.
            If this parameter is left to None, it is assumed that no (operating)
            cooling system is present in the zone, i.e. the zone is truly
            unconditioned.
        """
        self.ebe_lst = ebe_lst
        self.F_rad = F_rad
        self.A_tsn = A_tsn
        self.C_tsn = C_tsn
        self.R_tsn = R_tsn
        self.A_zan = A_zan
        self.C_zan = C_zan
        self.windows = windows
        self.ext_doors = ext_doors
        self.int_build_elems = int_build_elems
        self.ventilation = ventilation
        self.Q_dot_sol_cv = (
            Q_dot_sol_cv if Q_dot_sol_cv is not None
            else lambda t_sol_sec: Q_(0.0, 'W')
        )
        self.Q_dot_sol_rd = (
            Q_dot_sol_rd if Q_dot_sol_rd is not None
            else lambda t_sol_sec: Q_(0.0, 'W')
        )
        self.Q_dot_ihg_cv = (
            Q_dot_ihg_cv if Q_dot_ihg_cv is not None
            else lambda t_sol_sec: Q_(0.0, 'W')
        )
        self.Q_dot_ihg_rd = (
            Q_dot_ihg_rd if Q_dot_ihg_rd is not None
            else lambda t_sol_sec: Q_(0.0, 'W')
        )
        self.Q_dot_sys = (
            Q_dot_sys if Q_dot_sys is not None
            else lambda t_sol_sec, _: Q_(0.0, 'W')
        )
        self.ebe_ltn_dict: dict[str, list[TemperatureNode]] = {}
        self.tsn: ThermalStorageNode | None = None
        self.zan: ZoneAirNode | None = None

    def _create_ebe_nodes(self) -> None:
        """Takes the `ConstructionAssembly` object of each exterior building
        element in `self.ebe_lst`, which is made of `ConstructionLayer` objects,
        and creates a list of `TemperatureNode` objects from it: (1) an exterior
        surface node, (2) one or more building mass nodes and (3) an interior
        surface node. This list is then added to `self.ebe_ltn_dict`.
        """
        for ebe in self.ebe_lst:
            constr_assem = ebe.constr_assem
            ltn_lst = NodalThermalZoneModelBuilder._compose(constr_assem)
            ltn_red = NodalThermalZoneModelBuilder._reduce(ltn_lst)
            ltn = self._transform(ltn_red, ebe.net_area, ebe.ID, ebe._ext_surf.T_sa)
            self.ebe_ltn_dict[ebe.ID] = ltn

    @staticmethod
    def _compose(constr_assem: ConstructionAssembly) -> list[Quantity]:
        """Creates a list of thermal resistors and capacitors with a common
        pattern like:
        [R_in(1), C(1), R_out(1), R_in(2), C(2), R_out(2), ..., R_in(k), C(k), R_out(k)]
        """
        layers: list[ConstructionLayer] = list(constr_assem.layers.values())
        ltn = []
        # Iterate over all layers except the last one, which is always the
        # interior surface film:
        for layer in layers[:-1]:
            n = layer.num_slices
            R_slice = layer.R / (2 * n)
            C_slice = layer.C / n
            slices = [R_slice, C_slice, R_slice] * n
            ltn.extend(slices)
        # The last node in the linear thermal network is the interior surface
        # node, situated on the interior face of the exterior building element.
        # Its thermal capacity is zero, and its outgoing resistance is the
        # interior surface film resistance:
        ltn.append(Q_(0.0, 'J / (K * m**2)'))
        ltn.append(layers[-1].R)
        return ltn

    @staticmethod
    def _reduce(ltn: list[Quantity]) -> list[Quantity]:
        """Reduces `ltn` by adding resistors in series and omitting
        capacitors with a zero value (except the last C, which is the
        interior surface node).
        """
        # [R_in(1), C(1), R_(1+2), C(2), R_(2+3), ..., R_(k-1,k), C(k), R_out(k)]
        R_dummy = Q_(0, 'K * m**2 / W')
        ltn_red = []
        R = ltn[0]
        for i in range(1, len(ltn)):
            if R_dummy.check(ltn[i].dimensionality):
                # element i is a resistor: add it to R
                R += ltn[i]
            else:
                # element i is a capacitor (only keep C if it is greater than
                # zero or if it is the last C in the list).
                if ltn[i].m > 0.0 or i == len(ltn) - 2:
                    ltn_red.append(R)  # add R to the reduced list
                    ltn_red.append(ltn[i])  # add C to the reduced list
                    R = Q_(0.0, 'K * m**2 / W')  # reset R
        if R.magnitude > 0.0:
            ltn_red.append(R)
        return ltn_red

    def _transform(
        self,
        ltn_red: list[Quantity],
        A: Quantity,
        ID: str,
        T_sa: Callable[[float], Quantity]
    ) -> list[TemperatureNode]:
        """Transforms `ltn_red` into an `ExteriorBuildingElementLTN` object,
        being the linear thermal network model of the exterior building element.
        The first element in the `ExteriorBuildingElementLTN` object is always
        an exterior surface node, followed by one or more building mass nodes,
        and finally the interior surface node.
        """
        n = len(ltn_red)
        esn_lst = ltn_red[:3]
        isn_lst = ltn_red[(n - 3):n]
        esn = ExteriorSurfaceNode.create(
            ID=f'ESN {ID}',
            A=A,
            C=esn_lst[1],
            R1=esn_lst[0],
            R2=esn_lst[2],
            T_sa=T_sa
        )
        isn = InteriorSurfaceNode.create(
            ID=f'ISN {ID}',
            A=A,
            R1=isn_lst[0],
            R2=isn_lst[2],
            F_rad=self.F_rad
        )
        bmn_lst = [
            BuildingMassNode.create(
                ID=f'BMN{k + 1} {ID}',
                A=A,
                C=ltn_red[i],
                R1=ltn_red[i - 1],
                R2=ltn_red[i + 1]
            ) for k, i in enumerate(range(3, n - 2, 2))
        ]
        ltn = [esn]
        ltn.extend(bmn_lst)
        ltn.append(isn)
        return ltn

    def _create_tsn_node(self) -> None:
        """After the temperature nodes of each exterior building element are
        created, the thermal storage node of the thermal zone can be created.
        The thermal storage node is connected to the interior surface node of
        each exterior building element that surrounds the thermal zone and also
        to the zone air node of the thermal zone.
        """
        # Get the R2-resistor of each interior surface node:
        R1_lst = [ebe_ltn[-1].R2 for ebe_ltn in self.ebe_ltn_dict.values()]
        # Get the surface area associated with each interior surface node:
        A_lst = [ebe_ltn[-1].A for ebe_ltn in self.ebe_ltn_dict.values()]
        # Add the surface area of the thermal storage node to this list:
        A_lst.append(self.A_tsn)
        # Create the thermal storage node:
        self.tsn = ThermalStorageNode.create(
            ID='TSN',
            C=self.C_tsn,
            A=A_lst,
            R1=R1_lst,
            R2=self.R_tsn,
            F_rad=self.F_rad,
            windows=self.windows,
            ext_doors=self.ext_doors,
            int_build_elems=self.int_build_elems,
            Q_dot_sol_rd=self.Q_dot_sol_rd,
            Q_dot_ihg_rd=self.Q_dot_ihg_rd
        )

    def _create_zan_node(self) -> None:
        """After the temperature nodes of each exterior building element are
        created, the zone air node of the thermal zone can be created.
        The zone air node is connected to the interior surface node of
        each exterior building element that surrounds the thermal zone and also
        to the thermal storage node of the thermal zone.
        """
        # Get the R2-resistor of each interior surface node:
        R1_lst = [ebe_ltn[-1].R2 for ebe_ltn in self.ebe_ltn_dict.values()]
        # Get the surface area associated with each interior surface node:
        A_lst = [ebe_ltn[-1].A for ebe_ltn in self.ebe_ltn_dict.values()]
        # Add the surface area of the thermal storage node to this list:
        A_lst.append(self.A_tsn)
        # Add the floor area of the zone also to this list:
        A_lst.append(self.A_zan)
        # Create the zone air node:
        self.zan = ZoneAirNode.create(
            ID='ZAN',
            C=self.C_zan,
            A=A_lst,
            R1=R1_lst,
            R2=self.R_tsn,
            F_rad=self.F_rad,
            windows=self.windows,
            ext_doors=self.ext_doors,
            int_build_elems=self.int_build_elems,
            ventilation=self.ventilation,
            Q_dot_sol_cv=self.Q_dot_sol_cv,
            Q_dot_ihg_cv=self.Q_dot_ihg_cv,
            Q_dot_sys=self.Q_dot_sys
        )

    def get_nodes(self) -> tuple[dict[str, TemperatureNode], TemperatureNode, TemperatureNode]:
        """Returns (1) a dictionary with for each exterior building element
        surrounding the zone a list with the temperature nodes, (2) the thermal
        storage node, and (3) the zone air node.
        """
        self._create_ebe_nodes()
        self._create_tsn_node()
        self._create_zan_node()
        return self.ebe_ltn_dict, self.tsn, self.zan


class NodalThermalZoneModel:
    """Represents a lumped-capacitance model for calculating the conduction heat
    gain through the exterior building elements of a thermal zone when the zone
    air temperature is free-floating, and instead the cooling system capacity is
    given as a function of time.
    """
    def __init__(self) -> None:
        self.ebe_node_dict: dict[str, list[TemperatureNode]] = {}
        self.tsn: ThermalStorageNode | None = None
        self.zan: ZoneAirNode | None = None
        self.num_nodes: int = 0
        self._A: np.ndarray | None = None
        self._B: np.ndarray | None = None
        self._T_node_table: np.ndarray | None = None
        self._Q_dot_table_dict: dict[str, list[np.ndarray]] = {}

    @classmethod
    def create(
        cls,
        ext_build_elems: list[ExteriorBuildingElement],
        F_rad: Quantity = Q_(0.46, 'frac'),
        A_tsn: Quantity = Q_(1.0, 'm**2'),
        C_tsn: Quantity = Q_(0.0, 'J / (K * m**2)'),
        R_tsn: Quantity = Q_(float('inf'), 'K * m**2 / W'),
        A_zan: Quantity = Q_(1.0, 'm**2'),
        C_zan: Quantity = Q_(0.0, 'J / (K * m**2)'),
        windows: list[Window] | None = None,
        ext_doors: list[ExteriorBuildingElement] | None = None,
        int_build_elems: list[InteriorBuildingElement] | None = None,
        ventilation: SpaceVentilation | None = None,
        Q_dot_sol_cv: Callable[[float], Quantity] | None = None,
        Q_dot_sol_rd: Callable[[float], Quantity] | None = None,
        Q_dot_ihg_cv: Callable[[float], Quantity] | None = None,
        Q_dot_ihg_rd: Callable[[float], Quantity] | None = None,
        Q_dot_sys: Callable[[float, ...], Quantity] | None = None
    ) -> NodalThermalZoneModel:
        """Creates a `NodalThermalZoneModel`.

        Parameters
        ----------
        ext_build_elems:
            List with the exterior building elements surrounding the thermal
            zone.
        F_rad:
            The radiative fraction of the conduction heat gain.
        A_tsn:
            The surface area associated with the zone's interior thermal mass.
        C_tsn:
            The unit thermal capacity associated with the zone's interior thermal
            mass.
        R_tsn:
            The unit thermal resistance between the zone's interior thermal mass
            and the zone air.
        A_zan:
            The floor area of the zone.
        C_zan:
            The thermal capacity of the zone air per unit floor area.
        windows:
            List of the windows that are present in the exterior envelope of
            the zone.
        ext_doors:
            List of the exterior doors that are present in the exterior envelope
            of the zone.
        int_build_elems:
            List of the interior building elements, including interior doors,
            that surround the zone.
        ventilation:
            `SpaceVentilation` object that encapsulates the
            ventilation/infiltration model of the zone.
        Q_dot_sol_cv:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the convective part of the solar heat
            gain.
        Q_dot_sol_rd:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the radiative part of the solar heat
            gain.
        Q_dot_ihg_cv:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the convective part of the internal heat
            gains in the zone.
        Q_dot_ihg_rd:
            Function that takes the solar time in seconds from midnight (0 s) on
            the design day and returns the radiative part of the internal heat
            gains in the zone.
        Q_dot_sys:
            A function with signature
                f(t_sol_sec: float, T_zone: float) -> Quantity
            which takes the time `t_sol_sec` in seconds from midnight (0 s) of
            the considered day and the zone air temperature in Kelvins (K). It
            returns the cooling capacity (`Quantity` object) of the cooling
            system (i.e. the heat rate extracted from the zone air by the
            cooling system). Note that heat extraction by the cooling system
            is associated with a positive sign. Should instead the system supply
            heat to the zone, this must be associated with a negative sign.
            If this parameter is left to None, it is assumed that no (operating)
            cooling system is present in the zone, i.e. the zone is truly
            unconditioned.
        """
        thz_model_builder = NodalThermalZoneModelBuilder(
            ebe_lst=ext_build_elems,
            F_rad=F_rad,
            A_tsn=A_tsn,
            C_tsn=C_tsn,
            R_tsn=R_tsn,
            A_zan=A_zan,
            C_zan=C_zan,
            windows=windows,
            ext_doors=ext_doors,
            int_build_elems=int_build_elems,
            ventilation=ventilation,
            Q_dot_sol_cv=Q_dot_sol_cv,
            Q_dot_sol_rd=Q_dot_sol_rd,
            Q_dot_ihg_cv=Q_dot_ihg_cv,
            Q_dot_ihg_rd=Q_dot_ihg_rd,
            Q_dot_sys=Q_dot_sys
        )
        tup = thz_model_builder.get_nodes()
        model = cls()
        model.ebe_node_dict = tup[0]
        model.tsn = tup[1]
        model.zan = tup[2]
        # Count the number of ebe-nodes and add 2 to this number for the thermal
        # storage node and zone air node:
        num_ebe_nodes = sum(
            len(ebe_node_lst)
            for ebe_node_lst in model.ebe_node_dict.values()
        )
        model.num_nodes = num_ebe_nodes + 2
        return model

    def _build_A_matrix(self):
        # Create a square zero-matrix with size `num_nodes` x `num_nodes`:
        self._A = np.zeros((self.num_nodes, self.num_nodes))
        # Fill the A-matrix with the a-coefficients of the nodes.
        # The second-to-last row and column are for the thermal storage node.
        # The last row and column are for the zone air node.
        r, c = 0, 0  # row and column index in A
        c_max = self.num_nodes - 1  # index of the last column in A (zan)
        # Add the a-coefficients of the ebe-nodes:
        isn_index_lst = []  # keep a list with the indexes of the interior surface nodes (isn)
        for ebe_node_lst in self.ebe_node_dict.values():
            i_max = len(ebe_node_lst) - 1
            for i, ebe_node in enumerate(ebe_node_lst):
                a_coeff_lst = ebe_node.get_a_coefficients()
                if i == 0:
                    # esn
                    self._A[r, c] = a_coeff_lst[0]
                    self._A[r, c + 1] = a_coeff_lst[1]
                elif 0 < i < i_max:
                    # bmn
                    self._A[r, c - 1] = a_coeff_lst[0]
                    self._A[r, c] = a_coeff_lst[1]
                    self._A[r, c + 1] = a_coeff_lst[2]
                elif i == i_max:
                    # isn
                    self._A[r, c - 1] = a_coeff_lst[0]
                    self._A[r, c] = a_coeff_lst[1]
                    self._A[r, c_max] = a_coeff_lst[2]
                    isn_index_lst.append(c)
                r += 1
                c += 1
        # Add the a-coefficients of the thermal storage node (tsn):
        a_coeff_lst = self.tsn.get_a_coefficients()
        for j, a in enumerate(a_coeff_lst[:-2]):
            c_isn = isn_index_lst[j]
            self._A[r, c_isn] = a
        self._A[r, c_max - 1] = a_coeff_lst[-2]
        self._A[r, c_max] = a_coeff_lst[-1]
        r += 1
        # Add the a-coefficients of the zone air node (zan):
        a_coeff_lst = self.zan.get_a_coefficients()
        for j, a in enumerate(a_coeff_lst[:-2]):
            c_isn = isn_index_lst[j]
            self._A[r, c_isn] = a
        self._A[r, c_max - 1] = a_coeff_lst[-2]
        self._A[r, c_max] = a_coeff_lst[-1]

    def _build_B_matrix(self, t_sol_sec: float) -> None:
        """Builds the input matrix of the thermal network model at solar time
        moment `t_sol_sec` in seconds from midnight (0 s) of the design day.
        """
        # Create a zero-valued 1D-array with size `num_nodes`:
        self._B = np.zeros((self.num_nodes,))
        # Fill the B-matrix with the b-coefficients of the ebe-nodes:
        c = 0  # column-index of `self._T_node_table`
        for ebe_node_lst in self.ebe_node_dict.values():
            i_max = len(ebe_node_lst) - 1
            for i, ebe_node in enumerate(ebe_node_lst):
                if i == 0:
                    # esn
                    b = ebe_node.get_b_value([
                        self._T_node_table[-2, c],
                        self._T_node_table[-1, c]],
                        t_sol_sec
                    )
                    self._B[c] = b
                elif 0 < i < i_max:
                    # bmn
                    b = ebe_node.get_b_value([
                        self._T_node_table[-2, c],
                        self._T_node_table[-1, c]
                    ])
                    self._B[c] = b
                elif i == i_max:
                    # isn
                    b = ebe_node.get_b_value()
                    self._B[c] = b
                c += 1
        # Fill in the b-coefficient of the thermal storage node:
        b = self.tsn.get_b_value([
            self._T_node_table[-2, c],
            self._T_node_table[-1, c]],
            t_sol_sec
        )
        self._B[c] = b
        c += 1
        # Fill in the b-coefficient of the zone air node:
        b = self.zan.get_b_value([
            self._T_node_table[-2, c],
            self._T_node_table[-1, c]],
            t_sol_sec
        )
        self._B[c] = b

    def _init(
        self,
        init_values: list[list[Quantity]] | Quantity | None = None,
        cycle: int | None = None
    ) -> None:
        """Initializes the temperature node table at each new calculation
        cycle.

        Parameters
        ---------
        init_values: optional
            2D-list with the initial temperatures (`Quantity` objects) of the
            nodes at time index -2 (first row) and at time index -1 (second
            row). The number of columns of this 2D-list must be equal to the
            number of nodes in the linear thermal network.
            If `init_values` is None, an array of zeros will be created.
        cycle:
            The index of the current cycle.
        """
        if cycle is None or cycle == 0:
            if init_values is None:
                n_cols = self.num_nodes
                n_rows = 2
                self._T_node_table = np.zeros((n_rows, n_cols))
            elif isinstance(init_values, Quantity):
                n_cols = self.num_nodes
                n_rows = 2
                self._T_node_table = np.full((n_rows, n_cols), init_values.to('K').m)
            else:
                self._T_node_table = np.array([
                    [T.to('K').m for T in row]
                    for row in init_values
                ])
        else:
            # `cycle` > 0.
            # Initial values for the next cycle are the last two rows in
            # `_T_node_table` from the previous cycle:
            self._T_node_table = self._T_node_table[-2:, :]

    def _setup_Q_dot_table_dict(self):
        self._Q_dot_table_dict = {
            ebe_ID: []
            for ebe_ID in self.ebe_node_dict.keys()
        }

    def solve_one_step(self, t_sol_sec: float) -> None:
        """Solves the nodal thermal network model for the node temperatures and
        heat flows at the current solar time moment `t_sol_sec` in seconds from
        midnight (0 s) of the design day.
        """
        self._build_B_matrix(t_sol_sec)
        T_node_row = np.linalg.solve(self._A, self._B)
        Q_dot_dict = self._calculate_heat_flows(t_sol_sec, T_node_row)
        T_node_row = T_node_row.reshape(1, len(T_node_row))
        self._T_node_table = np.append(self._T_node_table, T_node_row, axis=0)
        self._update_Q_dot_table_dict(Q_dot_dict)

    def _calculate_heat_flows(
        self,
        t_sol_sec: float,
        T_node_row: np.ndarray
    ) -> dict[str, np.ndarray]:
        """Calculates the conduction heat transfer in the exterior building
        elements surrounding the zone at the solar time moment `t_sol_sec`
        expressed in seconds from midnight (0 s) on the design day.

        Other Parameters
        ----------------
        T_node_row: array of floats
            Numpy 1D-array of the node temperatures in units of Kelvin at time
            moment `t_sol_sec`.

        Returns
        -------
        A dictionary of which the keys are the IDs of the exterior building
        elements surrounding the zone and the values are Numpy arrays with
        the heat flows between the nodes.
        """
        Q_dot_dict = {}
        t = [0, 0]  # initial start and end position-index in T_node_row
        for ebe_ID, ebe_node_lst in self.ebe_node_dict.items():
            n = len(ebe_node_lst)
            t[1] = t[0] + n  # set the new end position-index
            Q_dot_row = np.zeros((n,))
            T_node_slice = T_node_row[t[0]:t[1]]
            i_max = n - 1
            for i, node in enumerate(ebe_node_lst):
                if i == 0:
                    # heat flow into exterior surface node
                    Q_dot_row[i] = (node._T_sa(t_sol_sec).to('K').m - T_node_slice[i]) / node._R1
                elif 0 < i < i_max:
                    # heat flow into building mass node
                    Q_dot_row[i] = (T_node_slice[i - 1] - T_node_slice[i]) / node._R1
                elif i == i_max:
                    # heat flow into interior surface node = heat flow out of
                    # interior surface node (as C = 0)
                    Q_dot_row[i] = (T_node_slice[i - 1] - T_node_slice[i]) / node._R1
            t[0] = t[1]  # move the start-position index to the current end-position index
            Q_dot_dict[ebe_ID] = Q_dot_row
        return Q_dot_dict

    def _update_Q_dot_table_dict(self, Q_dot_dict: dict[str, np.ndarray]):
        for ebe_ID, Q_dot_row in Q_dot_dict.items():
            self._Q_dot_table_dict[ebe_ID].append(Q_dot_row)

    def solve(
        self,
        num_steps: int,
        init_values: list[list[Quantity]] | Quantity | None = None,
        dt_hr: float = 1.0,
        num_cycles: int = 1
    ) -> None:
        """Solves the nodal thermal network model for the node temperatures in
        the course of time.

        Parameters
        ----------
        num_steps:
            The number of time steps in one cycle at which the node temperatures
            are calculated.
        init_values: optional
            2D-list with the initial temperatures (`Quantity` objects) of the
            nodes at time index -2 (first row) and at time index -1 (second
            row) before the first cycle. The number of columns of this 2D-list
            must be equal to the number of nodes in the thermal model.
            If `init_values` is None, an array of zeros will be created.
        dt_hr: optional
            The time step width in hours between two successive time moments
            at which the thermal model is solved. The default value
            is 1 hr. The product of the number of time steps per cycle and the
            time step width determines the duration (period) of one cycle.
        num_cycles:
            The number of times the cycle of `num_steps` calculations is repeated.
            Each new cycle starts with the node temperatures from the last two
            time indexes k-2 and k-1 of the previous cycle as the initial values.
            Only the last cycle is kept.

        Returns
        -------
        None
        """
        set_time_step(dt_hr)
        self._build_A_matrix()
        for cycle in range(num_cycles):
            self._init(init_values, cycle)
            self._setup_Q_dot_table_dict()
            for k in range(num_steps):
                self.solve_one_step(k * dt_sec)

    def get_node_temperatures(self, unit: str = 'degC') -> dict[str, pd.DataFrame]:
        """Returns for each exterior building element surrounding the zone, for
        the zone air, and for the zone's interior thermal mass a Pandas DataFrame
        of which the column titles are the ID's of the nodes in the nodal thermal
        network, ordered from the exterior surface node on the left to the
        interior surface node on the right in the case of an exterior building
        element. Each row contains the values of the node temperatures
        calculated at each time index in the last calculation cycle (see
        parameter `num_cycles` of method `solve()`).
        """
        d = {ebe_ID: [] for ebe_ID in self.ebe_node_dict.keys()}
        for T_node_row in self._T_node_table[2:, :-2]:
            # The first two rows are from the second-to-last cycle (these were
            # the initial values to start the last calculation cycle).
            t = [0, 0]  # initial start and end position-index in T_node_row
            for ebe_ID, ebe_node_lst in self.ebe_node_dict.items():
                n = len(ebe_node_lst)
                t[1] = t[0] + n  # set the new end position-index
                T_node_slice = T_node_row[t[0]:t[1]]
                d[ebe_ID].append(T_node_slice)
                t[0] = t[1]
        df_dict = {}
        for ebe_ID, ebe_node_lst in self.ebe_node_dict.items():
            data = Q_(d[ebe_ID], 'K').to(unit).m
            columns = [ebe_node.ID for ebe_node in ebe_node_lst]
            df_dict[ebe_ID] = pd.DataFrame(data=data, columns=columns)
        tsn_data = self._T_node_table[2:, -2]
        df_dict['TSN'] = pd.DataFrame(data=Q_(tsn_data, 'K').to(unit).m, columns=['TSN'])
        zan_data = self._T_node_table[2:, -1]
        df_dict['ZAN'] = pd.DataFrame(data=Q_(zan_data, 'K').to(unit).m, columns=['ZAN'])
        return df_dict

    def get_heat_flows(self, unit: str = 'W') -> dict[str, pd.DataFrame]:
        """Returns a dictionary of which the keys are the IDs of the exterior
        building elements. The keys map to "heat flow tables" (Pandas
        `DataFrame` objects) of the exterior building elements.
        (A heat flow table shows the heat flow rates into the nodes at each time
        index of the last calculation cycle of the design day.)
        The last column shows the heat that flows into the interior surface node
        of the exterior building element. This is also the heat that flows out
        to the zone air node (convective component) and to the thermal storage
        node (radiative component).
        """
        df_dict = {}
        for key, data in self._Q_dot_table_dict.items():
            ebe_node_lst = self.ebe_node_dict[key]
            columns = [ebe_node.ID for ebe_node in ebe_node_lst]
            df_dict[key] = pd.DataFrame(Q_(data, 'W').to(unit).m, columns=columns)
        return df_dict

    def conductive_heat_gain(self, unit: str = 'W') -> dict[str, Quantity]:
        """Returns a dictionary of which the keys are the IDs of the exterior
        building elements and the values are the conduction heat gain as
        function of time (`Quantity` array) at the interior side of these
        exterior building elements.
        """
        df_dict = self.get_heat_flows(unit)
        d = {k: Q_(df.iloc[:, -1].values, unit) for k, df in df_dict.items()}
        return d

    def zone_air_temperature(self, unit: str = 'degC') -> Quantity:
        """Returns the zone air temperature as function of time
        (`Quantity` array).
        """
        df_dict = self.get_node_temperatures(unit)
        T_zone = Q_(df_dict['ZAN'].iloc[:, -1].values, unit)
        return T_zone

    def thermal_storage_temperature(self, unit: str = 'degC') -> Quantity:
        """Returns the interior thermal mass temperature as function of time
        (`Quantity` array).
        """
        df_dict = self.get_node_temperatures(unit)
        T_tsn = Q_(df_dict['TSN'].iloc[:, -1].values, unit)
        return T_tsn
