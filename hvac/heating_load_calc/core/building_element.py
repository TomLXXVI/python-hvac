from __future__ import annotations

from typing import TypeVar
from abc import ABC, abstractmethod
from hvac import Quantity
from hvac.cooling_load_calc.building import (
    ConstructionAssembly, 
    WindowThermalProperties
)


Q_ = Quantity


TBuildingElement = TypeVar(
    'TBuildingElement', 
    bound='AbstractBuildingElement'
)


class AbstractBuildingElement(ABC):

    def __init__(self):
        self.ID: str = ''
        self._area: Quantity | None = None
        self.constr_assem: ConstructionAssembly | None = None
        self.T_int_d: Quantity | None = None
        self.T_ext_d: Quantity | None = None
        self.T_adj: Quantity | None = None
        self.T_int_surf: Quantity | None = None
        self.f1: Quantity | None = None
        self.f2: Quantity | None = None
        self.f_T: Quantity | None = None  # temperature adjustment factor
        self.building_elements: dict[str, TBuildingElement] = {}

    @staticmethod
    def _create(
        obj: TBuildingElement,
        ID: str,
        area: Quantity | tuple[Quantity, Quantity],
        constr_assem: ConstructionAssembly | WindowThermalProperties,
        T_int_d: Quantity,
        T_int_surf: Quantity,
        T_ext_d: Quantity,
        T_adj: Quantity | None = None,
        f1: Quantity | None = None
    ) -> TBuildingElement:
        """Creates a building element.

        Parameters
        ----------
        obj: 
            Empty object of type `TBuildingElement` to which the parameters 
            below will be assigned.
        ID: str
            Name to identify building element.
        area: Quantity or tuple of 2 Quantities
            Area of the building element. In case a tuple is given, the elements
            are considered to be the width and length or height of the building
            element.
        constr_assem: ConstructionAssembly
            The construction assembly that constitutes the building element.
        T_int_d: Quantity
            Design value of indoor air temperature of heated space.
        T_int_surf: Quantity
            Interior surface temperature of building element.
        T_ext_d: Quantity
            Design value of outdoor air temperature.
        T_adj: Quantity, default None
            Temperature on the other side of the building element.
            If other side is...
            --> the exterior:
                T_ext_d (i.e. the default, you can skip setting parameter `T_adj`)
            --> an adjacent heated space within the same building entity:
                T_int_d of adjacent space.
            --> an adjacent building entity within the same building:
                T_ext_an (acc. to EN 12831-1 ANB NA.5.4)
            --> the ground:
                T_ext_an (acc. to EN 12831-1 §6.3.2.5 - Table 7)
            --> an adjacent unheated space:
                set parameter `f1` acc. to EN 12831-1, B.2.4 - Table B.2 or
                use EN 12831-1 ANB NA.5.5 - Table NA.4
            --> an adjacent building:
                max(T_ext_an, 5 °C) or use EN 12831-1 ANB NA.5.6 - Table NA.5
        f1: Quantity, default None
            Adjustment factor for differences between the temperature of an
            adjacent space and the external design temperature.
            (EN 12831-1, B.2.4 - Table B.2)
        """
        obj.ID = ID
        obj._area = area
        obj.constr_assem = constr_assem
        obj.T_int_d = T_int_d
        obj.T_int_surf = T_int_surf
        obj.T_ext_d = T_ext_d
        obj.T_adj = T_adj if T_adj is not None else T_ext_d
        obj.f1 = f1 or (obj.T_int_d - obj.T_adj) / (obj.T_int_d - obj.T_ext_d)
        obj.f2 = (obj.T_int_surf - obj.T_int_d) / (obj.T_int_d - obj.T_ext_d)
        obj.f_T = obj.f1 + obj.f2
        return obj

    @property
    def area_gross(self) -> Quantity:
        """Gross area of the building element (including windows and doors)."""
        A = self._area[0] * self._area[1] if isinstance(self._area, tuple) else self._area
        return A

    @property
    def area(self) -> Quantity:
        """Net area of the building element (without windows and doors)."""
        A = self.area_gross
        A -= sum(be.area for be in self.building_elements.values()) or Q_(0.0, 'm ** 2')
        return A

    @property
    def A(self) -> Quantity:
        """Net area of the building element (without windows and doors)."""
        return self.area

    @property
    def U(self) -> Quantity:
        """Thermal transmittance of the building element."""
        return self.constr_assem.U

    @property
    @abstractmethod
    def H(self) -> Quantity:
        """Heat transfer coefficient of building element."""
        ...


class ExteriorBuildingElement(AbstractBuildingElement):

    def __init__(self):
        super().__init__()
        self.dU_tb: Quantity | None = None
        self.f_U: Quantity | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        area: Quantity | tuple[Quantity, Quantity],
        constr_assem: ConstructionAssembly | WindowThermalProperties,
        T_int_d: Quantity,
        T_int_surf: Quantity,
        T_ext_d: Quantity,
        dU_tb: Quantity = Q_(0.1, 'W / (m ** 2 * K)'),
        f_U: Quantity = Q_(1.0, 'frac')
    ) -> ExteriorBuildingElement:
        """Creates a building element in contact with the exterior.

        Parameters
        ----------
        ID: str
            Name to identify building element.
        area: Quantity or tuple of 2 Quantities
            Area of the building element. In case a tuple is given, the elements
            are considered to be the width and length or height of the building
            element.
        constr_assem: ConstructionAssembly
            The construction assembly that constitutes the building element.
        T_int_d: Quantity
            Design value of indoor air temperature of heated space.
        T_int_surf: Quantity
            Interior surface temperature of building element.
        T_ext_d: Quantity
            Design value of outdoor air temperature.
        dU_tb: Quantity, default 0.1 W/(m².K)
            Blanket additional thermal transmittance for thermal bridges
            (see EN 12831-1, B.2.1 - Table B.1).
        f_U: Quantity, default 1.0 frac:
            Correction factor allowing for the influence of building element
            properties and meteorological conditions (see EN 12831-1, B.2.2).
        """
        self = cls()
        self = cls._create(
            self, ID, area,
            constr_assem,
            T_int_d, T_int_surf, T_ext_d
        )
        self.dU_tb = dU_tb
        self.f_U = f_U
        return self

    def add_building_element(
        self,
        ID: str,
        area: Quantity | tuple[Quantity, Quantity],
        constr_assem: ConstructionAssembly | WindowThermalProperties,
        dU_tb: Quantity = Q_(0.1, 'W / (m ** 2 * K)'),
        f_U: Quantity = Q_(1.0, 'frac')
    ) -> None:
        """Adds a "sub-building element" (window or door) to the exterior
        building element.
        """
        ext_build_elem = ExteriorBuildingElement.create(
            ID=ID,
            area=area,
            constr_assem=constr_assem,
            T_int_d=self.T_int_d,
            T_int_surf=self.T_int_surf,
            T_ext_d=self.T_ext_d,
            dU_tb=dU_tb,
            f_U=f_U
        )
        self.building_elements[ext_build_elem.ID] = ext_build_elem

    @property
    def H(self) -> Quantity:
        H = self.A * (self.U + self.dU_tb) * self.f_U * self.f_T
        H += sum(be.H for be in self.building_elements.values()) or Q_(0.0, 'W / K')
        return H.to('W / K')


class AdjacentBuildingElement(AbstractBuildingElement):

    def __init__(self):
        super().__init__()
        self.kind_of_adj_space: str = 'unheated'

    @classmethod
    def create(
        cls,
        ID: str,
        area: Quantity | tuple[Quantity, Quantity],
        constr_assem: ConstructionAssembly,
        T_int_d: Quantity,
        T_int_surf: Quantity,
        T_ext_d: Quantity,
        kind_of_adjacent_space: str = 'unheated',
        T_adj: Quantity | None = None,
        f1: Quantity | None = None
    ) -> AdjacentBuildingElement:
        """Creates an adjacent building element, which separates the space from
        an adjacent space.

        Parameters
        ----------
        ID: str
            Name to identify building element.
        area: Quantity or tuple of 2 Quantities
            Area of the building element. In case a tuple is given, the elements
            are considered to be the width and length or height of the building
            element.
        constr_assem: ConstructionAssembly
            The construction assembly that constitutes the building element.
        T_int_d: Quantity
            Design value of indoor air temperature of heated space.
        T_int_surf: Quantity
            Interior surface temperature of building element.
        T_ext_d: Quantity
            Design value of outdoor air temperature.
        kind_of_adjacent_space: {'heated', 'unheated', 'building_entity'},
            default 'unheated'
            The kind of adjacent space. In case the adjacent space is in another
            building, use 'unheated'. In case the adjacent space is in another
            building entity of the same building, use 'building_entity'.
        T_adj: Quantity
            Temperature on the other side of the building element.
            If other side is...
            - an adjacent heated space within the same building entity:
                T_int_d of adjacent space.
            - an adjacent building entity within the same building:
                T_ext_an (acc. to EN 12831-1 ANB NA.5.4)
            - an adjacent unheated space:
                set parameter `f1` acc. to EN 12831-1, B.2.4 - Table B.2 or
                use EN 12831-1 ANB NA.5.5 - Table NA.4
            - an adjacent building:
                max(T_ext_an, 5 °C) or use EN 12831-1 ANB NA.5.6 - Table NA.5
        f1: Quantity, default None
            Adjustment factor for differences between the temperature of an
            adjacent space and the external design temperature.
            (EN 12831-1, B.2.4 - Table B.2)
        """
        self = cls()
        self = cls._create(
            self, ID, area,
            constr_assem,
            T_int_d, T_int_surf,
            T_ext_d, T_adj, f1
        )
        self.kind_of_adj_space = kind_of_adjacent_space
        return self

    def add_building_element(
        self,
        ID: str,
        area: Quantity | tuple[Quantity, Quantity],
        constr_assem: ConstructionAssembly | WindowThermalProperties,
    ) -> None:
        """Adds a "sub-building element" (door or window) to the adjacent
        building element.
        """
        adj_build_elem = AdjacentBuildingElement.create(
            ID=ID,
            area=area,
            constr_assem=constr_assem,
            T_int_d=self.T_int_d,
            T_int_surf=self.T_int_surf,
            T_ext_d=self.T_ext_d,
            kind_of_adjacent_space=self.kind_of_adj_space,
            T_adj=self.T_adj,
            f1=self.f1
        )
        self.building_elements[adj_build_elem.ID] = adj_build_elem

    @property
    def H(self) -> Quantity:
        H = self.A * self.U * self.f_T
        H += sum(be.H for be in self.building_elements.values()) or Q_(0.0, 'W / K')
        return H.to('W / K')


class GroundBuildingElement(AbstractBuildingElement):

    def __init__(self):
        super().__init__()
        self.A_slab: Quantity | None = None
        self.P_slab: Quantity | None = None
        self.dU_tb: Quantity | None = None
        self.z: Quantity | None = None
        self.U_equiv: Quantity | None = None
        self.f_dT_an: Quantity | None = None
        self.f_gw: Quantity | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        area: Quantity,
        constr_assem: ConstructionAssembly,
        T_int_d: Quantity,
        T_int_surf: Quantity,
        T_ext_d: Quantity,
        A_slab: Quantity,
        P_slab: Quantity,
        z: Quantity,
        dU_tb: Quantity = Q_(0.1, 'W / (m ** 2 * K)'),
        f_dT_an: Quantity = Q_(1.45, 'frac'),
        f_gw: Quantity = Q_(1.15, 'frac')
    ) -> GroundBuildingElement:
        """Creates a building element in contact with ground.

        Parameters
        ----------
        ID: str
            Name to identify building element.
        area: Quantity
            Floor area of the heated space.
        constr_assem: ConstructionAssembly
            The construction assembly that constitutes the building element.
        T_int_d: Quantity
            Design value of indoor air temperature of heated space.
        T_int_surf: Quantity
            Interior surface temperature of building element.
        T_ext_d: Quantity
            Design value of outdoor air temperature.
        A_slab: Quantity
            Area of the floor slab (see EN 12831-1, annex E).
        P_slab: Quantity
            Exposed perimeter of the floor slab (see EN 12831-1, annex E -
            Figure E.2: examples).
        z: Quantity
            Depth of the top edge of the floor slab below ground level.
        dU_tb: Quantity
            Blanket additional thermal transmittance for thermal bridges
            (see EN 12831-1, B.2.1 - Table B.1).
        f_dT_an: Quantity
            Correction factor taking into account the annual variation of the
            external temperature (see EN 12831-1, B.2.3).
        f_gw: Quantity
            Correction factor taking into account the influence of groundwater
            (see EN 12831-1, B.2.3).
        """
        self = cls()
        self = cls._create(
            self, ID, area,
            constr_assem,
            T_int_d, T_int_surf, T_ext_d
        )
        self.A_slab = A_slab
        self.P_slab = P_slab
        self.z = z
        self.dU_tb = dU_tb
        self.U_equiv = self._calculate_U_equiv()
        self.f_dT_an = f_dT_an
        self.f_gw = f_gw
        return self

    def _calculate_U_equiv(self):
        if self.z == 0.0:  # floor
            a = 0.9671
            b = -7.455
            c = (10.76, 9.773, 0.0265)
            d = -0.0203
            n = (0.5532, 0.6027, -0.9296)
        else:  # basement wall (z > 0)
            a = 0.93328
            b = -2.1552
            c = (0.0, 1.466, 0.1006)
            d = -0.0692
            n = (0.0, 0.45325, -1.0068)
        B = self.A_slab / (0.5 * self.P_slab)
        t_B = (c[0] + B.to('m').m) ** n[0]
        t_z = (c[1] + self.z.to('m').m) ** n[1]
        t_U = (
            c[2] + self.U.to('W / (m ** 2 * K)').m +
            self.dU_tb.to('W / (m ** 2 * K)').m
        )
        U_equiv = Q_(a / (b + t_B + t_z + t_U) + d, 'W / (m ** 2 * K)')
        return U_equiv

    @property
    def H(self) -> Quantity:
        H = self.f_dT_an.m * self.A * self.U_equiv * self.f_T * self.f_gw.m
        return H.to('W / K')
