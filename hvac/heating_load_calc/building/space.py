from __future__ import annotations

from typing import TYPE_CHECKING
from hvac import Quantity
from hvac.cooling_load_calc import ConstructionAssembly
from ..core import (
    TBuildingElement,
    ExteriorBuildingElement,
    AdjacentBuildingElement,
    GroundBuildingElement
)

if TYPE_CHECKING:
    from .ventilation_zone import VentilationZone

Q_ = Quantity


class UnheatedSpace:

    def __init__(self):
        self.ID: str = ''
        self.height: Quantity | None = None
        self.area: Quantity | None = None
        self.volume: Quantity | None = None
        self.T_int_d: Quantity = Q_(20.0, 'degC')
        self.T_ext_d: Quantity = Q_(-10.0, 'degC')
        self.grad_T_air: Quantity = Q_(1.0, 'K / m')
        self.height_occ_zone: Quantity = Q_(1.0, 'm')
        self.dT_surf: Quantity = Q_(0.0, 'K')
        self.dT_rad: Quantity = Q_(0.0, 'K')
        self.n_min: Quantity = Q_(0.5, '1 / hr')
        self.V_atd_d: Quantity | None = None
        self.V_exh: Quantity | None = None
        self.V_comb: Quantity | None = None
        self.V_sup: Quantity | None = None
        self.V_trf: Quantity | None = None
        self.V_open: Quantity | None = None
        self.T_trf: Quantity | None = None
        self.vz: VentilationZone | None = None
        self.building_elements: dict[str, list[TBuildingElement]] = {
            'exterior': [],
            'adjacent_heated_space': [],
            'adjacent_unheated_space': [],
            'adjacent_building_entity': [],
            'ground': []
        }

    @classmethod
    def create(
        cls,
        ID: str,
        height: Quantity,
        area: Quantity,
        volume: Quantity | None = None,
        T_int_d: Quantity = Q_(20.0, 'degC'),
        T_ext_d: Quantity = Q_(-10.0, 'degC'),
        grad_T_air: Quantity = Q_(1.0, 'K / m'),
        height_occ_zone: Quantity = Q_(1.0, 'm'),
        dT_surf: Quantity = Q_(0.0, 'K'),
        dT_rad: Quantity = Q_(0.0, 'K'),
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_atd_d: Quantity | None = None,
        V_exh: Quantity | None = None,
        V_comb: Quantity | None = None,
        V_sup: Quantity | None = None,
        V_trf: Quantity | None = None,
        V_open: Quantity | None = None,
        T_trf: Quantity | None = None,
        vz: VentilationZone | None = None
    ) -> UnheatedSpace:
        """Creates an unheated space.

        Parameters
        ----------
        ID: str
            Name of heated space.
        height: Quantity
            Mean height of the heated space.
        area: Quantity
            Floor area of the heated space.
        volume: Quantity, default None
            Volume of the heated space. If not given, will be calculated as the
            product of `height` and `area`.
        T_int_d: Quantity, default 20 degC
            Design value of indoor air temperature of the heated space.
        T_ext_d: Quantity, default -10 degC
            Design value of outdoor air temperature.
        grad_T_air: Quantity, default 1 K/m
            Air temperature gradient of the heat emission system used in the room.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        height_occ_zone: Quantity, default 1 m
            Height of the occupied zone in the heated space.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        dT_surf: Quantity, default 0.0 K
            Correction term for the influence of the heat emission system on
            surface temperatures.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        dT_rad: Quantity, default 0.0 K
            Difference between air and operative temperature (which can be
            approximated as the arithmetic average of air and mean radiant
            temperature).
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        n_min: Quantity, default 0.5 1/hr
            Minimum air change rate required for the heated space for reasons of
            air quality/hygiene and comfort (EN 12831-1, B.2.10 - Table B.7).
            Default value applies to permanent dwelling areas (living rooms,
            offices) and a ceiling height less than 3 m.
        V_atd_d: Quantity, default None
            Design air volume flow rate of the ATDs in the room
            (EN 12831-1, B.2.12).
            Only required if ATDs are used for ventilation.
        V_exh: Quantity, default None
            Exhaust ventilation air volume flow rate from the heated space.
        V_comb: Quantity, default None
            Air volume flow rate exhausted from the heated space that has not
            been included in the exhaust air volume flow of the ventilation
            system (typically, but not necessarily, combustion air if an open
            flue heater is present in the heated space).
        V_sup: Quantity, default None
            Supply air volume flow rate from the ventilation system into the
            heated space.
        V_trf: Quantity, default None
            Transfer air volume flow rate into the heated space from adjacent
            spaces.
        V_open: Quantity, default None
            External air volume flow rate into the heated space through large
            openings (EN 12831-1, Annex G).
        T_trf: Quantity, default None
            Temperature of the transfer air volume flow into the heated space
            from another space. In case the room height of the other space is
            less than 4 m, it is equal to the internal design temperature of
            the other space; otherwise, it is equal to mean air temperature of
            the other space (see EN 12831-1 §6.3.8.3).
        vz: VentilationZone, default None
            The ventilation zone of which the heated space is part of.
        """
        self = cls()
        self.ID = ID
        self.height = height
        self.area = area
        self.volume = volume or self.height * self.area
        self.T_int_d = T_int_d
        self.T_ext_d = T_ext_d
        self.grad_T_air = grad_T_air
        self.height_occ_zone = height_occ_zone
        self.dT_surf = dT_surf
        self.dT_rad = dT_rad
        self.n_min = n_min
        self.V_atd_d = V_atd_d or Q_(0.0, 'm ** 3 / hr')
        self.V_exh = V_exh or Q_(0.0, 'm ** 3 / hr')
        self.V_comb = V_comb or Q_(0.0, 'm ** 3 / hr')
        self.V_sup = V_sup or Q_(0.0, 'm ** 3 / hr')
        self.V_trf = V_trf or Q_(0.0, 'm ** 3 / hr')
        self.V_open = V_open or Q_(0.0, 'm ** 3 / hr')
        self.T_trf = T_trf if T_trf is not None else self.T_int_air
        self.vz = vz
        return self

    # noinspection PyIncorrectDocstring
    def add_exterior_building_element(
        self,
        ID: str,
        area: Quantity | tuple[Quantity, Quantity],
        _, **kwargs
    ) -> ExteriorBuildingElement:
        """Adds a new exterior building element to the heated space.

        Parameters
        ----------
        ID: str
            Name to identify building element.
        area: Quantity or tuple of 2 Quantities
            Area of the building element. In case a tuple is given, the elements
            are considered to be the width and length or height of the building
            element.

        Returns
        -------
        The newly added exterior building element (to allow for adding doors
        and/or windows to the building element).
        """
        ext_build_elem = ExteriorBuildingElement.create(
            ID=ID,
            area=area,
            constr_assem=ConstructionAssembly(),
            T_int_d=self.T_int_d,
            T_int_surf=self.T_int_surf,
            T_ext_d=self.T_ext_d,
            dU_tb=Q_(0.0, 'W / (m ** 2 * K)'),
            f_U=Q_(1.0, 'frac')
        )
        self.building_elements['exterior'].append(ext_build_elem)
        return ext_build_elem

    @property
    def T_int_air(self):
        """Mean internal air temperature."""
        if self.height.to('m').m >= 4.0:
            return (
                self.T_int_d.to('K') +
                self.grad_T_air * (self.height / 2 - self.height_occ_zone) +
                self.dT_rad
            )
        else:
            return self.T_int_d

    @property
    def T_int_surf(self):
        """Mean internal surface temperature."""
        if self.height.to('m').m >= 4.0:
            return (
                self.T_int_d.to('K') +
                self.grad_T_air * (self.height - self.height_occ_zone) +
                self.dT_surf
            )
        else:
            return self.T_int_d

    @property
    def A_env(self) -> Quantity:
        """Envelope area of the heated space."""
        A_env = (
            sum(be.area_gross for be in self.building_elements['exterior'])
            or Q_(0.0, 'm ** 2')
        )
        return A_env

    @property
    def V_leak_atd(self) -> Quantity:
        """External air volume flow rate into the heated space through leakages
        and ATDs.
        """
        if self.vz.V_atd_d:
            V_leak_atd = (
                self.vz.V_leak * (self.A_env / self.vz.A_env) +
                self.vz.V_atd * (self.V_atd_d / self.vz.V_atd_d)
            )
        else:
            V_leak_atd = self.vz.V_leak * (self.A_env / self.vz.A_env)
        return V_leak_atd

    @property
    def V_env(self) -> Quantity:
        """External air volume flow rate into the heated space through the
        building envelope.
        """
        V_env = (
            (self.vz.V_inf_add / self.vz.V_env) *
            min(self.vz.V_env, self.V_leak_atd * self.vz.f_dir) +
            (self.vz.V_env - self.vz.V_inf_add) * self.V_leak_atd /
            self.vz.V_env
        )
        return V_env

    @property
    def V_tech(self) -> Quantity:
        """Technical air volume flow rate into the heated space."""
        V_tech = max(self.V_sup + self.V_trf, self.V_exh + self.V_comb)
        return V_tech

    @property
    def V_min(self) -> Quantity:
        """Minimum air volume flow required for reasons of air quality/hygiene
        and comfort.
        """
        V_min = self.n_min * self.volume
        return V_min

    @property
    def T_sup(self) -> Quantity:
        return self.vz.T_sup or self.T_int_air

    def get_ventilation_heat_loss(self) -> Quantity:
        """Get ventilation heat loss of the heated space (including outdoor
        air infiltration).
        """
        rho_cp = Q_(0.34, 'W * hr / (m ** 3 * K)')
        V_ext = max(self.V_env + self.V_open, self.V_min - self.V_tech)
        Q_ven = rho_cp * (
            V_ext * (self.T_int_air - self.T_ext_d)
            + self.V_sup * (self.T_int_air - self.T_sup)
            + self.V_trf * (self.T_int_air - self.T_trf)
        )
        return Q_ven


class HeatedSpace(UnheatedSpace):

    def __init__(self):
        super().__init__()
        self.q_hu: Quantity | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        height: Quantity,
        area: Quantity,
        volume: Quantity | None = None,
        T_int_d: Quantity = Q_(20.0, 'degC'),
        T_ext_d: Quantity = Q_(-10.0, 'degC'),
        grad_T_air: Quantity = Q_(1.0, 'K / m'),
        height_occ_zone: Quantity = Q_(1.0, 'm'),
        dT_surf: Quantity = Q_(0.0, 'K'),
        dT_rad: Quantity = Q_(0.0, 'K'),
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_atd_d: Quantity | None = None,
        V_exh: Quantity | None = None,
        V_comb: Quantity | None = None,
        V_sup: Quantity | None = None,
        V_trf: Quantity | None = None,
        V_open: Quantity | None = None,
        T_trf: Quantity | None = None,
        q_hu: Quantity | None = None,
        vz: VentilationZone | None = None
    ) -> HeatedSpace:
        """Creates a heated space.

        Parameters
        ----------
        ID: str
            Name of heated space.
        height: Quantity
            Mean height of the heated space.
        area: Quantity
            Floor area of the heated space.
        volume: Quantity, default None
            Volume of the heated space. If not given, will be calculated as the
            product of `height` and `area`.
        T_int_d: Quantity, default 20 degC
            Design value of indoor air temperature of the heated space.
        T_ext_d: Quantity, default -10 degC
            Design value of outdoor air temperature.
        grad_T_air: Quantity, default 1 K/m
            Air temperature gradient of the heat emission system used in the room.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        height_occ_zone: Quantity, default 1 m
            Height of the occupied zone in the heated space.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        dT_surf: Quantity, default 0.0 K
            Correction term for the influence of the heat emission system on
            surface temperatures.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        dT_rad: Quantity, default 0.0 K
            Difference between air and operative temperature (which can be
            approximated as the arithmetic average of air and mean radiant
            temperature).
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        n_min: Quantity, default 0.5 1/hr
            Minimum air change rate required for the heated space for reasons of
            air quality/hygiene and comfort (EN 12831-1, B.2.10 - Table B.7).
            Default value applies to permanent dwelling areas (living rooms,
            offices) and a ceiling height less than 3 m.
        V_atd_d: Quantity, default None
            Design air volume flow rate of the ATDs in the room
            (EN 12831-1, B.2.12).
            Only required if ATDs are used for ventilation.
        V_exh: Quantity, default None
            Exhaust ventilation air volume flow rate from the heated space.
        V_comb: Quantity, default None
            Air volume flow rate exhausted from the heated space that has not
            been included in the exhaust air volume flow of the ventilation
            system (typically, but not necessarily, combustion air if an open
            flue heater is present in the heated space).
        V_sup: Quantity, default None
            Supply air volume flow rate from the ventilation system into the
            heated space.
        V_trf: Quantity, default None
            Transfer air volume flow rate into the heated space from adjacent
            spaces.
        V_open: Quantity, default None
            External air volume flow rate into the heated space through large
            openings (EN 12831-1, Annex G).
        T_trf: Quantity, default None
            Temperature of the transfer air volume flow into the heated space
            from another space. In case the room height of the other space is
            less than 4 m, it is equal to the internal design temperature of
            the other space; otherwise, it is equal to mean air temperature of
            the other space (see EN 12831-1 §6.3.8.3).
        q_hu: Quantity, default None
            Additional heating-up power per unit floor area
            (EN 12831-1, Annex F).
        vz: VentilationZone, default None
            The ventilation zone of which the heated space is part of.
        """
        self = cls()
        self.ID = ID
        self.height = height
        self.area = area
        self.volume = volume or self.height * self.area
        self.T_int_d = T_int_d
        self.T_ext_d = T_ext_d
        self.grad_T_air = grad_T_air
        self.height_occ_zone = height_occ_zone
        self.dT_surf = dT_surf
        self.dT_rad = dT_rad
        self.n_min = n_min
        self.V_atd_d = V_atd_d or Q_(0.0, 'm ** 3 / hr')
        self.V_exh = V_exh or Q_(0.0, 'm ** 3 / hr')
        self.V_comb = V_comb or Q_(0.0, 'm ** 3 / hr')
        self.V_sup = V_sup or Q_(0.0, 'm ** 3 / hr')
        self.V_trf = V_trf or Q_(0.0, 'm ** 3 / hr')
        self.V_open = V_open or Q_(0.0, 'm ** 3 / hr')
        self.T_trf = T_trf if T_trf is not None else self.T_int_air
        self.vz = vz
        self.q_hu = q_hu or Q_(0.0, 'W / m ** 2')
        return self

    def add_exterior_building_element(
        self,
        ID: str,
        area: Quantity | tuple[Quantity, Quantity],
        constr_assem: ConstructionAssembly,
        dU_tb: Quantity = Q_(0.1, 'W / (m ** 2 * K)'),
        f_U: Quantity = Q_(1.0, 'frac')
    ) -> ExteriorBuildingElement:
        """Adds a new exterior building element to the heated space.

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
        dU_tb: Quantity, default 0.1 W/(m².K)
            Blanket additional thermal transmittance for thermal bridges
            (see EN 12831-1, B.2.1 - Table B.1).
        f_U: Quantity, default 1.0 frac:
            Correction factor allowing for the influence of building element
            properties and meteorological conditions (see EN 12831-1, B.2.2).

        Returns
        -------
        The newly added exterior building element (to allow for adding doors
        and/or windows to the building element).
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
        self.building_elements['exterior'].append(ext_build_elem)
        return ext_build_elem

    def add_adjacent_building_element(
        self,
        ID: str,
        area: Quantity | tuple[Quantity, Quantity],
        constr_assem: ConstructionAssembly,
        kind_of_adjacent_space: str = 'unheated',
        T_adj: Quantity | None = None,
        f1: Quantity | None = None
    ) -> AdjacentBuildingElement:
        """Adds a new adjacent building element to the heated space.

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

        Returns
        -------
        The newly added adjacent building element (to allow for adding doors to
        the building element).
        """
        adj_build_elem = AdjacentBuildingElement.create(
            ID=ID,
            area=area,
            constr_assem=constr_assem,
            T_int_d=self.T_int_d,
            T_int_surf=self.T_int_surf,
            T_ext_d=self.T_ext_d,
            kind_of_adjacent_space=kind_of_adjacent_space,
            T_adj=T_adj,
            f1=f1
        )
        if adj_build_elem.kind_of_adj_space == 'heated':
            self.building_elements['adjacent_heated_space'].append(adj_build_elem)
        elif adj_build_elem.kind_of_adj_space == 'unheated':
            self.building_elements['adjacent_unheated_space'].append(adj_build_elem)
        elif adj_build_elem.kind_of_adj_space == 'building_entity':
            self.building_elements['adjacent_building_entity'].append(adj_build_elem)
        return adj_build_elem

    def add_ground_building_element(
        self,
        ID: str,
        area: Quantity,
        constr_assem: ConstructionAssembly,
        A_slab: Quantity,
        P_slab: Quantity,
        z: Quantity,
        dU_tb: Quantity = Q_(0.1, 'W / (m ** 2 * K)'),
        f_dT_an: Quantity = Q_(1.45, 'frac'),
        f_gw: Quantity = Q_(1.15, 'frac')
    ) -> None:
        """Adds a new building element in contact with ground to the heated
        space.

        Parameters
        ----------
        ID: str
            Name to identify building element.
        area: Quantity
            Floor area of the heated space.
        constr_assem: ConstructionAssembly
            The construction assembly that constitutes the building element.
        A_slab: Quantity
            Area of the floor slab.
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

        Returns
        -------
        None
        """
        grd_build_elem = GroundBuildingElement.create(
            ID=ID,
            area=area,
            constr_assem=constr_assem,
            T_int_d=self.T_int_d,
            T_int_surf=self.T_int_surf,
            T_ext_d=self.T_ext_d,
            A_slab=A_slab,
            P_slab=P_slab,
            z=z,
            dU_tb=dU_tb,
            f_dT_an=f_dT_an,
            f_gw=f_gw
        )
        self.building_elements['ground'].append(grd_build_elem)

    @property
    def H_Tie(self) -> Quantity:
        """Heat transfer coefficient from the heated space directly to the
        exterior."""
        return (
            sum(be.H for be in self.building_elements['exterior'])
            or Q_(0.0, 'W / K')
        )

    @property
    def H_Tia(self) -> Quantity:
        """Heat transfer coefficient from the heated space to adjacent heated
        spaces."""
        return (
            sum(be.H for be in self.building_elements['adjacent_heated_space'])
            or Q_(0.0, 'W / K')
        )

    @property
    def H_TiaBE(self) -> Quantity:
        """Heat transfer coefficient from the heated space to adjacent building
        entities."""
        return (
            sum(be.H for be in self.building_elements['adjacent_building_entity'])
            or Q_(0.0, 'W / K')
        )

    @property
    def H_Tiae(self) -> Quantity:
        """Heat transfer coefficient from the heated space to adjacent unheated
        spaces or building."""
        return (
            sum(be.H for be in self.building_elements['adjacent_unheated_space'])
            or Q_(0.0, 'W / K')
        )

    @property
    def H_Tig(self) -> Quantity:
        """Heat transfer coefficient from the heated space to ground."""
        return (
            sum(be.H for be in self.building_elements['ground'])
            or Q_(0.0, 'W / K')
        )

    def get_transmission_heat_loss(self, x: str = 'total') -> Quantity:
        """Returns the design transmission heat loss from the heated space to
        another space.

        Parameters
        ----------
        x: {'total', 'exterior', 'adjacent_heated_spaces',
            'adjacent_building_entities', 'adjacent_unheated_spaces',
           'ground'}, default 'total'
           Specifies which other space is meant.

        Returns
        -------
        Quantity
            Transmission heat loss from the heated space.
        """
        dT_d = self.T_int_d - self.T_ext_d
        if x == 'total':
            H_T = self.H_Tie + self.H_Tia + self.H_TiaBE + self.H_Tiae + self.H_Tig
            return H_T * dT_d
        elif x == 'exterior':
            return self.H_Tie * dT_d
        elif x == 'adjacent_heated_spaces':
            return self.H_Tia * dT_d
        elif x == 'adjacent_building_entities':
            return self.H_TiaBE * dT_d
        elif x == 'adjacent_unheated_spaces':
            return self.H_Tiae * dT_d
        elif x == 'ground':
            return self.H_Tig * dT_d
        return None

    def get_additional_heating_up_power(self) -> Quantity:
        """Returns the additional heating-up power for the heated space."""
        Q_hu = self.area * self.q_hu
        return Q_hu

    def get_heat_load(self) -> Quantity:
        """Returns the total heating loss of the heated space."""
        Q_trm = self.get_transmission_heat_loss()
        Q_ven = self.get_ventilation_heat_loss()
        Q_hu = self.get_additional_heating_up_power()
        return Q_trm + Q_ven + Q_hu
