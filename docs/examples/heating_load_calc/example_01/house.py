"""
Example of the heating load calculation of a house.
"""
import pandas as pd
from hvac import Quantity
from hvac.heating_load_calc import (
    ClimateDesignData,
    Building,
    BuildingEntity,
    VentilationZone,
    HeatedSpace,
    UnheatedSpace,
    WindowPropertiesShelf,
    ConstructionAssembly,
    HeatFlowDirection,
    Geometry,
    SolidLayer,
    SurfaceFilm,
    MaterialShelf
)
import hvac.cooling_load_calc.construction_data.wtcb.exterior_walls as ew
import hvac.cooling_load_calc.construction_data.wtcb.interior_walls as iw
import hvac.cooling_load_calc.construction_data.wtcb.floors as fl
import hvac.cooling_load_calc.construction_data.wtcb.roofs as rf

Q_ = Quantity


class ConstructionAssemblies:
    # prepare the construction assemblies used for the building

    def __init__(self):
        # exterior walls
        self.ext_wall_wtcb_f1 = ew.create_ext_wall_F1(t_ins=Q_(12, 'cm'))

        # interior walls
        self.int_wall_wtcb_f1 = iw.create_int_wall_F1(t_ins=Q_(6, 'cm'))
        self.int_wall_proj1 = self._create_int_wall_proj1()

        # floors
        self.floor_wtcb_f4 = fl.create_floor_F4(
            t_ins=Q_(12, 'cm'),
            T_grd=Q_(10, 'degC')
        )

        # ceilings
        self.ceiling_wtcb_f13 = rf.create_ceiling_F13(t_ins=Q_(12, 'cm'))

        # windows
        self.window_ashrae_5a = WindowPropertiesShelf.load('window-5a-operable-wood/vinyl')

        # interior doors
        self.int_door = ConstructionAssembly.create(
            ID='door_int',
            U=Q_(4.0, 'W / (m ** 2 * K)'),
            geometry=Geometry(t=Q_(4, 'cm'))
        )

        # exterior doors
        self.ext_door = ConstructionAssembly.create(
            ID='door_ext',
            U=Q_(3.0, 'W / (m ** 2 * K)'),
            geometry=Geometry(t=Q_(8, 'cm'))
        )

    @staticmethod
    def _create_int_wall_proj1():
        # create the interior wall construction assembly "à la carte"
        ext_surf_film = SurfaceFilm.create(
            ID='ext_surf_film',
            geometry=Geometry(),
            heat_flow_dir=HeatFlowDirection.HORIZONTAL,
            T_mn=Q_(10, 'degC')
        )
        gypsum_layer_ext = SolidLayer.create(
            ID='gypsum_layer_ext',
            geometry=Geometry(t=Q_(1.5, 'cm')),
            material=MaterialShelf.load('gypsum-plaster')
        )
        brick_layer = SolidLayer.create(
            ID='brick_layer',
            geometry=Geometry(t=Q_(9, 'cm')),
            material=MaterialShelf.load('terracotta-block-1200kg/m3')
        )
        gypsum_layer_int = SolidLayer.create(
            ID='gypsum_layer_int',
            geometry=Geometry(t=Q_(1.5, 'cm')),
            material=MaterialShelf.load('gypsum-plaster')
        )
        int_surf_film = SurfaceFilm.create(
            ID='int_surf_film',
            geometry=Geometry(),
            heat_flow_dir=HeatFlowDirection.HORIZONTAL,
            T_mn=Q_(20, 'degC')
        )
        int_wall_proj1 = ConstructionAssembly.create(
            ID='int_wall_proj1',
            layers=[
                ext_surf_film,
                gypsum_layer_ext,
                brick_layer,
                gypsum_layer_int,
                int_surf_film
            ]
        )
        return int_wall_proj1


class House:

    def __init__(self):
        # Get the construction assemblies:
        self.constr_assem = ConstructionAssemblies()

        # Declare the building attributes:
        self.building = Building()
        self.building_entity = BuildingEntity()
        self.ventilation_zone = VentilationZone()
        self.kitchen_and_dining_room = HeatedSpace()
        self.living_room = HeatedSpace()
        self.bedroom_1 = HeatedSpace()
        self.bedroom_2 = HeatedSpace()
        self.bathroom = HeatedSpace()
        self.hall_way = UnheatedSpace()
        self.toilet = UnheatedSpace()

        # Create the building and configure all of its heated spaces:
        self._create_building()
        self._config_kitchen_and_dining_room()
        self._config_living_room()
        self._config_bedroom_1()
        self._config_bedroom_2()
        self._config_bathroom()

    def _create_building(self):
        # Create the building:
        self.building = Building.create(
            ID='house',
            climate_data=ClimateDesignData(
                T_ext_d=Q_(-7.0, 'degC'),
                T_ext_an=Q_(10.0, 'degC'),
                T_ext_min=Q_(0.0, 'degC')
            )
        )
        # Add a building entity to the building:
        self.building_entity = self.building.add_building_entity(ID='house')
        # Add a ventilation zone to the building entity:
        self.ventilation_zone = self.building_entity.add_ventilation_zone(ID='house')
        # Add the heated spaces to the ventilation zone:
        self.kitchen_and_dining_room = self.ventilation_zone.add_heated_space(
            ID='kitchen_and_dining_room',
            height=Q_(3.0, 'm'),
            area=Q_(27.0, 'm ** 2'),
            T_int_d=Q_(20, 'degC'),
            V_exh=Q_(75.0, 'm ** 3 /hr'),
            V_sup=Q_(100.0, 'm ** 3 / hr')
        )
        self.living_room = self.ventilation_zone.add_heated_space(
            ID='living_room',
            height=Q_(3.0, 'm'),
            area=Q_(14.0, 'm ** 2'),
            T_int_d=Q_(20, 'degC'),
            V_sup=Q_(75.0, 'm ** 3 / hr')
        )
        self.bedroom_1 = self.ventilation_zone.add_heated_space(
            ID='bedroom_1',
            height=Q_(3, 'm'),
            area=Q_(18, 'm ** 2'),
            T_int_d=Q_(18.0, 'degC'),
            V_sup=Q_(72, 'm ** 3 / hr')
        )
        self.bedroom_2 = self.ventilation_zone.add_heated_space(
            ID='bedroom_2',
            height=Q_(3, 'm'),
            area=Q_(11, 'm ** 2'),
            T_int_d=Q_(18.0, 'degC'),
            V_sup=Q_(72, 'm ** 3 / hr')
        )
        self.bathroom = self.ventilation_zone.add_heated_space(
            ID='bathroom',
            height=Q_(3, 'm'),
            area=Q_(7.5, 'm ** 2'),
            T_int_d=Q_(24.0, 'degC'),
            V_exh=Q_(50.0, 'm ** 3 / hr'),
            V_trf=Q_(50.0, 'm ** 3 / hr'),
            T_trf=Q_(10.0, 'degC')
        )
        self.hall_way = self.ventilation_zone.add_unheated_space(
            ID='hallway',
            height=Q_(6, 'm'),
            area=Q_(11.5, 'm ** 2'),
            T_int_d=Q_(10, 'degC'),
            V_exh=Q_(169, 'm ** 3 / hr'),
            V_trf=Q_(75, 'm ** 3 / hr')
        )
        self.toilet = self.ventilation_zone.add_unheated_space(
            ID='toilet',
            height=Q_(3, 'm'),
            area=Q_(1.75, 'm ** 2'),
            T_int_d=Q_(10, 'degC'),
            V_exh=Q_(25, 'm ** 3 / hr'),
            V_trf=Q_(25, 'm ** 3 / hr')
        )

    def _config_kitchen_and_dining_room(self) -> None:
        # Add the exterior wall at the north side of the kitchen and dining room:
        self.kitchen_and_dining_room.add_exterior_building_element(
            ID='ext_wall_north',
            area=(Q_(5.6, 'm'), Q_(3.0, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        # Add the wall adjacent to the neighboring house:
        self.kitchen_and_dining_room.add_adjacent_building_element(
            ID='adj_wall_west',
            area=(Q_(4.65, 'm'), Q_(3.0, 'm')),
            constr_assem=self.constr_assem.int_wall_wtcb_f1,
            kind_of_adjacent_space='unheated',
            T_adj=Q_(10, 'degC')
        )
        # Add the exterior wall at the east side of the kitchen and dining room:
        ext_wall_east = self.kitchen_and_dining_room.add_exterior_building_element(
            ID='ext_wall_east',
            area=(Q_(6.225, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        # Add the window and the backdoor to the east wall:
        ext_wall_east.add_building_element(
            ID='window_east',
            area=(Q_(3.7, 'm'), Q_(2.7, 'm')),
            constr_assem=self.constr_assem.window_ashrae_5a
        )
        ext_wall_east.add_building_element(
            ID='door_east',
            area=(Q_(0.8, 'm'), Q_(2.7, 'm')),
            constr_assem=self.constr_assem.ext_door
        )
        # Add the interior wall adjacent to the toilet:
        int_wall_toilet = self.kitchen_and_dining_room.add_adjacent_building_element(
            ID='int_wall_toilet',
            area=(Q_(3.025, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.int_wall_proj1,
            T_adj=Q_(10, 'degC')
        )
        # Add the door to the interior wall:
        int_wall_toilet.add_building_element(
            ID='door_toilet',
            area=(Q_(0.8, 'm'), Q_(2.1, 'm')),
            constr_assem=self.constr_assem.int_door
        )
        # Add the floor to the kitchen and dining room:
        self.kitchen_and_dining_room.add_ground_building_element(
            ID='floor',
            area=Q_(27.0, 'm ** 2'),
            constr_assem=self.constr_assem.floor_wtcb_f4,
            A_slab=Q_(59.36, 'm ** 2'),
            P_slab=Q_(32.4, 'm'),
            z=Q_(0, 'm')
        )
        # Add ceiling to the kitchen and the dining room:
        self.kitchen_and_dining_room.add_adjacent_building_element(
            ID='ceiling',
            area=Q_(27.0, 'm ** 2'),
            constr_assem=self.constr_assem.ceiling_wtcb_f13,
            kind_of_adjacent_space='heated',
            T_adj=Q_(18.60, 'degC')  # area-weighted average adjacent space temperature
        )

    def _config_living_room(self) -> None:
        ext_wall_east = self.living_room.add_exterior_building_element(
            ID='ext_wall_east',
            area=(Q_(4.375, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        ext_wall_east.add_building_element(
            ID='window_east',
            area=(Q_(3.0, 'm'), Q_(2.7, 'm')),
            constr_assem=self.constr_assem.window_ashrae_5a
        )
        ext_wall_south = self.living_room.add_exterior_building_element(
            ID='ext_wall_south',
            area=(Q_(4.0, 'm'), Q_(3.0, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        ext_wall_south.add_building_element(
            ID='window_south',
            area=(Q_(2, 'm'), Q_(1.8, 'm')),
            constr_assem=self.constr_assem.window_ashrae_5a
        )
        int_wall_hallway = self.living_room.add_adjacent_building_element(
            ID='int_wall_hallway',
            area=(Q_(4.075, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.int_wall_proj1,
            T_adj=Q_(10, 'degC')
        )
        int_wall_hallway.add_building_element(
            ID='door_hallway',
            area=(Q_(0.8, 'm'), Q_(2.1, 'm')),
            constr_assem=self.constr_assem.int_door
        )
        self.living_room.add_ground_building_element(
            ID='floor',
            area=Q_(27.0, 'm ** 2'),
            constr_assem=self.constr_assem.floor_wtcb_f4,
            A_slab=Q_(59.36, 'm ** 2'),
            P_slab=Q_(32.4, 'm'),
            z=Q_(0, 'm')
        )
        self.living_room.add_adjacent_building_element(
            ID='ceiling',
            area=Q_(14.0, 'm ** 2'),
            constr_assem=self.constr_assem.ceiling_wtcb_f13,
            kind_of_adjacent_space='heated',
            T_adj=Q_(16.67, 'degC')  # area-weighted average adjacent space temperature
        )

    def _config_bedroom_1(self) -> None:
        ext_wall_east = self.bedroom_1.add_exterior_building_element(
            ID='ext_wall_east',
            area=(Q_(5.905, 'm'), Q_(3.0, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        ext_wall_east.add_building_element(
            ID='window_east',
            area=(Q_(3.0, 'm'), Q_(1.5, 'm')),
            constr_assem=self.constr_assem.window_ashrae_5a
        )
        ext_wall_south = self.bedroom_1.add_exterior_building_element(
            ID='ext_wall_south',
            area=(Q_(5.6, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        ext_wall_south.add_building_element(
            ID='window_south',
            area=(Q_(3, 'm'), Q_(1.5, 'm')),
            constr_assem=self.constr_assem.window_ashrae_5a
        )
        self.bedroom_1.add_adjacent_building_element(
            ID='int_wall_west',
            area=(Q_(1.95, 'm'), Q_(3.0, 'm')),
            constr_assem=self.constr_assem.int_wall_wtcb_f1,
            T_adj=Q_(10, 'degC')
        )
        int_wall_hallway = self.bedroom_1.add_adjacent_building_element(
            ID='int_wall_hallway',
            area=(Q_(2.39 + 4.045, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.int_wall_proj1,
            T_adj=Q_(10, 'degC')
        )
        int_wall_hallway.add_building_element(
            ID='door_hallway',
            area=(Q_(0.8, 'm'), Q_(2.1, 'm')),
            constr_assem=self.constr_assem.int_door
        )
        self.bedroom_1.add_adjacent_building_element(
            ID='floor',
            area=Q_(18, 'm ** 2'),
            constr_assem=self.constr_assem.ceiling_wtcb_f13,
            kind_of_adjacent_space='heated',
            T_adj=Q_(20, 'degC')
        )
        self.bedroom_1.add_adjacent_building_element(
            ID='ceiling',
            area=Q_(18, 'm ** 2'),
            constr_assem=self.constr_assem.ceiling_wtcb_f13,
            T_adj=Q_(0, 'degC')
        )

    def _config_bedroom_2(self) -> None:
        ext_wall_east = self.bedroom_2.add_exterior_building_element(
            ID='ext_wall_east',
            area=(Q_(4.695, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        ext_wall_east.add_building_element(
            ID='window_east',
            area=(Q_(3.7, 'm'), Q_(1.5, 'm')),
            constr_assem=self.constr_assem.window_ashrae_5a
        )
        self.bedroom_2.add_exterior_building_element(
            ID='ext_wall_north',
            area=(Q_(2.955, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        self.bedroom_2.add_adjacent_building_element(
            ID='int_wall_bathroom',
            area=(Q_(3.35, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.int_wall_proj1,
            kind_of_adjacent_space='heated',
            T_adj=Q_(24.0, 'degC')
        )
        int_wall_hallway = self.bedroom_2.add_adjacent_building_element(
            ID='int_wall_hallway',
            area=(Q_(1.045, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.int_wall_proj1,
            T_adj=Q_(10.0, 'degC')
        )
        int_wall_hallway.add_building_element(
            ID='door_bedroom_2',
            area=(Q_(0.8, 'm'), Q_(2.1, 'm')),
            constr_assem=self.constr_assem.int_door
        )
        self.bedroom_2.add_adjacent_building_element(
            ID='floor',
            area=Q_(11, 'm ** 2'),
            constr_assem=self.constr_assem.ceiling_wtcb_f13,
            kind_of_adjacent_space='heated',
            T_adj=Q_(20, 'degC')
        )
        self.bedroom_2.add_adjacent_building_element(
            ID='ceiling',
            area=Q_(11, 'm ** 2'),
            constr_assem=self.constr_assem.ceiling_wtcb_f13,
            T_adj=Q_(0.0, 'degC')
        )

    def _config_bathroom(self):
        self.bathroom.add_adjacent_building_element(
            ID='int_wall_west',
            area=(Q_(3.65, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.int_wall_wtcb_f1,
            T_adj=Q_(10.0, 'degC')
        )
        ext_wall_north = self.bathroom.add_exterior_building_element(
            ID='ext_wall_north',
            area=(Q_(2.645, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.ext_wall_wtcb_f1
        )
        ext_wall_north.add_building_element(
            ID='window_north',
            area=(Q_(0.6, 'm'), Q_(0.6, 'm')),
            constr_assem=self.constr_assem.window_ashrae_5a
        )
        self.bathroom.add_adjacent_building_element(
            ID='int_wall_bedroom_2',
            area=(Q_(3.35, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.int_wall_proj1,
            kind_of_adjacent_space='heated',
            T_adj=Q_(18, 'degC')
        )
        int_wall_hallway = self.bathroom.add_adjacent_building_element(
            ID='int_wall_hallway',
            area=(Q_(2.3, 'm'), Q_(3, 'm')),
            constr_assem=self.constr_assem.int_wall_proj1,
            T_adj=Q_(10, 'degC')
        )
        int_wall_hallway.add_building_element(
            ID='door_hallway',
            area=(Q_(0.8, 'm'), Q_(2.1, 'm')),
            constr_assem=self.constr_assem.int_door
        )
        self.bathroom.add_adjacent_building_element(
            ID='floor',
            area=Q_(7.5, 'm ** 2'),
            constr_assem=self.constr_assem.ceiling_wtcb_f13,
            kind_of_adjacent_space='heated',
            T_adj=Q_(20, 'degC')
        )
        self.bathroom.add_adjacent_building_element(
            ID='ceiling',
            area=Q_(7.5, 'm ** 2'),
            constr_assem=self.constr_assem.ceiling_wtcb_f13,
            T_adj=Q_(0.0, 'degC')
        )


def main():
    house = House()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', 800
    ):
        print(house.ventilation_zone.get_summary())
        print()
        print(house.building_entity.get_summary())


if __name__ == '__main__':
    main()
