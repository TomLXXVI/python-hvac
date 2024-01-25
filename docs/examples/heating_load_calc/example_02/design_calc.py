"""EXAMPLE 02 (part 1)
----------------------
HEATING LOAD DESIGN CALCULATION

In this example the heating load under design conditions is calculated for a
very simple, small single-zone building. Based on this heat load, the heating
system can be sized. In part 2 of this example, we will then simulate this
building together with its heating system on a winter day. The heating system
will consist of panel radiators. With the simulation we can analyze the effect
of the water volume flow rate through the radiators and of the water supply
temperature on the space air temperature.
"""
from hvac import Quantity
from hvac.heating_load_calc import (
    ClimateDesignData,
    Building,
    ConstructionAssembly,
    HeatFlowDirection,
)
from hvac.cooling_load_calc import wtcb


Q_ = Quantity


def main():
    """Creates the model of a simplified single-zone building and prints it
    heating load on the peak-winter design day.
    """
    # Specify the climatic data that applies to the geographic location where the
    # building is situated for the peak-winter design day:
    climatic_design_data = ClimateDesignData(
        T_ext_d=Q_(-5.7, 'degC'),   # outdoor temperature on the peak-winter design day
        T_ext_an=Q_(11.5, 'degC'),  # annual average outdoor temperature
        T_ext_min=Q_(4.4, 'degC')   # average of the daily minimum temperatures in the coldest month
    )

    # Set the desired indoor temperature:
    T_int_des = Q_(20, 'degC')

    # Create the model of our building:
    building = SomeSingleZoneBuilding(climatic_design_data, T_int_des)

    # Print the heating load of this building:
    print(building.get_heat_load())


class ConstructionAssemblyFactory:
    """Contains functions to create the construction assemblies of the
    building elements of our building. It uses the `wtcb` subpackage to quickly
    configure construction assemblies that have a predefined composition
    of the construction layers.
    """
    def __init__(
        self,
        climatic_design_data: ClimateDesignData,
        T_int_des: Quantity
    ) -> None:
        self.climatic_design_data = climatic_design_data
        self.T_int_des = T_int_des

    def create_exterior_wall(self) -> ConstructionAssembly:
        ew = wtcb.exterior_walls.create_ext_wall_wtcb_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.climatic_design_data.T_ext_min,
            T_int=self.T_int_des
        )
        return ew

    def create_roof(self) -> ConstructionAssembly:
        rf = wtcb.roofs.create_roof_wtcb_F1(
            t_ins=Q_(18, 'cm'),
            heat_flow_dir=HeatFlowDirection.UPWARDS,
            T_ext=self.climatic_design_data.T_ext_d,
            T_int=self.T_int_des
        )
        return rf

    def create_floor(self) -> ConstructionAssembly:
        fl = wtcb.floors.create_floor_wtcb_F3(
            t_ins=Q_(8, 'cm'),
            heat_flow_dir=HeatFlowDirection.DOWNWARDS,
            T_ext=self.climatic_design_data.T_ext_min,
            T_int=self.T_int_des
        )
        return fl


class SomeSingleZoneBuilding:
    """Models a very simple single-zone building with an interior floor area
    of 100 m² and a ceiling height of 3 m. The building has four exterior walls,
    a roof, and a floor. All the walls have the same interior dimensions:
    length is 10 m and height is equal to the ceiling height. The building has
    a skylight on the roof.
    """
    def __init__(
        self,
        climatic_design_data: ClimateDesignData,
        T_int_des: Quantity
    ) -> None:
        """Creates the single-zone building."""
        # Create the building hierarchy:
        self.building = Building.create(
            ID='my-test-building',
            climate_data=climatic_design_data
        )
        self.building_entity = self.building.add_building_entity(
            ID='my-test-building-entity',
        )
        self.ventilation_zone = self.building_entity.add_ventilation_zone(
            ID='my-test-ventilation-zone'
        )
        self.space = self.ventilation_zone.add_heated_space(
            ID='my-test-space',
            height=Q_(3, 'm'),
            area=Q_(100, 'm**2')
        )

        # Instantiate the `ConstructionAssemblyFactory` to create the
        # construction assemblies:
        self.constr_assem_factory = ConstructionAssemblyFactory(
            climatic_design_data=climatic_design_data,
            T_int_des=T_int_des
        )
        constr_assem_ext_wall = self.constr_assem_factory.create_exterior_wall()
        constr_assem_roof = self.constr_assem_factory.create_roof()
        constr_assem_floor = self.constr_assem_factory.create_floor()

        # Determine the exterior dimensions of the exterior building elements:
        ext_wall_thick = constr_assem_ext_wall.thickness
        roof_thick = constr_assem_roof.thickness
        L_ext_wall = Q_(10, 'm') + 2 * ext_wall_thick
        H_ext_wall = self.space.height + roof_thick
        A_ext_wall = L_ext_wall * H_ext_wall

        # Add the exterior building elements to the space:
        for ID in ['S-wall', 'N-wall', 'E-wall', 'W-wall']:
            self.space.add_exterior_building_element(
                ID=ID,
                area=A_ext_wall,
                constr_assem=constr_assem_ext_wall
            )
        self.roof = self.space.add_exterior_building_element(
            ID='roof',
            area=L_ext_wall**2,
            constr_assem=constr_assem_roof
        )
        self.space.add_ground_building_element(
            ID='floor',
            area=self.space.area,
            constr_assem=constr_assem_floor,
            A_slab=L_ext_wall**2,
            P_slab=L_ext_wall * 4,
            z=Q_(0, 'm')
        )

        self.roof.add_building_element(
            ID='skylight',
            area=Q_(25, 'm**2'),
            constr_assem=wtcb.WindowPropertiesShelf.load('window-5a-operable-wood/vinyl'),
        )

    def get_heat_load(self) -> Quantity:
        """Returns the heating load of the single-zone building."""
        return self.space.get_heat_load()


if __name__ == '__main__':
    main()
