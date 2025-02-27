"""EXAMPLE 02 (part 1)
----------------------
HEATING LOAD DESIGN CALCULATION

In this example the heating load under design conditions is calculated for a
(very) simple, small single-zone building. Based on this heat load, the heating
system can be sized.

In part 2 of this example, we will estimate the heating energy consumption of 
this building using the bin table method with TMY data valid for the geographic 
location where the building is situated.
"""
import warnings
from dataclasses import dataclass
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from hvac.heating_load_calc import (
    ClimateDesignData,
    Building,
    ConstructionAssembly,
)
from hvac.cooling_load_calc.construction_data import wtcb, shelves


warnings.filterwarnings('ignore', category=CoolPropWarning)
Q_ = Quantity
Air = Fluid('Air')


@dataclass
class DesignInfo:
    Q_dot_tot: Quantity
    Q_dot_trm: Quantity
    Q_dot_ven: Quantity
    H_trm: Quantity
    V_dot_ven: Quantity
    T_int_des: Quantity
    T_ext_des: Quantity

    def __str__(self) -> str:
        s = (
            f"outdoor/indoor design temperature = "
            f"{self.T_ext_des.to('degC'):~P.1f}/{self.T_int_des.to('degC'):~P.1f}\n"
            f"total design heat loss = {self.Q_dot_tot.to('kW'):~P.3f}\n"
            f"transmission heat loss = {self.Q_dot_trm.to('kW'):~P.3f}\n"
            f"ventilation/infiltration heat loss = {self.Q_dot_ven.to('kW'):~P.3f}\n"
            f"transmission heat loss coefficient = {self.H_trm.to('kW / K'):~P.3f}\n"
            f"ventilation/infiltration air volume flow rate = {self.V_dot_ven.to('m**3/hr'):~P.3f}"
        )
        return s


def main():
    """
    Creates the model of a simplified single-zone building and prints it
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

    #
    print(building.get_design_load_info())


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
        ew = wtcb.exterior_walls.create_ext_wall_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.climatic_design_data.T_ext_min,
            T_int=self.T_int_des
        )
        return ew

    def create_roof(self) -> ConstructionAssembly:
        rf = wtcb.roofs.create_roof_F1(
            t_ins=Q_(18, 'cm'),
            T_ext=self.climatic_design_data.T_ext_d,
            T_int=self.T_int_des
        )
        return rf

    def create_floor(self) -> ConstructionAssembly:
        fl = wtcb.floors.create_floor_F3(
            t_ins=Q_(8, 'cm'),
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
        """
        Creates the single-zone building.
        """
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

        # Instantiate the `ConstructionAssemblyFactory` and create the
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

        # Add the skylight to the roof:
        self.roof.add_building_element(
            ID='skylight',
            area=Q_(25, 'm**2'),
            constr_assem=shelves.WindowPropertiesShelf.load('window-5a-operable-wood/vinyl'),
        )

    def get_design_load_info(self) -> DesignInfo:
        """Returns the heating load of the single-zone building."""
        delta_T_des = self.space.T_int_d - self.space.T_ext_d
        Q_dot_trm = self.space.get_transmission_heat_loss()
        Q_dot_ven = self.space.get_ventilation_heat_loss()
        Q_dot_tot = self.space.get_heat_load()
        H_trm = Q_dot_trm / delta_T_des
        H_ven = Q_dot_ven / delta_T_des
        outdoor_air = Air(T=self.space.T_ext_d, P=Q_(101.325, 'kPa'))
        V_dot_ven = H_ven / (outdoor_air.rho * outdoor_air.cp)
        return DesignInfo(
            Q_dot_tot, Q_dot_trm, Q_dot_ven,
            H_trm, V_dot_ven, self.space.T_int_d,
            self.space.T_ext_d
        )


if __name__ == '__main__':
    main()
