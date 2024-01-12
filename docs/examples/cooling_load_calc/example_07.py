"""
EXAMPLE 7
---------
CREATING A UNCONDITIONED ZONE AND SOLVING FOR THE ZONE AIR TEMPERATURE.
In an unconditioned zone the zone air temperature is not predetermined as in a
conditioned zone. Its value depends on a global heat balance of the zone
air (the sum of all heat gains to the zone air and the heat extracted from the
zone air by the cooling system should be zero).
"""
import pandas as pd

from hvac import Quantity
from hvac.sun import Location, ClimateType, ReferenceDates
from hvac.cooling_load_calc import (
    WeatherData,
    wtcb,
    ExteriorBuildingElement,
    ConstructionAssembly,
    HeatFlowDirection,
    VentilationZone
)
from hvac.cooling_load_calc import UnconditionedZone


Q_ = Quantity


class ExteriorBuildingElements:
    """Class that holds the exterior building elements of the unconditioned
    zone.
    """

    def __init__(self):
        """Creates the exterior building elements."""
        self.location = Location(
            fi=Q_(51.183, 'deg'),
            L_loc=Q_(3.8, 'deg'),
            altitude=Q_(8, 'm'),
            climate_type=ClimateType.MID_LATITUDE_SUMMER,
            timezone='Etc/GMT-1'
        )
        self.weather_data = WeatherData.create_from_climatic_design_data(
            location=self.location,
            date=ReferenceDates.get_date_for('Jul'),
            T_db_des=Q_(26.7, 'degC'),
            T_db_rng=Q_(11.3, 'K'),
            T_wb_mc=Q_(19.2, 'degC'),
            T_wb_rng=Q_(4.7, 'K')
        )
        # reference zone air temperature for creating the construction assemblies:
        self.T_zone = Q_(22, 'degC')

        self.south_wall = self._create_south_wall()
        self.west_wall = self._create_west_wall()
        self.north_wall = self._create_north_wall()
        self.east_wall = self._create_east_wall()
        self.roof = self._create_roof()

    def _create_south_wall(self):
        # We use a construction assembly factory function from the wtcb
        # subpackage to create a construction assembly with a predefined
        # configuration of the construction layers:
        constr_assem = wtcb.exterior_walls.create_ext_wall_wtcb_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.weather_data.T_db(12),
            T_int=self.T_zone,
            T_asp=self.weather_data.T_db(0),
            dT_asp=Q_(10, 'K'),
        )
        # Create the exterior wall:
        ext_wall = ExteriorBuildingElement.create(
            ID='south_wall',
            T_zone=self.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(10.6, 'm') * Q_(3.5, 'm'),
            weather_data=self.weather_data,
            gamma=Q_(0, 'deg'),
            beta=Q_(90, 'deg')
        )
        # Add a window to the exterior wall. We load the window properties
        # from the wtcb window properties shelf:
        window_props = wtcb.WindowPropertiesShelf.load(
            ID='window-5a-operable-wood/vinyl'
        )
        ext_wall.add_window(
            ID='window_south',
            width=Q_(8, 'm'),
            height=Q_(2.15, 'm'),
            props=window_props
        )
        return ext_wall

    def _create_west_wall(self):
        constr_assem = wtcb.exterior_walls.create_ext_wall_wtcb_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.weather_data.T_db(12),
            T_int=self.T_zone,
            T_asp=self.weather_data.T_db(0),
            dT_asp=Q_(10, 'K'),
        )
        ext_wall = ExteriorBuildingElement.create(
            ID='west_wall',
            T_zone=self.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(10.6, 'm') * Q_(3.5, 'm'),
            weather_data=self.weather_data,
            gamma=Q_(90, 'deg'),
            beta=Q_(90, 'deg')
        )
        ext_wall.add_door(
            ID='ext_door_west',
            width=Q_(1, 'm'),
            height=Q_(2.2, 'm'),
            constr_assem=ConstructionAssembly.create(
                ID='ext_door_assem',
                U=Q_(4.0, 'W / (K * m**2)')
            )
        )
        return ext_wall

    def _create_north_wall(self):
        constr_assem = wtcb.exterior_walls.create_ext_wall_wtcb_F1(
            t_ins=Q_(12, 'cm'),
            T_asp=self.T_zone,
            dT_asp=Q_(10, 'K'),
        )
        ext_wall = ExteriorBuildingElement.create(
            ID='north_wall',
            T_zone=self.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(10.6, 'm') * Q_(3.5, 'm'),
            weather_data=self.weather_data,
            gamma=Q_(0, 'deg'),
            beta=Q_(90, 'deg')
        )
        return ext_wall

    def _create_east_wall(self):
        constr_assem = wtcb.exterior_walls.create_ext_wall_wtcb_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.weather_data.T_db(12),
            T_int=self.T_zone,
            T_asp=self.weather_data.T_db(0),
            dT_asp=Q_(10, 'K'),
        )
        ext_wall = ExteriorBuildingElement.create(
            ID='east_wall',
            T_zone=self.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(10.6, 'm') * Q_(3.5, 'm'),
            weather_data=self.weather_data,
            gamma=Q_(-90, 'deg'),
            beta=Q_(90, 'deg')
        )
        window_props = wtcb.WindowPropertiesShelf.load(
            ID='window-5a-operable-wood/vinyl'
        )
        ext_wall.add_window(
            ID='window_east',
            width=Q_(3, 'm'),
            height=Q_(2.15, 'm'),
            props=window_props
        )
        return ext_wall

    def _create_roof(self):
        constr_assem = wtcb.roofs.create_roof_wtcb_F1(
            t_ins=Q_(12, 'cm'),
            heat_flow_dir=HeatFlowDirection.DOWNWARDS,
            T_ext=self.weather_data.T_db(12),
            T_int=self.T_zone,
            T_asp=self.weather_data.T_db(0),
            dT_asp=Q_(10, 'K')
        )
        roof = ExteriorBuildingElement.create(
            ID='roof',
            T_zone=self.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(10.6, 'm') * Q_(10.6, 'm'),
            weather_data=self.weather_data,
            gamma=Q_(0, 'deg'),
            beta=Q_(0, 'deg')
        )
        return roof


def main():
    # Get the exterior building elements of the zone:
    ext_build_elems = ExteriorBuildingElements()

    # Create a ventilation zone where the thermal zone is part of:
    vez = VentilationZone.create(ID='office-ventilation')

    # Create the unconditioned zone:
    unconditioned_zone = UnconditionedZone.create(
        ID='unconditioned-office',
        weather_data=ext_build_elems.weather_data,
        floor_area=Q_(100, 'm**2'),
        height=Q_(3.5, 'm'),
        ventilation_zone=vez,
        C_tsn=Q_(200, 'kJ / (m**2 * K)'),
        A_tsn=Q_(100, 'm**2'),
        R_tsn=Q_(0.13, 'K * m**2 / W')
    )

    # Add the exterior building elements to the zone:
    unconditioned_zone.add_ext_build_elem([
        ext_build_elems.south_wall,
        ext_build_elems.west_wall,
        ext_build_elems.north_wall,
        ext_build_elems.east_wall,
        ext_build_elems.roof
    ])

    # Add space ventilation to the unconditioned zone:
    unconditioned_zone.add_ventilation()

    # Solve the unconditioned zone for the zone air temperature:
    unconditioned_zone.solve(
        F_rad=Q_(0.46, 'frac'),
        Q_dot_sys_fun=lambda t_sol_sec: Q_(0, 'W'),
        dt_hr=1.0,
        num_cycles=25
    )

    df = unconditioned_zone.temperature_heat_gain_table()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print(df)


if __name__ == '__main__':
    main()
