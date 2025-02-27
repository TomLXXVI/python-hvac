"""
EXAMPLE 6
---------
CREATING A CONDITIONED ZONE AND SOLVING FOR THE COOLING LOAD OF THE ZONE
(i.e., the rate at which heat needs to be extracted by the cooling system to
maintain the zone air temperature at its setpoint value).
"""
from typing import Callable
from datetime import datetime, time
import pytz
import pandas as pd
from hvac import Quantity
from hvac.sun import Location, ClimateType, ReferenceDates
from hvac.sun.time import solar_to_clock_time
from hvac.cooling_load_calc_old import (
    WeatherData,
    FixedTemperatureZone,
    ExteriorBuildingElement,
    InteriorBuildingElement,
    ConstructionAssembly,
    VentilationZone,
    OfficeEquipment,
    EquipmentHeatGain,
    PeopleHeatGain
)
from hvac.cooling_load_calc_old import wtcb
from hvac.cooling_load_calc_old import shelves
from hvac.cooling_load_calc_old.core.utils import convert_to_solar_seconds
from hvac.charts import LineChart

Q_ = Quantity


def main():
    # Set the geographic location:
    location = Location(
        fi=Q_(51.183, 'deg'),
        L_loc=Q_(3.8, 'deg'),
        altitude=Q_(8, 'm'),
        climate_type=ClimateType.MID_LATITUDE_SUMMER,
        timezone='Etc/GMT-1'
    )

    # Create the weather data model for the given location based on the clear-sky
    # model for solar radiation:
    weather_data = WeatherData.create_from_climatic_design_data(
        location=location,
        date=ReferenceDates.get_date_for('Jul'),  # this is the "design day"
        T_db_des=Q_(26.7, 'degC'),
        T_db_rng=Q_(11.3, 'K'),
        T_wb_mc=Q_(19.2, 'degC'),
        T_wb_rng=Q_(4.7, 'K')
    )

    # Print sunrise at solar and local time
    sr_st = datetime.combine(
        date=weather_data.date,
        time=location.sunrise,
        tzinfo=pytz.timezone(location.timezone)
    )
    sr_lt = solar_to_clock_time(sr_st, location.L_loc.to('deg').m)
    ss_st = datetime.combine(
        date=weather_data.date,
        time=location.sunset,
        tzinfo=pytz.timezone(location.timezone)
    )
    ss_lt = solar_to_clock_time(ss_st, location.L_loc.to('deg').m)

    print(
        f"design day = {weather_data.date}",
        f"sunrise at {sr_st} solar time / {sr_lt} local time",
        f"sunset at {ss_st} solar time / {ss_lt} local time",
        sep='\n'
    )

    # Create a ventilation zone (keep all the default values):
    vez = VentilationZone.create(ID='office-ventilation')

    # Create the thermal zone (an office space as example, see class `Office`
    # below):
    office = Office(
        weather_data=weather_data,
        T_zone=Q_(22, 'degC'),  # controlled setpoint temperature of zone air
        thermal_mass=(          # used to model the interior thermal mass of the zone
            Q_(300, 'kJ / (K * m**2)'),  # unit thermal capacity
            Q_(100, 'm**2'),             # surface area of interior thermal mass
            Q_(0.13, 'K * m**2 / W')     # unit thermal resistance between interior thermal mass and zone air
        ),
        vez=vez  # the ventilation zone to which the thermal zone belongs
    )

    # Display a table with the heat gains and the resulting zone load:
    df = office.heat_gains()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print(df)

    # Plot a line chart of the different heat gains and the resulting sensible
    # zone load:
    chart = LineChart()
    chart.add_xy_data(
        label='conduction HG',
        x1_values=df.index,
        y1_values=df['Q_dot_cnd'],
        style_props={'marker': 'o'}
    )
    chart.add_xy_data(
        label='solar HG',
        x1_values=df.index,
        y1_values=df['Q_dot_sol'],
        style_props={'marker': 'o'}
    )
    chart.add_xy_data(
        label='sens. internal HG',
        x1_values=df.index,
        y1_values=df['Q_dot_sen_ihg'],
        style_props={'marker': 'o'}
    )
    chart.add_xy_data(
        label='sens. vent. HG',
        x1_values=df.index,
        y1_values=df['Q_dot_sen_vent'],
        style_props={'marker': 'o'}
    )
    chart.add_xy_data(
        label='sens. zone load',
        x1_values=df.index,
        y1_values=df['Q_dot_sen_zone'],
        style_props={'marker': 'o'}
    )
    chart.add_legend()
    chart.x1.add_title('solar time index')
    chart.y1.add_title('heat flow, W')
    chart.y1.scale(-1000, 6500, 500)
    chart.show()


class Office:

    def __init__(
        self,
        weather_data: WeatherData,
        T_zone: Quantity,
        thermal_mass: tuple[Quantity, Quantity, Quantity] | None = None,
        vez: VentilationZone | None = None
    ) -> None:
        """Creates the thermal zone of the office space."""
        self.weather_data = weather_data

        # We assume the controlled zone air temperature remains constant in time
        # during the design day (the setpoint isn't modified):
        self.T_zone = lambda t: T_zone

        # Create a `ConditionedZone` object:
        self.zone = FixedTemperatureZone.create(
            ID='office',
            weather_data=weather_data,
            T_zone=self.T_zone,
            floor_area=Q_(100, 'm**2'),
            height=Q_(3, 'm'),
            ventilation_zone=vez
        )

        # Create the exterior and interior building elements that surround the
        # zone and add them to the `ConditionedZone` object:
        self.south_wall = self._create_south_wall()
        self.west_wall = self._create_west_wall()
        self.north_wall = self._create_north_wall()
        self.east_wall = self._create_east_wall()
        self.roof = self._create_roof()

        # Create the interior thermal mass of the zone:
        if thermal_mass is not None:
            self.zone.add_thermal_storage_node(
                C=thermal_mass[0],
                A=thermal_mass[1],
                R_tz=thermal_mass[2]
            )

        # Add space ventilation to the thermal zone (keep all the default values):
        if vez is not None:
            self.zone.add_ventilation()

        # Add internal heat gains to the thermal zone:
        self._add_office_equipment()
        self._add_people()

    def _create_south_wall(self):
        # We use a construction assembly factory function from the wtcb
        # subpackage to create a construction assembly with a predefined
        # configuration of the construction layers:
        constr_assem = wtcb.exterior_walls.create_ext_wall_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.weather_data.T_db(12),
            T_int=self.T_zone(0),
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
        window_props = shelves.WindowPropertiesShelf.load(
            ID='window-5a-operable-wood/vinyl'
        )
        ext_wall.add_window(
            ID='window_1',
            width=Q_(8, 'm'),
            height=Q_(2.15, 'm'),
            props=window_props
        )
        # Add the exterior wall to the `ConditionedZone` object:
        self.zone.add_ext_build_elem(ext_wall)
        return ext_wall

    def _create_west_wall(self):
        constr_assem = wtcb.exterior_walls.create_ext_wall_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.weather_data.T_db(12),
            T_int=self.T_zone(0),
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
        self.zone.add_ext_build_elem(ext_wall)
        return ext_wall

    def _create_north_wall(self):
        constr_assem = wtcb.interior_walls.create_int_wall_F1(
            t_ins=Q_(12, 'cm'),
            T_adj=self.T_zone(0),
            T_int=self.T_zone(0),
        )
        int_wall = InteriorBuildingElement.create(
            ID='north_wall',
            T_zone=self.T_zone,
            T_adj=self.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(10.6, 'm') * Q_(3.5, 'm'),
        )
        int_wall.add_door(
            ID='int_door_north',
            width=Q_(1, 'm'),
            height=Q_(2.2, 'm'),
            constr_assem=ConstructionAssembly.create(
                ID='door_assem',
                U=Q_(3.0, 'W / (m**2 * K)'),
            )
        )
        self.zone.add_int_build_elem(int_wall)
        return int_wall

    def _create_east_wall(self):
        constr_assem = wtcb.exterior_walls.create_ext_wall_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.weather_data.T_db(12),
            T_int=self.T_zone(0),
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
        window_props = shelves.WindowPropertiesShelf.load(
            ID='window-5a-operable-wood/vinyl'
        )
        ext_wall.add_window(
            ID='window_1',
            width=Q_(3, 'm'),
            height=Q_(2.15, 'm'),
            props=window_props
        )
        self.zone.add_ext_build_elem(ext_wall)
        return ext_wall

    def _create_roof(self):
        constr_assem = wtcb.roofs.create_roof_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.weather_data.T_db(12),
            T_int=self.T_zone(0),
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
        self.zone.add_ext_build_elem(roof)
        return roof

    def _add_office_equipment(self):
        """Creates an `EquipmentHeatGain` object with the internal heat gain
        from the office equipment and adds it to the thermal zone.
        """
        # Create a schedule function that returns the fraction of operating
        # office equipment as a function of time. All the office equipment
        # (100%) is running during the office opening hours from 9:00 a.m. until
        # 5:00 p.m. Outside the office opening hours, all equipment is off.
        # The schedule function to be created must take solar time in seconds
        # from midnight (0 s) and return a float between 0 and 1 to indicate
        # the fraction of equipment that is running at a particular time.
        def _create_office_schedule() -> Callable[[float], float]:
            t_sol_sec_start = convert_to_solar_seconds(
                clock_time=time(9),
                date=self.weather_data.date,
                L_loc=self.weather_data.location.L_loc,
                tz_loc='Europe/Brussels'
            )
            t_sol_sec_end = convert_to_solar_seconds(
                clock_time=time(17),
                date=self.weather_data.date,
                L_loc=self.weather_data.location.L_loc,
                tz_loc='Europe/Brussels'
            )

            def schedule(t_sol_sec: float) -> float:
                if t_sol_sec_start <= t_sol_sec <= t_sol_sec_end:
                    return 1.0
                else:
                    return 0.0
            return schedule

        # Create an `OfficeEquipment` object. The heat gain density is
        # estimated to be 7.79 W/m² (see ASHRAE Fundamentals 2017, Chapter 18,
        # 2.4 Appliances). The radiative fraction of the heat gain is estimated
        # to be 10%, which means that 90% of the heat gain is directly released
        # by convection to the zone air:
        office_equip = OfficeEquipment.create(
            ID='office-equipment',
            heat_gain_density=Q_(7.79, 'W / m**2'),
            floor_area=Q_(100, 'm**2'),
            F_rad=Q_(10, 'pct'),
            schedule=_create_office_schedule()
        )
        # Create an `EquipmentHeatGain` object and add the `OfficeEquipment`
        # object to this object:
        equip_hg = EquipmentHeatGain(ID='office-equipment')
        equip_hg.add_equipment(office_equip)
        # Note: it is also possible to add several `Equipment` objects to a
        # single `EquipmentHeatGain` object.
        # Add the `EquipmentHeatGain` object to the internal heat gains of the
        # space:
        self.zone.add_internal_heat_gain(equip_hg)

    def _add_people(self):
        """Creates a `PeopleHeatGain` object with the internal heat gain
        from the people in the office and adds it to the thermal zone.
        """
        def _create_occupancy_schedule() -> Callable[[float], float]:
            t_sol_sec_start = convert_to_solar_seconds(
                clock_time=time(9),
                date=self.weather_data.date,
                L_loc=self.weather_data.location.L_loc,
                tz_loc='Europe/Brussels'
            )
            t_sol_sec_end = convert_to_solar_seconds(
                clock_time=time(17),
                date=self.weather_data.date,
                L_loc=self.weather_data.location.L_loc,
                tz_loc='Europe/Brussels'
            )

            def schedule(t_sol_sec: float) -> int:
                if t_sol_sec_start <= t_sol_sec <= t_sol_sec_end:
                    return 9  # 9 people in the office during opening hours
                else:
                    return 0
            return schedule

        people_hg = PeopleHeatGain.create(
            ID='office-people',
            Q_dot_sen_person=Q_(70, 'W'),
            Q_dot_lat_person=Q_(45, 'W'),
            F_rad=Q_(60, 'pct'),
            schedule=_create_occupancy_schedule()
        )
        self.zone.add_internal_heat_gain(people_hg)

    def heat_gains(self) -> pd.DataFrame:
        """Get the heat gains and the resulting sensible and latent cooling
        load of the space."""
        df = self.zone.solve(
            dt_hr=1.0,
            num_cycles=6,
            unit='W'
        )
        return df


if __name__ == '__main__':
    main()
