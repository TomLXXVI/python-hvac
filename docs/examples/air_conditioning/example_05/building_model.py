"""
COOLING LOAD CALCULATION OF A SINGLE-ZONE CONDITIONED BUILDING
---------------------------------------------------------------
In this script the cooling load calculation is programmed of a single-zone
conditioned exposition hall. The exposition hall has three exterior walls
(north, west, and east oriented). On the roof there are two sky-lights.
The exterior walls have no windows. Heat transfer through interior partition
walls can be assumed zero, as from a thermal view the building is conditioned as
a single space with a single zone air setpoint temperature.
The total floor area is 18,68 m x 28,02 m = 523,4136 m² and the building height
is 6 m.
"""
from typing import Callable
from datetime import time
import pandas as pd
from hvac import Quantity
from hvac.charts import LineChart
from hvac.sun import Location, ReferenceDates
from hvac.cooling_load_calc_old import (
    Geometry,
    HeatFlowDirection,
    Material,
    SurfaceFilm,
    SolidLayer,
    AirLayer,
    ConstructionAssembly,
    WindowThermalProperties,
    FixedTemperatureZone,
    ExteriorBuildingElement,
    VentilationZone,
    WeatherData,
    PeopleHeatGain,
)
from hvac.cooling_load_calc_old.core.utils import (
    AirLayerTemperatureSolver,
    convert_to_solar_seconds
)


Q_ = Quantity


# First, we create the different construction assemblies of the building
# envelope.

# The construction assemblies are created inside functions that we group inside
# a class `ConstructionAssemblyFactory`. Through the function parameters, it will be
# possible to adapt our construction assemblies to different conditions, e.g.,
# we can still choose the insulation thickness inside a wall when we create
# the building envelope.


class ConstructionAssemblies:
    """Class that groups the construction assemblies from which the building
    elements of the building envelope can be made.
    """
    @staticmethod
    def exterior_wall_light(
        T_ext: Quantity,
        T_int: Quantity,
        v_wind: Quantity,
        t_ins: Quantity
    ) -> ConstructionAssembly:
        """Creates the construction assembly of a light-weight exterior wall.

        Parameters
        ----------
        T_ext:
            outdoor air temperature
        T_int:
            indoor air temperature
        v_wind:
            wind speed
        t_ins:
            thickness of the insulation layer

        Returns
        -------
        ConstructionAssembly
        """
        # Create the layers from outdoors to indoors that compose the
        # construction assembly:
        ext_surf_film = SurfaceFilm.create(
            ID='ext_surf_film',
            geometry=Geometry(),
            heat_flow_dir=HeatFlowDirection.HORIZONTAL,
            T_mn=T_ext,
            is_internal_surf=False,
            wind_speed=v_wind
        )
        cement_plaster = SolidLayer.create(
            ID='cement_plaster',
            geometry=Geometry(t=Q_(10, 'mm')),
            material=Material(
                k=Q_(0.950, 'W / (m * K)'),
                rho=Q_(1900, 'kg / m ** 3'),
                c=Q_(840, 'J / (kg * K)')
            )
        )
        cement_plaster.num_slices = 2
        # The cement plaster layer is divided into 2 slices from which the
        # linear thermal network model of the exterior wall will be derived.
        # (Each slice will represent a temperature node in the linear thermal
        # network model of the wall. By default, a layer always has one
        # temperature node in which all the heat capacity of the layer is
        # thought to be concentrated.)
        chipboard = SolidLayer.create(
            ID='chipboard',
            geometry=Geometry(t=Q_(25, 'mm')),
            material=Material(
                k=Q_(0.100, 'W / (m * K)'),
                rho=Q_(450, 'kg / m ** 3'),
                c=Q_(1880, 'J / (kg * K)')
            )
        )
        chipboard.num_slices = 2
        # air gap --> see below
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=Material(
                k=Q_(0.035, 'W / (m * K)'),
                rho=Q_(35, 'kg / m ** 3'),
                c=Q_(840, 'J / (kg * K)')
            )
        )
        insulation.num_slices = int(round((t_ins.to('cm') / Q_(4, 'cm')).m))
        plasterboard = SolidLayer.create(
            ID='plasterboard',
            geometry=Geometry(t=Q_(13, 'mm')),
            material=Material(
                k=Q_(0.230, 'W / (m * K)'),
                rho=Q_(900, 'kg / m ** 3'),
                c=Q_(840, 'J / (kg * K)')
            )
        )
        plasterboard.num_slices = 2
        int_surf_film = SurfaceFilm.create(
            ID='int_surf_film',
            geometry=Geometry(),
            heat_flow_dir=HeatFlowDirection.HORIZONTAL,
            T_mn=T_int
        )
        # Create the air gap:
        # To determine the thermal resistance of the air gap the temperature
        # difference across the air gap and the mean temperature in the air
        # gap are needed. These values can be calculated with the class
        # `AirLayerTemperatureSolver`.
        air_layer_temp_solver = AirLayerTemperatureSolver(
            T_ext=T_ext,
            T_int=T_int,
            R_ea=ext_surf_film.R + cement_plaster.R + chipboard.R,
            R_ai=insulation.R + plasterboard.R + int_surf_film.R
        )
        *_, dT, T_mn = air_layer_temp_solver.solve(
            T_ae_guess=T_ext.to('K') - Q_(3, 'K'),
            T_ai_guess=T_int.to('K') + Q_(5, 'K')
        )
        air_gap = AirLayer.create(
            ID='air_gap',
            geometry=Geometry(t=Q_(40, 'mm'), w=Q_(float('inf'), 'm')),
            # we set the width of the air gap to infinity, to distinguish it
            # from an air void
            dT=dT,
            heat_flow_dir=HeatFlowDirection.HORIZONTAL,
            T_mn=T_mn
        )
        # Create the construction assembly of the exterior wall:
        ext_wall = ConstructionAssembly.create(
            ID=f"ext-wall-light-t_ins={t_ins.to('mm'):~P.0f}",
            # The layers must be ordered from the outdoor to the indoor:
            layers=[
                ext_surf_film,
                cement_plaster,
                chipboard,
                air_gap,
                insulation,
                plasterboard,
                int_surf_film
            ]
        )
        return ext_wall

    @staticmethod
    def exterior_wall_heavy(
        T_ext: Quantity,
        T_int: Quantity,
        v_wind: Quantity,
        t_ins: Quantity
    ) -> ConstructionAssembly:
        """Creates the construction assembly of a heavy-weight exterior wall.

        Parameters
        ----------
        See `create_exterior_wall_light`.

        Returns
        -------
        ConstructionAssembly
        """
        # Create the layers from outdoors to indoors that compose the
        # construction assembly:
        ext_surf_film = SurfaceFilm.create(
            ID='ext_surf_film',
            geometry=Geometry(),
            heat_flow_dir=HeatFlowDirection.HORIZONTAL,
            T_mn=T_ext,
            is_internal_surf=False,
            wind_speed=v_wind
        )
        concrete_wall = SolidLayer.create(
            ID='concrete_wall',
            geometry=Geometry(t=Q_(14, 'cm')),
            material=Material(
                k=Q_(1.19, 'W / (m * K)'),
                rho=Q_(1700, 'kg / m ** 3'),
                c=Q_(1000, 'J / (kg * K)')
            )
        )
        concrete_wall.num_slices = 7
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=Material(
                k=Q_(0.04, 'W / (m * K)'),
                rho=Q_(22.0, 'kg / m ** 3'),
                c=Q_(1030, 'J / (kg * K)')
            )
        )
        insulation.num_slices = int(round((t_ins.to('cm') / Q_(4, 'cm')).m))
        int_surf_film = SurfaceFilm.create(
            ID='int_surf_film',
            geometry=Geometry(),
            heat_flow_dir=HeatFlowDirection.HORIZONTAL,
            T_mn=T_int
        )
        # Create the construction assembly of the exterior wall:
        ext_wall = ConstructionAssembly.create(
            ID=f"ext-wall-heavy-t_ins={t_ins.to('mm'):~P.0f}",
            layers=[
                ext_surf_film,
                concrete_wall,
                insulation,
                int_surf_film
            ]
        )
        return ext_wall

    @staticmethod
    def roof(
        T_ext: Quantity,
        T_int: Quantity,
        v_wind: Quantity,
        t_ins: Quantity,
        heat_flow_dir: HeatFlowDirection
    ) -> ConstructionAssembly:
        """Creates the construction assembly of the roof.

        Parameters
        ----------
        `heat_flow_dir`:
            Indication of the direction of heat flow through the construction
            assembly.

        Other Parameters
        ----------------
        See `create_exterior_wall_light`.
        """
        # Create the layers from outdoors to indoors that compose the
        # construction assembly:
        ext_surf_film = SurfaceFilm.create(
            ID='ext_surf_film',
            geometry=Geometry(),
            heat_flow_dir=heat_flow_dir,
            T_mn=T_ext,
            is_internal_surf=False,
            wind_speed=v_wind
        )
        roofing_felt = SolidLayer.create(
            ID='roofing_felt',
            geometry=Geometry(t=Q_(3, 'mm')),
            material=Material(
                k=Q_(0.170, 'W / (m * K)'),
                rho=Q_(1200, 'kg / m ** 3'),
                c=Q_(1470, 'J / (kg * K)')
            )
        )
        roofing_felt.num_slices = 2
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=Material(
                k=Q_(0.035, 'W / (m * K)'),
                rho=Q_(15, 'kg / m ** 3'),
                c=Q_(1470, 'J / (kg * K)')
            )
        )
        insulation.num_slices = int(round((t_ins.to('cm') / Q_(4, 'cm')).m))
        concrete = SolidLayer.create(
            ID='concrete',
            geometry=Geometry(t=Q_(100, 'mm')),
            material=Material(
                k=Q_(1.9, 'W / (m * K)'),
                rho=Q_(2500, 'kg / m ** 3'),
                c=Q_(840, 'J / (kg * K)')
            )
        )
        concrete.num_slices = 5
        int_surf_film = SurfaceFilm.create(
            ID='int_surf_film',
            geometry=Geometry(),
            heat_flow_dir=heat_flow_dir,
            T_mn=T_int
        )
        # Create the construction assembly of the roof:
        roof = ConstructionAssembly.create(
            ID='roof',
            layers=[
                ext_surf_film,
                roofing_felt,
                insulation,
                concrete,
                int_surf_film
            ]
        )
        return roof

    @staticmethod
    def sky_light() -> WindowThermalProperties:
        """Returns the thermal conductive and radiative properties of the
        skylights on the roof.
        """
        wnd_props = WindowThermalProperties(
            ID='sky-light',
            U=Q_(1.6, 'W / (m ** 2 * K)'),
            SHGC_cog_dir={
                0.0: 0.4,
                40.0: 0.4,
                50.0: 0.4,
                60.0: 0.4,
                70.0: 0.4,
                80.0: 0.4
            },
            SHGC_cog_dif=0.4,
            SHGC_wnd=0.4
        )
        return wnd_props


# After the construction assemblies have been defined, we can construct the
# building elements. This will be done with the class `BuildingConstructor`.

class BuildingConstructor:
    """Class for creating the building envelope and put it around the
    single-zone building.
    """
    # The construction assemblies of the building envelope:
    light_wall: ConstructionAssembly
    heavy_wall: ConstructionAssembly
    roof: ConstructionAssembly

    @classmethod
    def construct(cls, space: FixedTemperatureZone, is_heavy: bool = False) -> None:
        """Takes the given space and adds the building envelope around it.
        With the boolean parameter `is_heavy`, we can still make a choice
        between either a heavy-weight (concrete exterior walls), or a
        light-weight wall construction.
        """
        # First, we need to instantiate the construction assemblies of the
        # exterior walls and the roof of the building. For this, we can take the
        # design indoor zone air temperature from our `ConditionedZone`
        # object and the average dry-bulb temperature of the outdoor air on the
        # design day. The other parameters needed to instantiate the construction
        # assemblies are selected inside the method `_create_construction_assemblies`.
        # We define three construction assemblies: (1) for a light-weight
        # exterior wall, (2) for a heavy-weight exterior wall, and (3) for the
        # roof.
        cls._create_construction_assemblies(
            T_int=space.T_zone_des,
            T_ext=space.weather_data.T_db_avg
        )
        # Depending on the choice made between a light-weight or heavy-weight
        # wall, we select the appropriate exterior wall construction assembly:
        if is_heavy:
            wall = cls.heavy_wall
        else:
            wall = cls.light_wall
        # Now we can construct all the exterior building elements of the building
        # envelope and add them to the space:
        cls._construct_north_wall(space, wall)
        cls._construct_west_wall(space, wall)
        cls._construct_east_wall(space, wall)
        cls._construct_roof(space, cls.roof)
        # The interior thermal mass of the building will affect the cooling
        # load of the building at a given hour of the day:
        cls._add_interior_thermal_mass(space)

    @classmethod
    def _create_construction_assemblies(
        cls,
        T_int: Quantity,
        T_ext: Quantity
    ) -> None:
        """Creates the construction assemblies of the space.

        Parameters
        ----------
        T_int:
            Design value of the interior zone air temperature.
        T_ext:
            Design value of the outdoor air temperature.
        """
        cls.light_wall = ConstructionAssemblies.exterior_wall_light(
            T_ext=T_ext,
            T_int=T_int,
            v_wind=Q_(3.4, 'm / s'),
            t_ins=Q_(140, 'mm')
        )
        cls.heavy_wall = ConstructionAssemblies.exterior_wall_heavy(
            T_ext=T_ext,
            T_int=T_int,
            v_wind=Q_(3.4, 'm / s'),
            t_ins=Q_(140, 'mm')
        )
        cls.roof = ConstructionAssemblies.roof(
            T_ext=T_ext,
            T_int=T_int,
            v_wind=Q_(3.4, 'm / s'),
            t_ins=Q_(169, 'mm'),
            heat_flow_dir=HeatFlowDirection.DOWNWARDS
        )

    @staticmethod
    def _construct_north_wall(
        space: FixedTemperatureZone,
        constr_assem: ConstructionAssembly
    ) -> None:
        """Adds the north wall to the building envelope of the space."""
        wall = ExteriorBuildingElement.create(
            ID='north-wall',
            T_zone=space.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(18.68 * 7.0, 'm**2'),
            weather_data=space.weather_data,
            gamma=Q_(180, 'deg'),  # azimuth angle, North = 180°
            beta=Q_(90, 'deg'),  # surface slope angle, vertical wall = 90°
            surface_color='dark'
        )
        space.add_ext_build_elem(wall)

    @staticmethod
    def _construct_west_wall(
        space: FixedTemperatureZone,
        constr_assem: ConstructionAssembly
    ) -> None:
        """Adds the west wall to the building envelope of the space."""
        wall = ExteriorBuildingElement.create(
            ID='west-wall',
            T_zone=space.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(28.02 * 7.0, 'm**2'),
            weather_data=space.weather_data,
            gamma=Q_(90, 'deg'),  # West = +90°
            beta=Q_(90, 'deg'),
            surface_color='dark'

        )
        space.add_ext_build_elem(wall)

    @staticmethod
    def _construct_east_wall(
        space: FixedTemperatureZone,
        constr_assem: ConstructionAssembly
    ) -> None:
        """Adds the east wall to the building envelope of the space."""
        wall = ExteriorBuildingElement.create(
            ID='east-wall',
            T_zone=space.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(28.02 * 7.0, 'm**2'),
            weather_data=space.weather_data,
            gamma=Q_(-90, 'deg'),  # East = -90°
            beta=Q_(90, 'deg'),
            surface_color='dark'
        )
        space.add_ext_build_elem(wall)

    @staticmethod
    def _construct_roof(
        space: FixedTemperatureZone,
        constr_assem: ConstructionAssembly
    ) -> None:
        """Adds the roof to the building envelope of the space."""
        roof = ExteriorBuildingElement.create(
            ID='roof',
            T_zone=space.T_zone,
            constr_assem=constr_assem,
            gross_area=Q_(28.02 * 18.68, 'm**2'),
            weather_data=space.weather_data,
            gamma=Q_(180, 'deg'),
            beta=Q_(0, 'deg'),
            surface_color='dark'
        )
        space.add_ext_build_elem(roof)
        # Add skylight 1:
        roof.add_window(
            ID='sky-light-1',
            width=Q_(3, 'm'),
            height=Q_(10, 'm'),
            props=ConstructionAssemblies.sky_light(),
        )
        # Add skylight 2:
        roof.add_window(
            ID='sky-light-2',
            width=Q_(3, 'm'),
            height=Q_(10, 'm'),
            props=ConstructionAssemblies.sky_light(),
        )

    @staticmethod
    def _add_interior_thermal_mass(space: FixedTemperatureZone) -> None:
        """Adds interior thermal mass to the space."""
        space.add_thermal_storage_node(
            A=space.floor_area,
            R_tz=Q_(0.015, 'm ** 2 * K / W'),
            C=Q_(300, 'kJ / (m ** 2 * K)')  # heavy-weight construction
        )


# Using our `BuildingConstructor` class, we can now create the thermal model of
# the single-zone building.

class BuildingModeler:
    """Class for creating the thermal model of the single-zone building."""
    @classmethod
    def create(
        cls,
        weather_data: WeatherData,
        T_comfort: Quantity = Q_(26, 'degC'),
        T_economy: Quantity = Q_(26, 'degC'),
        num_people_max: int = 200,
        num_people_min: int = 100,
        is_heavy_construction: bool = False
    ) -> FixedTemperatureZone:
        """
        Creates and returns the thermal model of our single-zone exposition hall.

        Parameters
        ----------
        weather_data:
            The weather design information valid for the peak summer design-day
            (or any other day of the year).
        T_comfort:
            The zone air setpoint temperature when "comfort cooling" is active.
        T_economy:
            The zone air setpoint temperature when "economy cooling" is active
            (this setpoint temperature should be equal to our higher than the
            comfort value).
        num_people_max:
            The number of people that may be present in the exposition hall
            during daytime hours (from 09:00 until 17:00).
        num_people_min:
            The number of people that may be present in the exposition hall
            during nighttime hours.
        is_heavy_construction:
            Indicates if the exterior walls of the exposition hall have a
            heavy-weight construction (True) or a light-weight construction
            (False).
        """
        setpoint_schedule = BuildingModeler._create_setpoint_schedule(
            T_comfort, T_economy, weather_data
        )
        people_heat_gain = BuildingModeler._create_people_heat_gain(
            num_people_max, num_people_min, weather_data
        )
        space = cls._create_single_zone(
            weather_data=weather_data,
            setpoint_schedule=setpoint_schedule,
            people_heat_gain=people_heat_gain,
            is_heavy_construction=is_heavy_construction
        )
        return space
    
    @staticmethod
    def _create_single_zone(
        weather_data: WeatherData,
        setpoint_schedule: Callable[[float], Quantity],
        people_heat_gain: PeopleHeatGain,
        is_heavy_construction: bool = False
    ) -> FixedTemperatureZone:
        """
        Creates and returns the single space of the building.

        Parameters
        ----------
        weather_data:
            The weather data valid for the peak summer design-day (or any
            other day of the year).
        setpoint_schedule:
            Time schedule to change the setpoint of the interior space air
            temperature between "comfort" and "economy".
        people_heat_gain:
            The heat gain from people in the space.
        """
        # Create the ventilation zone to which the space belongs:
        vez = VentilationZone.create(
            ID='ventilation-zone'
            # keep all the default values
        )
        # Create the space:
        space = FixedTemperatureZone.create(
            ID='exposition_hall',
            weather_data=weather_data,
            T_zone=setpoint_schedule,
            floor_area=Q_(28.02 * 18.68, 'm**2'),
            height=Q_(6.0, 'm'),
            T_zone_des=Q_(26, 'degC'),
            ventilation_zone=vez
        )
        # Create the building envelope around the space:
        BuildingConstructor.construct(space, is_heavy=is_heavy_construction)
        # Add the people heat gain to the space:
        space.add_internal_heat_gain(people_heat_gain)
        return space
    
    @staticmethod
    def _create_setpoint_schedule(
        T_comfort: Quantity,
        T_economy: Quantity,
        weather_data: WeatherData
    ) -> Callable[[float], float]:
        """Returns a function with signature f(t_sol_sec: float) -> Quantity
        that takes the solar time in seconds from midnight on the design day
        and returns the setpoint of the zone air temperature at that time.
        A day is divided in two time blocks, one block from 09:00 until 17:00 
        and the other hours of the day are all in the second block. In the
        first block the comfort temperature `T_comfort` is set. Outside this
        time block the economy setpoint temperature `T_economy` is active.
        """
        # Convert 09:00 standard time and 17:00 standard time to solar seconds:
        t_sol_sec_start = convert_to_solar_seconds(
            clock_time=time(9),
            date=weather_data.date,
            L_loc=weather_data.location.L_loc,
            tz_loc='Europe/Brussels'
        )
        t_sol_sec_end = convert_to_solar_seconds(
            clock_time=time(17),
            date=weather_data.date,
            L_loc=weather_data.location.L_loc,
            tz_loc='Europe/Brussels'
        )

        def schedule(t_sol_sec: float) -> Quantity:
            if t_sol_sec_start <= t_sol_sec <= t_sol_sec_end:
                return T_comfort
            else:
                return T_economy

        return schedule

    @staticmethod
    def _create_people_heat_gain(
        people_max: int,
        people_min: int,
        weather_data: WeatherData
    ) -> PeopleHeatGain:
        """Adds the internal heat gain from people in the space.
        It is assumed that between 09:00 and 17:00 the number of people in the 
        space is equal to the value of `people_max`. During the other hours of 
        the day, the number of people is equal to the value of `people_min`.
        """
        # Convert 09:00 standard time and 17:00 standard time to solar seconds:
        def _create_occupancy_schedule() -> Callable[[float], float]:
            t_sol_sec_start = convert_to_solar_seconds(
                clock_time=time(9),
                date=weather_data.date,
                L_loc=weather_data.location.L_loc,
                tz_loc='Europe/Brussels'
            )
            t_sol_sec_end = convert_to_solar_seconds(
                clock_time=time(17),
                date=weather_data.date,
                L_loc=weather_data.location.L_loc,
                tz_loc='Europe/Brussels'
            )

            def schedule(t_sol_sec: float) -> int:
                if t_sol_sec_start <= t_sol_sec <= t_sol_sec_end:
                    return people_max
                else:
                    return people_min
            return schedule

        people_hg = PeopleHeatGain.create(
            ID='people',
            Q_dot_sen_person=Q_(70, 'W'),
            Q_dot_lat_person=Q_(45, 'W'),
            F_rad=Q_(60, 'pct'),
            schedule=_create_occupancy_schedule()
        )
        return people_hg


def main():
    """Creates the thermal model of the single-zone exposition hall and runs the
    cooling load calculation routine on this model.
    """
    # Set the geographic location and the weather design data:
    weather_data = WeatherData.create_from_climatic_design_data(
        location=Location(
            fi=Q_(51.183, 'deg'),
            L_loc=Q_(3.8, 'deg'),
            altitude=Q_(8, 'm'),
            timezone='Etc/GMT-1'
        ),
        date=ReferenceDates.get_date_for('Jul'),
        # Data taken from ASHRAE's climatic design data tables
        # (ASHRAE Handbook 2017):
        T_db_des=Q_(26.7, 'degC'),
        T_db_rng=Q_(11.3, 'K'),
        T_wb_mc=Q_(19.2, 'degC'),
        T_wb_rng=Q_(4.7, 'K')
    )

    # Pass the weather data and set the other parameters to create our thermal
    # model of the single-zone building:
    hall = BuildingModeler.create(
        weather_data,
        T_comfort=Q_(26, 'degC'),
        T_economy=Q_(32, 'degC'),
        num_people_max=200,
        num_people_min=100,
        is_heavy_construction=False
    )

    # Run the cooling load calculation routine on the thermal model and get the
    # results back in a Pandas `DataFrame` object.
    df = hall.solve(unit='kW')

    # Display a table with the heat gains and the cooling loads on an hourly
    # basis during the selected summer peak design-day:
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', 1000
    ):
        print(df)

    # Display the maximum values in the table:
    print(
        "Maximum sensible cooling load = "
        f"{df['Q_dot_sen_zone'].max():.3f} kW "
        f"({df['Q_dot_sen_zone'].idxmax()})",
        "Maximum latent cooling load = "
        f"{df['Q_dot_lat_zone'].max():.3f} kW "
        f"({df['Q_dot_lat_zone'].idxmax()})",
        "Maximum total cooling load = "
        f"{df['Q_dot_zone'].max():.3f} kW "
        f"({df['Q_dot_zone'].idxmax()})",
        sep='\n', end='\n\n'
    )

    # Display a diagram of the hourly cooling load and the time course of the
    # indoor air temperature:
    chart = LineChart()
    chart.add_xy_data(
        label='cooling load',
        x1_values=df.index,
        y1_values=df['Q_dot_zone']
    )
    chart.x1.add_title('hour of the day')
    chart.x1.scale(0, 24, 1)
    chart.y1.add_title('cooling load, kW')
    chart.y1.scale(0, 55, 5)
    chart.show()


if __name__ == '__main__':
    main()
