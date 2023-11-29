"""
COOLING LOAD CALCULATION OF AN EXPOSITION HALL
----------------------------------------------

The exposition hall has three exterior walls (north, west, and east). The
roof has two skylights. The exterior walls have no windows. Heat flows
through interior partition walls are assumed zero. The building can be
considered as one single space. The building's floor area is 18,68 m x 28,02 m
= 523,4136 m². The building height is 6 m.

(More information about the cooling load calculation API can also be found in
the file "info.md" in docs/examples/heating_load_calc.)
"""
from datetime import date
import pandas as pd
from hvac import Quantity
from hvac.charts import LineChart
from hvac.cooling_load_calc import (
    Geometry,
    HeatFlowDirection,
    Material,
    SurfaceLayer,
    BuildingComponent,
    AirSpace,
    ConstructionAssembly,
    WindowThermalProperties,
    Space,
    VentilationZone,
    BuildingEntity,
    Building,
    Location,
    ClimateData,
    TemperatureSchedule,
    OnOffSchedule,
    OccupancySchedule,
    InternalHeatGain,
    PeopleHeatGain
)

Q_ = Quantity


# First, we will define the different construction assemblies the building
# envelope is (or can be) made of. The construction assemblies are assembled
# inside functions that we group inside a class `ConstructionAssemblies`.
# Through the function parameters, it will be possible to adapt our construction
# assemblies to different conditions, e.g., we can still choose the insulation
# thickness inside a wall when we define later on the construction of the building
# envelope.


class ConstructionAssemblies:
    """
    Class that groups the construction assemblies from which the elements of
    the building envelope are made.
    """
    @staticmethod
    def exterior_wall_light(
        T_ext: Quantity,
        T_int: Quantity,
        T_asp: Quantity,
        dT_asp: Quantity,
        v_wind: Quantity,
        t_ins: Quantity
    ) -> ConstructionAssembly:
        """
        Creates the construction assembly of a light-weight exterior wall.

        Parameters
        ----------
        T_ext:
            outdoor air temperature
        T_int:
            indoor air temperature
        T_asp:
            mean airspace temperature
        dT_asp:
            temperature difference across airspace
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
        ext_surf_film = SurfaceLayer.create(
            ID='ext_surf_film',
            geometry=Geometry(),
            heat_flow_direction=HeatFlowDirection.HORIZONTAL,
            Tmn=T_ext,
            internal_surface=False,
            wind_speed=v_wind
        )
        cement_plaster = BuildingComponent.create(
            ID='cement_plaster',
            geometry=Geometry(t=Q_(10, 'mm')),
            material=Material(
                k=Q_(0.950, 'W / (m * K)'),
                rho=Q_(1900, 'kg / m ** 3'),
                c=Q_(840, 'J / (kg * K)')
            )
        )
        chipboard = BuildingComponent.create(
            ID='chipboard',
            geometry=Geometry(t=Q_(25, 'mm')),
            material=Material(
                k=Q_(0.100, 'W / (m * K)'),
                rho=Q_(450, 'kg / m ** 3'),
                c=Q_(1880, 'J / (kg * K)')
            )
        )
        air_gap = AirSpace.create(
            ID='air_gap',
            geometry=Geometry(t=Q_(40, 'mm')),
            dT=dT_asp,
            heat_flow_direction=HeatFlowDirection.HORIZONTAL,
            Tmn=T_asp
        )
        insulation = BuildingComponent.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=Material(
                k=Q_(0.035, 'W / (m * K)'),
                rho=Q_(35, 'kg / m ** 3'),
                c=Q_(840, 'J / (kg * K)')
            )
        )
        insulation.slices = 7
        # The insulation layer is divided in a number of slices from which the
        # linear thermal network model of the exterior wall will be derived.
        # (Each slice will represent a temperature node in the linear thermal
        # network model of the wall. By default, a layer has always only one
        # temperature node in which all the heat capacity of the layer is
        # thought to be concentrated.)
        plasterboard = BuildingComponent.create(
            ID='plasterboard',
            geometry=Geometry(t=Q_(13, 'mm')),
            material=Material(
                k=Q_(0.230, 'W / (m * K)'),
                rho=Q_(900, 'kg / m ** 3'),
                c=Q_(840, 'J / (kg * K)')
            )
        )
        int_surf_film = SurfaceLayer.create(
            ID='int_surf_film',
            geometry=Geometry(),
            heat_flow_direction=HeatFlowDirection.HORIZONTAL,
            Tmn=T_int
        )
        # Create the construction assembly of the exterior wall:
        ext_wall = ConstructionAssembly.create(
            ID=f"ext_wall_type1 (t_ins = {t_ins.to('mm'):~P.0f})",
            # Layers must be ordered from the outdoor to the indoor:
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
        """
        Creates the construction assembly of a heavy-weight exterior wall.

        Parameters
        ----------
        See `create_exterior_wall_light`.

        Returns
        -------
        ConstructionAssembly
        """
        # Create the layers from outdoors to indoors that compose the
        # construction assembly:
        ext_surf_film = SurfaceLayer.create(
            ID='ext_surf_film',
            geometry=Geometry(),
            heat_flow_direction=HeatFlowDirection.HORIZONTAL,
            Tmn=T_ext,
            internal_surface=False,
            wind_speed=v_wind
        )
        concrete_wall = BuildingComponent.create(
            ID='concrete_wall',
            geometry=Geometry(t=Q_(14, 'cm')),
            material=Material(
                k=Q_(1.19, 'W / (m * K)'),
                rho=Q_(1700, 'kg / m ** 3'),
                c=Q_(1000, 'J / (kg * K)')
            )
        )
        concrete_wall.slices = 7
        insulation = BuildingComponent.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=Material(
                k=Q_(0.04, 'W / (m * K)'),
                rho=Q_(22.0, 'kg / m ** 3'),
                c=Q_(1030, 'J / (kg * K)')
            )
        )
        insulation.slices = 7
        int_surf_film = SurfaceLayer.create(
            ID='int_surf_film',
            geometry=Geometry(),
            heat_flow_direction=HeatFlowDirection.HORIZONTAL,
            Tmn=T_int
        )
        # Create the construction assembly of the exterior wall:
        ext_wall = ConstructionAssembly.create(
            ID=f"ext_wall_type1 (t_ins = {t_ins.to('mm'):~P.0f})",
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
        """
        Creates the construction assembly of the roof.

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
        ext_surf_film = SurfaceLayer.create(
            ID='ext_surf_film',
            geometry=Geometry(),
            heat_flow_direction=heat_flow_dir,
            Tmn=T_ext,
            internal_surface=False,
            wind_speed=v_wind
        )
        roofing_felt = BuildingComponent.create(
            ID='roofing_felt',
            geometry=Geometry(t=Q_(3, 'mm')),
            material=Material(
                k=Q_(0.170, 'W / (m * K)'),
                rho=Q_(1200, 'kg / m ** 3'),
                c=Q_(1470, 'J / (kg * K)')
            )
        )
        roofing_felt.slices = 3
        insulation = BuildingComponent.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=Material(
                k=Q_(0.035, 'W / (m * K)'),
                rho=Q_(15, 'kg / m ** 3'),
                c=Q_(1470, 'J / (kg * K)')
            )
        )
        insulation.slices = 8
        concrete = BuildingComponent.create(
            ID='concrete',
            geometry=Geometry(t=Q_(100, 'mm')),
            material=Material(
                k=Q_(1.9, 'W / (m * K)'),
                rho=Q_(2500, 'kg / m ** 3'),
                c=Q_(840, 'J / (kg * K)')
            )
        )
        concrete.slices = 10
        int_surf_film = SurfaceLayer.create(
            ID='int_surf_film',
            geometry=Geometry(),
            heat_flow_direction=heat_flow_dir,
            Tmn=T_int
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
        """
        Returns the thermal and radiative properties of the skylights on the
        roof.
        """
        wnd_props = WindowThermalProperties(
            ID='sky_light',
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
# building. For this, we write the class `BuildingConstructor`.

class BuildingConstructor:
    """
    Class for creating the building construction and put it around the
    single-zone building.
    """
    # The construction assemblies of the building envelope:
    light_wall: ConstructionAssembly
    heavy_wall: ConstructionAssembly
    roof: ConstructionAssembly

    @classmethod
    def construct(cls, space: Space, is_heavy: bool = False) -> None:
        """
        Constructs the building envelope around the given space.

        With boolean `is_heavy`, we can still make a choice between either a
        heavy-weight (concrete exterior walls), or a light-weight wall
        construction.
        """
        # First, we need to instantiate the construction assemblies of the
        # exterior walls and roof of the building. For this, we can take the
        # indoor air and outdoor air temperature from our `Space` object. The
        # other parameters needed to instantiate the construction assemblies are
        # selected inside the method `_create_construction_assemblies`.
        # We define three construction assemblies: (1) for a light-weight
        # exterior wall, (2) for a heavy-weight exterior wall, and (3) for the
        # roof.
        cls._create_construction_assemblies(
            T_int=space.T_int_fun.base_value,
            T_ext=space.climate_data.Tdb_avg
        )
        # Depending on the choice made between a light-weight or heavy-weight
        # wall, we select the appropriate exterior wall construction assembly:
        if is_heavy:
            wall = cls.heavy_wall
        else:
            wall = cls.light_wall
        # Now we can construct all the building elements of the building
        # envelope:
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
        """
        Creates the construction assemblies of the space.

        Parameters
        ----------
        T_int:
            Interior temperature.
        T_ext:
            Exterior temperature.
        """
        cls.light_wall = ConstructionAssemblies.exterior_wall_light(
            T_ext=T_ext,
            T_int=T_int,
            T_asp=Q_(32.2, 'degC'),
            dT_asp=Q_(5.6, 'K'),
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
        space: Space,
        constr_assem: ConstructionAssembly
    ) -> None:
        """Adds the north wall to the building envelope of the space."""
        space.add_ext_building_element(
            ID='north_wall',
            azimuth=Q_(0, 'deg'),
            tilt=Q_(90.0, 'deg'),
            width=space.width,
            height=space.height,
            construction_assembly=constr_assem,
            surface_color='light-colored'
        )

    @staticmethod
    def _construct_west_wall(
        space: Space,
        constr_assem: ConstructionAssembly
    ) -> None:
        """Adds the west wall to the building envelope of the space."""
        space.add_ext_building_element(
            ID='west_wall',
            azimuth=Q_(270, 'deg'),
            tilt=Q_(90, 'deg'),
            width=space.length,
            height=space.height,
            construction_assembly=constr_assem,
            surface_color='light-colored'
        )

    @staticmethod
    def _construct_east_wall(
        space: Space,
        constr_assem: ConstructionAssembly
    ) -> None:
        """Adds the east wall to the building envelope of the space."""
        space.add_ext_building_element(
            ID='east_wall',
            azimuth=Q_(90, 'deg'),
            tilt=Q_(90, 'deg'),
            width=space.length,
            height=space.width,
            construction_assembly=constr_assem,
            surface_color='light-colored'
        )

    @staticmethod
    def _construct_roof(
        space: Space,
        constr_assem: ConstructionAssembly
    ) -> None:
        """Adds the roof to the building envelope of the space."""
        roof = space.add_ext_building_element(
            ID='roof',
            azimuth=Q_(0, 'deg'),
            tilt=Q_(0, 'deg'),
            width=space.width,
            height=space.length,
            construction_assembly=constr_assem,
            surface_color='dark-colored'
        )
        # Add skylight 1:
        roof.add_window(
            ID='sky_light_1',
            width=Q_(3, 'm'),
            height=Q_(10, 'm'),
            therm_props=ConstructionAssemblies.sky_light(),
        )
        # Add skylight 2:
        roof.add_window(
            ID='sky_light_2',
            width=Q_(3, 'm'),
            height=Q_(10, 'm'),
            therm_props=ConstructionAssemblies.sky_light(),
        )

    @staticmethod
    def _add_interior_thermal_mass(space: Space) -> None:
        """Adds interior thermal mass to the space."""
        space.add_internal_thermal_mass(
            A=space.floor_area,
            R=Q_(0.015, 'm ** 2 * K / W'),
            C=Q_(300, 'kJ / (m ** 2 * K)')  # heavy-weight construction
        )


# Using our `BuildingConstructor` class, we can now create the thermal model of
# our single-zone building, that we'll call `ExpositionHall`.

class ExpositionHall:
    """Class for modeling the single-zone exposition hall."""
    @classmethod
    def create(
        cls,
        climate: ClimateData,
        T_comfort: Quantity = Q_(26, 'degC'),
        T_economy: Quantity = Q_(26, 'degC'),
        num_people_max: int = 200,
        num_people_min: int = 100,
        is_heavy_construction: bool = False
    ) -> Space:
        """
        Creates the thermal model of our single-zone exposition hall and
        returns its single space.

        Parameters
        ----------
        climate:
            The climate data applicable to the peak summer design-day (or any
            other day of the year).
        T_comfort:
            The zone air setpoint temperature when "comfort cooling" is active.
        T_economy:
            The zone air setpoint temperature when "economy cooling" is active
            (this setpoint temperature should be equal to our higher than the
            comfort value).
        num_people_max:
            The number of people that may be present in the exposition hall
            during daytime hours (from 8 until 18 o'clock).
        num_people_min:
            The number of people that may be present in the exposition hall
            during nighttime hours (from 18 until 0 o'clock and from 0 until
            8 o'clock).
        is_heavy_construction:
            Indicate if the exterior walls of the exposition hall have a
            heavy-weight construction (True) or a light-weight construction
            (False).
        """
        # Create a time schedule for turning the cooling system on/off:
        cooling_schedule = None  # not used - cooling is permanently active
        # Set a time schedule to change the setpoint of the interior
        # temperature between "comfort" and "economy":
        int_temp_schedule = cls._set_T_setpoint_schedule(
            T_comfort=T_comfort,
            T_economy=T_economy
        )
        # Add the internal heat gains inside the building:
        internal_heat_gains = cls._add_internal_heat_gains(
            people_max=num_people_max,
            people_min=num_people_min
        )
        # Create the thermal model of the single-zone building:
        space = cls._create_single_zone_building(
            climate=climate,
            interior_temperature_schedule=int_temp_schedule,
            cooling_schedule=cooling_schedule,
            internal_heat_gains=internal_heat_gains,
            is_heavy_construction=is_heavy_construction
        )
        return space

    @staticmethod
    def _create_single_zone_building(
        climate: ClimateData,
        interior_temperature_schedule: TemperatureSchedule,
        cooling_schedule: OnOffSchedule | None,
        internal_heat_gains: list[InternalHeatGain],
        is_heavy_construction: bool = False
    ) -> Space:
        """
        Creates and configures the building hierarchy and returns the single
        space of the building.

        Parameters
        ----------
        climate:
            The climate data applicable to the peak summer design-day (or any
            other day of the year).
        interior_temperature_schedule:
            Time schedule to change the setpoint of the interior space air
            temperature between "comfort" and "economy" (see method
            _set_T_setpoint_schedule)
        cooling_schedule:
            Time schedule for turning the cooling system on or off. If None,
            the cooling system is permanently active.
        internal_heat_gains:
            List with the internal heat gains in the single-zone building
            (see method _add_internal_heat_gains).
        """
        # Create the space:
        space = Space.create(
            ID='exposition_hall',
            height=Q_(6.0, 'm'),
            length=Q_(1.5 * 18.68, 'm'),
            width=Q_(18.68, 'm'),
            climate_data=climate,
            T_int_fun=interior_temperature_schedule,
            cooling_schedule=cooling_schedule,
        )
        # Add the building envelope around the space:
        BuildingConstructor.construct(space, is_heavy=is_heavy_construction)

        # Add the internal heat gains to the space:
        for ihg in internal_heat_gains:
            space.add_internal_heat_gain(ihg)

        # Define the building hierarchy of the single-zone building:
        # Create a ventilation zone, and add the space to this zone:
        vz = VentilationZone.create(ID='exposition_hall')
        vz.add_space(space)
        # Configure the ventilation characteristics of the space (keep all the
        # default values):
        space.configure_ventilation()
        # Create a building entity, and add the ventilation zone to this
        # building entity:
        be = BuildingEntity.create(ID='exposition_hall')
        be.add_ventilation_zone(vz)
        # Create a building, and add the building entity to this building:
        bu = Building.create(ID='exposition_hall')
        bu.add_building_entity(be)
        # Only return the single space of the building:
        return space
   
    @staticmethod
    def _set_T_setpoint_schedule(
        T_comfort: Quantity,
        T_economy: Quantity
    ) -> TemperatureSchedule:
        """
        Creates a time schedule for changing the setpoint of the space air
        temperature.

        We assume that between 18 h and 8 h, the space air temperature setpoint
        is set to the "economy" value. During the other hours of the day, the
        setpoint is set to the "comfort" value.
        """
        # Create a temperature schedule with a base value of `T_comfort`:
        schedule = TemperatureSchedule.create(
            ID='sch_space_air',
            base_value=T_comfort
        )

        # Change the setpoint to `T_economy` between 0 and 8 h:
        for t in range(8):
            schedule.set_value(t, T_economy)

        # Change the setpoint again to `T_economy` between 18 and 24 h:
        for t in range(18, 24):
            schedule.set_value(t, T_economy)
        return schedule

    @staticmethod
    def _create_cooling_schedule() -> OnOffSchedule:
        """
        Creates a time schedule for turning the cooling system on or off.

        Between 18 h and 8 h the cooling system is turned off. During the other
        hours of the day, the cooling system is turned on.
        """
        schedule = OnOffSchedule.create(
            ID='cooling_schedule',
            base_value=True
        )

        for t in range(8):
            schedule.set_value(t, False)

        for t in range(18, 24):
            schedule.set_value(t, False)
        return schedule

    @staticmethod
    def _add_internal_heat_gains(
        people_max: int,
        people_min: int
    ) -> list[InternalHeatGain]:
        """Adds the internal heat gains into the space.

        Only the internal heat gain by people will be considered here.

        We assume that between 18 h and 8 h the number of people in the space
        is equal to the value of `people_min`. During the other hours of the
        day, the number of people is equal to the value of `people_max`.
        """
        occupancy_schedule = OccupancySchedule.create(
            ID='occupancy_schedule',
            base_value=people_max
        )
        for t in range(8):
            occupancy_schedule.set_value(t, people_min)
        for t in range(18, 24):
            occupancy_schedule.set_value(t, people_min)
        people = PeopleHeatGain.create(
            ID='people',
            Q_sen_person=Q_(75, 'W'),
            Q_lat_person=Q_(55, 'W'),
            F_rad=Q_(58, 'pct'),
            schedule=occupancy_schedule
        )
        return [people]


def main():
    """
    Creates the thermal model of the single-zone exposition hall and runs the
    cooling load calculation routine on this model.
    """
    # Set the geographic location and the climate design data:
    climate = ClimateData.create_from_ASHRAE_weather_data(
        location=Location(
            name='Ghent',
            lat=Q_(51.183, 'deg'),
            lon=Q_(3.8, 'deg'),
            alt=Q_(8.0, 'm'),
            tz='Europe/Brussels'
        ),
        design_day=date(2022, 7, 21),  # summer peak design-day
        # Data taken from ASHRAE's climatic design data tables
        # (ASHRAE Handbook 2017):
        Tdb_avg=Q_(26.7, 'degC'),
        Tdb_range=Q_(11.3, 'K'),
        Twb_mc=Q_(19.2, 'degC'),
        tau_beam=0.426,
        tau_dif=2.247
    )

    # Pass the climate data and set the other parameters to create our thermal
    # model of the single-zone building:
    hall = ExpositionHall.create(
        climate,
        T_comfort=Q_(26, 'degC'),
        T_economy=Q_(26, 'degC'),
        num_people_max=200,
        num_people_min=100,
        is_heavy_construction=False
    )

    # Run the cooling load calculation routine on the thermal model and get the
    # results back in a Pandas `DataFrame` object.
    Q_gains = hall.get_heat_gains(unit='kW')

    # Display a table with the heat gains and the cooling loads on an hourly
    # basis during the selected summer peak design-day:
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', 1000
    ):
        print(Q_gains)

    # Display the maximum values in the table:
    print(
        "Maximum sensible cooling load = "
        f"{Q_gains['Q_sen_load'].max():.3f} kW "
        f"({Q_gains['Q_sen_load'].idxmax()})",
        "Maximum latent cooling load = "
        f"{Q_gains['Q_lat_load'].max():.3f} kW "
        f"({Q_gains['Q_lat_load'].idxmax()})",
        "Maximum total cooling load = "
        f"{Q_gains['Q_tot_load'].max():.3f} kW "
        f"({Q_gains['Q_tot_load'].idxmax()})",
        sep='\n', end='\n\n'
    )

    # Display a diagram of the hourly cooling load and the time course of the
    # indoor air temperature:
    T_int_df = hall.get_space_air_temperatures()

    chart = LineChart()
    chart.add_y2_axis()
    chart.add_xy_data(
        label='cooling load',
        x1_values=[time.hour for time in Q_gains['time']],
        y1_values=Q_gains['Q_tot_load']
    )
    chart.add_xy_data(
        label='space air temperature',
        x1_values=[time.hour for time in Q_gains['time']],
        y2_values=T_int_df['T_space'],
        style_props={'color': 'orange'}
    )
    chart.x1.add_title('hour of the day')
    chart.x1.scale(0, 24, 1)
    chart.y1.add_title('cooling load, W')
    chart.y1.scale(0, 55, 5)
    chart.y2.add_title('space air temperature, °C')
    chart.y2.scale(0, 35, 5)
    chart.show()

    # Show a table of the heat flows into and out from the interior thermal
    # mass:
    Q_stor = hall.get_thermal_storage_heat_flows(Q_unit='kW')
    with pd.option_context('display.width', None):
        print(Q_stor)

    # Show a diagram of the hourly dry-bulb outdoor air temperature during the
    # cooling design day:
    chart2 = LineChart()
    chart2.add_xy_data(
        label='dry-bulb outdoor air temperature',
        x1_values=[time.hour for time in climate.Tdb_profile['t']],
        y1_values=[Tdb.to('degC').m for Tdb in climate.Tdb_profile['T']]
    )
    chart2.x1.add_title('hour of the day')
    chart2.x1.scale(0, 24, 1)
    chart2.y1.add_title('dry-bulb outdoor air temperature, °C')
    chart2.show()


if __name__ == '__main__':
    main()
