"""Cooling load calculation of an exposition hall."""
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


class ConstructionAssemblies:

    @staticmethod
    def create_ext_wall_type1(
        T_ext: Quantity,
        T_int: Quantity,
        T_asp: Quantity,
        dT_asp: Quantity,
        v_wind: Quantity,
        t_ins: Quantity
    ) -> ConstructionAssembly:
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
        ext_wall = ConstructionAssembly.create(
            ID=f"ext_wall_type1 (t_ins = {t_ins.to('mm'):~P.0f})",
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
    def create_ext_wall_type2(
        T_ext: Quantity,
        T_int: Quantity,
        v_wind: Quantity,
        t_ins: Quantity
    ) -> ConstructionAssembly:
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
    def create_roof(
        T_ext: Quantity,
        T_int: Quantity,
        v_wind: Quantity,
        t_ins: Quantity,
        heat_flow_dir: HeatFlowDirection
    ) -> ConstructionAssembly:
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
    def create_sky_light() -> WindowThermalProperties:
        sky_light_props = WindowThermalProperties(
            ID='sky_light_props',
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
        return sky_light_props


class Hall:

    @staticmethod
    def add_north_wall(space: Space, ca: ConstructionAssembly) -> None:
        space.add_ext_building_element(
            ID='north_wall',
            azimuth=Q_(0, 'deg'),
            tilt=Q_(90.0, 'deg'),
            width=Q_(18.68, 'm'),
            height=Q_(6.0, 'm'),
            construction_assembly=ca,
            surface_color='light-colored'
        )

    @staticmethod
    def add_west_wall(space: Space, ca: ConstructionAssembly) -> None:
        space.add_ext_building_element(
            ID='west_wall',
            azimuth=Q_(270, 'deg'),
            tilt=Q_(90, 'deg'),
            width=Q_(1.5 * 18.68, 'm'),
            height=Q_(6.0, 'm'),
            construction_assembly=ca,
            surface_color='light-colored'
        )

    @staticmethod
    def add_east_wall(space: Space, ca: ConstructionAssembly) -> None:
        space.add_ext_building_element(
            ID='east_wall',
            azimuth=Q_(90, 'deg'),
            tilt=Q_(90, 'deg'),
            width=Q_(1.5 * 18.68, 'm'),
            height=Q_(6.0, 'm'),
            construction_assembly=ca,
            surface_color='light-colored'
        )

    @staticmethod
    def add_roof(space: Space, ca: ConstructionAssembly) -> None:
        roof = space.add_ext_building_element(
            ID='roof',
            azimuth=Q_(0, 'deg'),
            tilt=Q_(0, 'deg'),
            width=Q_(18.68, 'm'),
            height=Q_(1.5 * 18.68, 'm'),
            construction_assembly=ca,
            surface_color='dark-colored'
        )
        roof.add_window(
            ID='sky_light_1',
            width=Q_(3, 'm'),
            height=Q_(10, 'm'),
            therm_props=ConstructionAssemblies.create_sky_light(),
        )
        roof.add_window(
            ID='sky_light_2',
            width=Q_(3, 'm'),
            height=Q_(10, 'm'),
            therm_props=ConstructionAssemblies.create_sky_light(),
        )

    @staticmethod
    def add_int_thermal_mass(space: Space) -> None:
        space.add_internal_thermal_mass(
            A=Q_(348.81, 'm ** 2'),
            R=Q_(0.015, 'm ** 2 * K / W'),
            C=Q_(300, 'kJ / (m ** 2 * K)')
        )

    @classmethod
    def create_building(
        cls,
        climate_data: ClimateData,
        internal_heat_gains: list[InternalHeatGain] | None,
        T_setpoint_schedule: TemperatureSchedule,
        cooling_schedule: OnOffSchedule | None
    ) -> Space:
        # create space -> ventilation zone -> building entity -> building
        space = Space.create(
            ID='industrial_hall',
            height=Q_(6.0, 'm'),
            length=Q_(1.5 * 18.68, 'm'),
            width=Q_(18.68, 'm'),
            climate_data=climate_data,
            T_int_fun=T_setpoint_schedule,
            cooling_schedule=cooling_schedule
        )
        vz = VentilationZone.create(
            ID='exposition_hall'
        )
        vz.add_space(space)
        space.add_ventilation(vz)
        be = BuildingEntity.create(ID='exposition_hall')
        be.add_ventilation_zone(vz)
        bu = Building.create(ID='exposition_hall')
        bu.add_building_entity(be)

        # create construction assemblies for the building elements of the space
        # ext_wall_type1 = ConstructionAssemblies.create_ext_wall_type1(
        #     T_ext=climate_data.Tdb_avg,
        #     T_int=T_setpoint_schedule.base_value,
        #     T_asp=Q_(32.2, 'degC'),
        #     dT_asp=Q_(5.6, 'K'),
        #     v_wind=Q_(3.4, 'm / s'),
        #     t_ins=Q_(140, 'mm')
        # )
        ext_wall_type2 = ConstructionAssemblies.create_ext_wall_type2(
            T_ext=climate_data.Tdb_avg,
            T_int=T_setpoint_schedule.base_value,
            v_wind=Q_(3.4, 'm / s'),
            t_ins=Q_(140, 'mm')
        )
        roof = ConstructionAssemblies.create_roof(
            T_ext=climate_data.Tdb_avg,
            T_int=T_setpoint_schedule.base_value,
            v_wind=Q_(3.4, 'm / s'),
            t_ins=Q_(169, 'mm'),
            heat_flow_dir=HeatFlowDirection.DOWNWARDS
        )

        # create and add the building elements to the space
        ext_wall = ext_wall_type2
        cls.add_north_wall(space, ext_wall)
        cls.add_west_wall(space, ext_wall)
        cls.add_east_wall(space, ext_wall)
        cls.add_roof(space, roof)
        cls.add_int_thermal_mass(space)

        # add internal heat gains to the space
        if internal_heat_gains:
            for ihg in internal_heat_gains:
                space.add_internal_heat_gain(ihg)

        return space

    @staticmethod
    def create_T_setpoint_schedule(
        T_comfort: Quantity,
        T_absence: Quantity
    ) -> TemperatureSchedule:
        schedule = TemperatureSchedule.create(
            ID='setpoint_temperature_schedule',
            base_value=T_comfort
        )
        for t in range(8):
            schedule.set_value(t, T_absence)
        for t in range(18, 24):
            schedule.set_value(t, T_absence)
        return schedule

    @staticmethod
    def create_cooling_schedule() -> OnOffSchedule:
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
    def create_internal_heat_gains() -> list[InternalHeatGain]:
        int_heat_gains = []
        occupancy_schedule = OccupancySchedule.create(
            ID='occupancy_schedule',
            base_value=200
        )
        for t in range(8):
            occupancy_schedule.set_value(t, 100)
        for t in range(18, 24):
            occupancy_schedule.set_value(t, 100)
        people = PeopleHeatGain.create(
            ID='people',
            Q_sen_person=Q_(75, 'W'),
            Q_lat_person=Q_(55, 'W'),
            F_rad=Q_(58, 'pct'),
            schedule=occupancy_schedule
        )
        int_heat_gains.append(people)
        return int_heat_gains


def main():
    # set location and climate design data
    climate_data = ClimateData.create(
        location=Location(
            name='Ghent',
            lat=Q_(51.183, 'deg'),
            lon=Q_(3.8, 'deg'),
            alt=Q_(8.0, 'm'),
            tz='Europe/Brussels'
        ),
        design_day=date(2022, 7, 21),
        Tdb_avg=Q_(26.7, 'degC'),
        Tdb_range=Q_(11.3, 'K'),
        Twb_mc=Q_(19.2, 'degC'),
        tau_beam=0.426,
        tau_dif=2.247
    )

    # define time schedule for space air temperature setpoint
    T_setpoint_schedule = Hall.create_T_setpoint_schedule(
        T_comfort=Q_(26.0, 'degC'),
        T_absence=Q_(26.0, 'degC')
    )

    # define time schedule for turning the cooling system on/off
    # cooling_schedule = Hall.create_cooling_schedule()

    # create the internal heat gains
    ihg_list = Hall.create_internal_heat_gains()

    # create the space
    hall = Hall.create_building(
        climate_data=climate_data,
        internal_heat_gains=ihg_list,
        T_setpoint_schedule=T_setpoint_schedule,
        cooling_schedule=None
    )

    # get the cooling loads on an hourly basis
    Q_gains = hall.get_heat_gains(unit='kW')
    with pd.option_context('display.width', None):
        print(Q_gains)

    # show diagram of hourly cooling load and indoor air temperature
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
        y2_values=T_int_df['T_int'],
        style_props={'color': 'orange'}
    )
    chart.x1.add_title('hour of the day')
    chart.x1.scale(0, 24, 1)
    chart.y1.add_title('cooling load, W')
    chart.y2.add_title('space air temperature, °C')
    chart.show()

    Q_stor = hall.get_thermal_storage_heat_flows(Q_unit='kW')
    with pd.option_context('display.width', None):
        print(Q_stor)

    # show diagram of dry-bulb outdoor air temperature
    # time_arr = climate_data.Tdb_profile['t']
    # Tdb_arr = climate_data.Tdb_profile['T']
    # chart2 = LineChart()
    # chart2.add_xy_data(
    #     label='dry-bulb outdoor air temperature',
    #     x1_values=[time.hour for time in time_arr],
    #     y1_values=[Tdb.to('degC').m for Tdb in Tdb_arr]
    # )
    # chart2.x1.add_title('hour of the day')
    # chart2.x1.scale(0, 24, 1)
    # chart2.y1.add_title('dry-bulb outdoor air temperature, °C')
    # chart2.show()


if __name__ == '__main__':
    main()
