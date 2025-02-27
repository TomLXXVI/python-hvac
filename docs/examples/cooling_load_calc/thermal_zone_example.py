"""
Modeling of a Building Thermal Zone.

In this module, the creation of a thermal zone model is programmed for 
demonstration purposes only. The process has been split up in a number of steps  
(using functions).

`create_exterior_wall(...)`:
    Creates an `ExteriorBuildingElement` object that models an exterior wall
    that surrounds the thermal zone.

`add_window(...)`:
    Adds a `Window` object to an `ExteriorBuildingElement` object.

`add_ext_door(...)`:
    Adds an exterior door (also an `ExteriorBuildingElement` object) to an
    exterior wall.

`create_interior_wall(...)`:
    Creates an `InteriorBuildingElement` object that models an interior wall
    that separates the thermal zone from an adjacent thermal zone.

`add_int_door(...)`:
    Adds an interior door (also an `InteriorBuildingElement` object) to an
    interior wall.

`create_roof(...)`:
    Creates a roof (`ExteriorBuildingElement` object) to put on the thermal 
    zone. 

`create_zone_envelope(...)`:
    Uses the functions above to develop all the exterior and interior building 
    elements that surround the thermal zone.

`create_internal_heat_gains(...)`:
    Creates the internal heat gains in the thermal zone (`InternalHeatGain`
    objects).

`create_zone(...)`:
    Creates a `ThermalZone` object and uses the functions 
    `create_zone_envelope(...)` and `create_internal_heat_gains(...)` to add
    the exterior and interior building elements around the zone and to add the
    internal heat gains to the zone. Also adds ventilation to the thermal zone.
    The function returns the fully configured `ThermalZone` object, which can
    then be used for running a simulation.


Notes
-----
This script uses building construction catalogs from package `wtcb`. However, 
before these catalogs can be used, the building materials shelf and window 
properties shelf must have been installed first (by default in a directory 
`wtcb-database` in the user's home directory, see `wtcb.setup.py`). To install
these shelves, you can run the scripts `materials.py` and `windows.py` in the 
`wtcb` package, or call the `main()` functions in these modules (as done in
notebook no. 5).
"""
from hvac import Quantity
from hvac.cooling_load_calc.construction_data import (
    ExteriorWallCatalog,
    InteriorWallCatalog,
    RoofCatalog,
    WindowPropertiesShelf
)
from hvac.cooling_load_calc import (
    ExteriorBuildingElement, 
    WeatherData,
    ConstructionAssembly,
    InteriorBuildingElement,
    PeopleHeatGain,
    InternalHeatGain,
    convert_to_clock_time,
    VentilationZone,
    ThermalZone
)


Q_ = Quantity
ext_wall_catalog = ExteriorWallCatalog()
int_wall_catalog = InteriorWallCatalog()
roof_catalog = RoofCatalog()


def create_exterior_wall(
    name: str,
    width: Quantity,
    height: Quantity,
    azimuth_angle: Quantity,
    weather_data: WeatherData,
    T_int: Quantity
) -> ExteriorBuildingElement:
    area = width * height
    constr_assem_F1 = ext_wall_catalog(
        'F1',
        t_ins=Q_(10, 'cm'),
        T_ext=weather_data.T_db_avg,
        T_int=T_int,
        v_wind=Q_(3, 'm / s')
    )
    exterior_wall = ExteriorBuildingElement.create(
        name=name,
        gross_area=area,
        constr_assem=constr_assem_F1,
        weather_data=weather_data,
        azimuth_angle=azimuth_angle,
        slope_angle=Q_(90, 'deg')
    )
    return exterior_wall


def add_window(
    ext_wall: ExteriorBuildingElement,
    name: str,
    width: Quantity,
    height: Quantity,
) -> None:
    wnd_props = WindowPropertiesShelf.load('window-5a-operable-wood/vinyl')
    ext_wall.add_window(
        name=name,
        width=width,
        height=height,
        props=wnd_props,
        F_rad=0.33
    )


def add_ext_door(
    ext_wall: ExteriorBuildingElement,
    name: str,
    width: Quantity,
    height: Quantity
) -> None:
    ext_wall.add_door(
        name=name,
        width=width,
        height=height,
        constr_assem=ConstructionAssembly.create(
            ID='ext-door-assem',
            layers=None,
            U=Q_(3.0, 'W / (m**2 * K)'),
        )
    )


def create_interior_wall(
    name: str,
    adj_zone_name: str,
    width: Quantity,
    height: Quantity,
    T_int: Quantity,
    T_adj: Quantity
) -> InteriorBuildingElement:
    area = width * height
    constr_assem_F3 = int_wall_catalog(
        'F3', 
        t_ins=None, 
        T_int=T_int,
        T_adj=T_adj
    )
    int_wall = InteriorBuildingElement.create(
        name=name,
        adjacent_zone_name=adj_zone_name,
        constr_assem=constr_assem_F3,
        gross_area=area,
        adjacent_zone_temperature=Q_(28, 'degC'),
        V_dot_trf=None,
        F_rad=0.46
    )
    return int_wall


def add_int_door(
    int_wall: InteriorBuildingElement,
    name: str,
    width: Quantity,
    height: Quantity,
) -> None:
    int_wall.add_door(
        name=name,
        width=width,
        height=height,
        constr_assem=ConstructionAssembly.create(
            ID='int-door-assem',
            layers=None,
            U=Q_(4.0, 'W / (m**2 * K)')
        ),
        F_rad=0.46
    )


def create_roof(
    name: str,
    width: Quantity,
    height: Quantity,
    azimuth_angle: Quantity,
    weather_data: WeatherData,
    T_int: Quantity
) -> ExteriorBuildingElement:
    area = width * height
    constr_assem_F1 = roof_catalog(
        ID='F1',
        t_ins=Q_(18, 'cm'),
        T_ext=weather_data.T_db_avg,
        T_int=T_int,
        v_wind=Q_(3, 'm / s')
    )
    roof = ExteriorBuildingElement.create(
        name=name,
        gross_area=area,
        constr_assem=constr_assem_F1,
        weather_data=weather_data,
        azimuth_angle=azimuth_angle,
        slope_angle=Q_(30, 'deg')
    )
    return roof


def create_zone_envelope(
    weather_data: WeatherData,
    T_int_des: Quantity
) -> tuple[tuple[ExteriorBuildingElement, ...], tuple[InteriorBuildingElement, ...]]:
    # Exterior South wall
    ext_wall_south = create_exterior_wall(
        name='ew_S',
        width=Q_(10, 'm'),
        height=Q_(4, 'm'),
        azimuth_angle=Q_(0, 'deg'),
        weather_data=weather_data,
        T_int=T_int_des
    )
    for window_name in ['wnd_S1', 'wnd_S2']:
        add_window(
            ext_wall=ext_wall_south,
            name=window_name,
            width=Q_(3, 'm'),
            height=Q_(2, 'm')
        )
    # Exterior East wall
    ext_wall_east = create_exterior_wall(
        name='ew_E',
        width=Q_(6, 'm'),
        height=Q_(4, 'm'),
        azimuth_angle=Q_(90, 'deg'),
        weather_data=weather_data,
        T_int=T_int_des
    )
    add_window(
        ext_wall=ext_wall_east,
        name='wnd_E1',
        width=Q_(2, 'm'),
        height=Q_(2, 'm')
    )
    add_ext_door(
        ext_wall=ext_wall_east,
        name='ed_E1',
        width=Q_(1, 'm'),
        height=Q_(2.2, 'm')
    )
    # Exterior West wall
    ext_wall_west = create_exterior_wall(
        name='ew_W',
        width=Q_(6, 'm'),
        height=Q_(4, 'm'),
        azimuth_angle=Q_(-90, 'deg'),
        weather_data=weather_data,
        T_int=T_int_des
    )
    # Interior North wall
    int_wall_north = create_interior_wall(
        name='iw_N',
        adj_zone_name='zone_2',
        width=Q_(10, 'm'),
        height=Q_(4, 'm'),
        T_int=T_int_des,
        T_adj=T_int_des.to('K') + Q_(5, 'K') 
    )
    add_int_door(
        int_wall=int_wall_north,
        name='id_N1',
        width=Q_(0.8, 'm'),
        height=Q_(2.1, 'm')
    )
    # Roof with slope oriented to the South
    roof = create_roof(
        name='rf_S',
        width=Q_(10, 'm'),
        height=Q_(6, 'm'),
        azimuth_angle=Q_(0, 'deg'),
        weather_data=weather_data,
        T_int=T_int_des
    )
    add_window(
        ext_wall=roof,
        name='wnd_S3',
        width=Q_(1, 'm'),
        height=Q_(1, 'm')
    )
    ext_build_elems = (ext_wall_west, ext_wall_south, ext_wall_east, roof)
    int_build_elems = (int_wall_north,)
    return ext_build_elems, int_build_elems


def create_internal_heat_gains(weather_data: WeatherData) -> tuple[InternalHeatGain, ...]:
    
    def occupancy_schedule(t_sol_sec: float) -> int:
        dt_clock, _ = convert_to_clock_time(
            t_sol_sec,
            date=weather_data.location.date,
            L_loc=weather_data.location.L_loc,
            tz_loc=weather_data.location.timezone
        )
        t_clock = dt_clock.time()
        if 8 <= t_clock.hour < 12:
            return 3
        elif 12 <= t_clock.hour < 18:
            return 5
        else:
            return 0 
    
    ihg_people = PeopleHeatGain.create(
        name='people',
        Q_dot_sen_person=Q_(75, 'W'),
        Q_dot_lat_person=Q_(55, 'W'),
        F_rad=0.58,
        schedule=occupancy_schedule
    )
    return (ihg_people,)


def create_zone(weather_data: WeatherData, T_int_des: Quantity) -> ThermalZone:
    zone = ThermalZone.create(
        name='zone_1',
        floor_area=Q_(60, 'm ** 2'),
        ceiling_height=Q_(3, 'm'),
        unit_R_im=Q_(0.37, 'm**2 * K / W'),
        unit_C_im=Q_(100.0, 'kJ / (m**2 * K)'),
        weather_data=weather_data,
        factor_surf_im=2.0
    )
    ext_build_elems, int_build_elems = create_zone_envelope(weather_data, T_int_des)
    int_heat_gains = create_internal_heat_gains(weather_data)
    ventilation_zone = VentilationZone.create(name='vent_zone1')  # use defaults
    zone.add_ext_build_elems(*ext_build_elems)
    zone.add_int_build_elems(*int_build_elems)
    zone.add_internal_heat_gains(*int_heat_gains)
    zone.add_ventilation(ventilation_zone)  # only minimal ventilation
    return zone
