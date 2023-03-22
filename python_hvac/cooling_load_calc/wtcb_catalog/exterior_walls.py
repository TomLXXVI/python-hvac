"""
02.A. CONSTRUCTION ASSEMBLIES: EXTERIOR WALLS

Create construction assemblies and store them on the construction assemblies
shelf.
"""
import pandas as pd
from hvac import Quantity
from hvac.cooling_load_calc import (
    Geometry,
    Material,
    HeatFlowDirection,
    SurfaceLayer,
    BuildingComponent,
    AirSpace,
    MechanicalFastening,
    ConstructionAssembly
)

from hvac.cooling_load_calc.wtcb_catalog import (
    MaterialsShelf,
    ConstructionAssembliesShelf,
    db_path
)

Q_ = Quantity


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF1(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(15.0, 'cm')),
        material=MaterialsShelf.load('blokken cellenbeton, gelijmd, 600 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF1 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F2
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF2(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('blokken gebakken aarde, 1200 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF2 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F3
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF3(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(19.0, 'cm')),
        material=MaterialsShelf.load('blokken gebakken aarde, 1200 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF3 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F4
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF4(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('volle betonblokken, geëxpandeerde klei')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF4 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F5
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF5(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('volle blokken halfzwaar beton, 1700 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF5 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F6
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF6(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('holle betonblokken, geëxpandeerde klei, 1200 kg/m3, t=14 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF6 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F7
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF7(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('polyurethaan (PUR), plaat, 150 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('holle betonblokken, geëxpandeerde klei, 1200 kg/m3, t=14 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF7 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F8
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF8(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(19.0, 'cm')),
        material=MaterialsShelf.load('holle betonblokken, geëxpandeerde klei, 1200 kg/m3, t=19 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF8 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F9
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF9(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('holle blokken zwaar beton, 1400 kg/m3, t=14 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF9 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F10
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF10(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('polyurethaan (PUR), plaat, 150 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('holle blokken zwaar beton, 1400 kg/m3, t=14 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF10 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F11
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF11(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialsShelf.load('volle blokken halfzwaar beton, 1800 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(15.0, 'cm')),
        material=MaterialsShelf.load('blokken cellenbeton, gelijmd, 600 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF11 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F12
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF12(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialsShelf.load('volle blokken halfzwaar beton, 1800 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('holle betonblokken, geëxpandeerde klei, 1200 kg/m3, t=14 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF12 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F13
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF13(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialsShelf.load('volle blokken halfzwaar beton, 1800 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('holle blokken zwaar beton, 1400 kg/m3, t=14 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF13 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F14
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF14(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    wood_siding = BuildingComponent.create(
        ID='wood_siding',
        geometry=Geometry(t=Q_(2.5, 'cm')),
        material=Material(
            k=Q_(0.18, 'W / (m * K)'),
            rho=Q_(540, 'kg / m ** 3'),
            c=Q_(1.17, 'kJ / (kg * K)')
        )
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(15.0, 'cm')),
        material=MaterialsShelf.load('blokken cellenbeton, gelijmd, 600 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF14 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            wood_siding,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1500, 'mm ** 2'),  # well ventilated
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    # ext_wall.layers['wood_siding'].slices = 10 --> removed from construction assembly as airspace is well ventilated
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F15
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF15(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    wood_siding = BuildingComponent.create(
        ID='wood_siding',
        geometry=Geometry(t=Q_(2.5, 'cm')),
        material=Material(
            k=Q_(0.18, 'W / (m * K)'),
            rho=Q_(540, 'kg / m ** 3'),
            c=Q_(1.17, 'kJ / (kg * K)')
        )
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=Q_(6, 'cm')),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('holle betonblokken, geëxpandeerde klei, 1200 kg/m3, t=14 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF15 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            wood_siding,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1500, 'mm ** 2'),  # well ventilated
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    # ext_wall.layers['wood_siding'].slices = 10 --> removed from construction assembly as airspace is well ventilated
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F16
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF16(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3')
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    OSB_board = BuildingComponent.create(
        ID='OSB_board',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=MaterialsShelf.load('OSB-plaat')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF16 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            OSB_board,
            insulation,
            gypsum_board,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(500, 'mm ** 2'),  # slightly ventilated
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F17
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF17(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    wood_siding = BuildingComponent.create(
        ID='wood_siding',
        geometry=Geometry(t=Q_(2.5, 'cm')),
        material=Material(
            k=Q_(0.18, 'W / (m * K)'),
            rho=Q_(540, 'kg / m ** 3'),
            c=Q_(1.17, 'kJ / (kg * K)')
        )
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    OSB_board = BuildingComponent.create(
        ID='OSB_board',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=MaterialsShelf.load('OSB-plaat')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF17(t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            wood_siding,
            air_space,
            OSB_board,
            insulation,
            gypsum_board,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1500, 'mm ** 2'),  # well ventilated
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    # ext_wall.layers['wood_siding'].slices = 10 --> removed from construction assembly as airspace is well ventilated
    ext_wall.layers['insulation'].slices = 5
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F18
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF18(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('blokken gebakken aarde, 1200 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF18 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F19
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF19(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('holle betonblokken, geëxpandeerde klei, 1200 kg/m3, t=14 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF19 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 10
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F20
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF20(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    cement_plastering = BuildingComponent.create(
        ID='cement_plastering',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=MaterialsShelf.load('cementpleister')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('blokken cellenbeton, gelijmd, 600 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF20 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            cement_plastering,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['cement_plastering'].slices = 2
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F21
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF21(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    cement_plastering = BuildingComponent.create(
        ID='cement_plastering',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=MaterialsShelf.load('cementpleister')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialsShelf.load('volle betonblokken, geëxpandeerde klei')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF21 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            cement_plastering,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['cement_plastering'].slices = 2
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['inner_leaf'].slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F22
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF22(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    vinyl_siding = BuildingComponent.create(
        ID='vinyl_siding',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=Material(
            k=Q_(0.16, 'W / (m * K)'),
            rho=Q_(1375.0, 'kg / m ** 3'),
            c=Q_(250.0, 'J / (kg * K)')
        )
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(29.0, 'cm')),
        material=MaterialsShelf.load('blokken gebakken aarde, 1500 kg/m3')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('geëxtrudeerd polystyreen (XPS), plaat')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF22 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            vinyl_siding,
            air_space,
            outer_leaf,
            insulation,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1500, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 15
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['gypsum_layer'].slices = 2
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F23
# ------------------------------------------------------------------------------

def create_ca_ext_wall_wtcbF23(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    vinyl_siding = BuildingComponent.create(
        ID='vinyl_siding',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=Material(
            k=Q_(0.16, 'W / (m * K)'),
            rho=Q_(1375.0, 'kg / m ** 3'),
            c=Q_(250.0, 'J / (kg * K)')
        )
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(29.0, 'cm')),
        material=MaterialsShelf.load('blokken gebakken aarde, 1500 kg/m3')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('polyurethaan (PUR), plaat, 150 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF22 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            vinyl_siding,
            air_space,
            outer_leaf,
            insulation,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1500, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    ext_wall.layers['outer_leaf'].slices = 15
    ext_wall.layers['insulation'].slices = 5
    ext_wall.layers['gypsum_layer'].slices = 2
    return ext_wall


# ------------------------------------------------------------------------------

def main():
    t_ins = Q_(12, 'cm')

    ca_ext_wall_wtcbF1 = create_ca_ext_wall_wtcbF1(t_ins)
    ca_ext_wall_wtcbF2 = create_ca_ext_wall_wtcbF2(t_ins)
    ca_ext_wall_wtcbF3 = create_ca_ext_wall_wtcbF3(t_ins)
    ca_ext_wall_wtcbF4 = create_ca_ext_wall_wtcbF4(t_ins)
    ca_ext_wall_wtcbF5 = create_ca_ext_wall_wtcbF5(t_ins)
    ca_ext_wall_wtcbF6 = create_ca_ext_wall_wtcbF6(t_ins)
    ca_ext_wall_wtcbF7 = create_ca_ext_wall_wtcbF7(t_ins)
    ca_ext_wall_wtcbF8 = create_ca_ext_wall_wtcbF8(t_ins)
    ca_ext_wall_wtcbF9 = create_ca_ext_wall_wtcbF9(t_ins)
    ca_ext_wall_wtcbF10 = create_ca_ext_wall_wtcbF10(t_ins)
    ca_ext_wall_wtcbF11 = create_ca_ext_wall_wtcbF11(t_ins)
    ca_ext_wall_wtcbF12 = create_ca_ext_wall_wtcbF12(t_ins)
    ca_ext_wall_wtcbF13 = create_ca_ext_wall_wtcbF13(t_ins)
    ca_ext_wall_wtcbF14 = create_ca_ext_wall_wtcbF14(t_ins)
    ca_ext_wall_wtcbF15 = create_ca_ext_wall_wtcbF15(t_ins)
    ca_ext_wall_wtcbF16 = create_ca_ext_wall_wtcbF16(t_ins)
    ca_ext_wall_wtcbF17 = create_ca_ext_wall_wtcbF17(t_ins)
    ca_ext_wall_wtcbF18 = create_ca_ext_wall_wtcbF18(t_ins)
    ca_ext_wall_wtcbF19 = create_ca_ext_wall_wtcbF19(t_ins)
    ca_ext_wall_wtcbF20 = create_ca_ext_wall_wtcbF20(t_ins)
    ca_ext_wall_wtcbF21 = create_ca_ext_wall_wtcbF21(t_ins)
    ca_ext_wall_wtcbF22 = create_ca_ext_wall_wtcbF22(t_ins)
    ca_ext_wall_wtcbF23 = create_ca_ext_wall_wtcbF23(t_ins)

    # print(ca_ext_wall_wtcbF1, end='\n\n')
    # print(ca_ext_wall_wtcbF2, end='\n\n')
    # print(ca_ext_wall_wtcbF3, end='\n\n')
    # print(ca_ext_wall_wtcbF4, end='\n\n')
    # print(ca_ext_wall_wtcbF5, end='\n\n')
    # print(ca_ext_wall_wtcbF6, end='\n\n')
    # print(ca_ext_wall_wtcbF7, end='\n\n')
    # print(ca_ext_wall_wtcbF8, end='\n\n')
    # print(ca_ext_wall_wtcbF9, end='\n\n')
    # print(ca_ext_wall_wtcbF10, end='\n\n')
    # print(ca_ext_wall_wtcbF11, end='\n\n')
    # print(ca_ext_wall_wtcbF12, end='\n\n')
    # print(ca_ext_wall_wtcbF13, end='\n\n')
    # print(ca_ext_wall_wtcbF14, end='\n\n')
    # print(ca_ext_wall_wtcbF15, end='\n\n')
    # print(ca_ext_wall_wtcbF16, end='\n\n')
    # print(ca_ext_wall_wtcbF17, end='\n\n')
    # print(ca_ext_wall_wtcbF18, end='\n\n')
    # print(ca_ext_wall_wtcbF19, end='\n\n')
    # print(ca_ext_wall_wtcbF20, end='\n\n')
    # print(ca_ext_wall_wtcbF21, end='\n\n')
    # print(ca_ext_wall_wtcbF22, end='\n\n')
    # print(ca_ext_wall_wtcbF23, end='\n\n')

    ConstructionAssembliesShelf.add(
        ca_ext_wall_wtcbF1,
        ca_ext_wall_wtcbF2,
        ca_ext_wall_wtcbF3,
        ca_ext_wall_wtcbF4,
        ca_ext_wall_wtcbF5,
        ca_ext_wall_wtcbF6,
        ca_ext_wall_wtcbF7,
        ca_ext_wall_wtcbF8,
        ca_ext_wall_wtcbF9,
        ca_ext_wall_wtcbF10,
        ca_ext_wall_wtcbF11,
        ca_ext_wall_wtcbF12,
        ca_ext_wall_wtcbF13,
        ca_ext_wall_wtcbF14,
        ca_ext_wall_wtcbF15,
        ca_ext_wall_wtcbF16,
        ca_ext_wall_wtcbF17,
        ca_ext_wall_wtcbF18,
        ca_ext_wall_wtcbF19,
        ca_ext_wall_wtcbF20,
        ca_ext_wall_wtcbF21,
        ca_ext_wall_wtcbF22,
        ca_ext_wall_wtcbF23
    )

    with pd.option_context(
            'display.max_rows', None,
            'display.max_columns', None,
            'display.width', None,
            'display.colheader_justify', 'center'
    ):
        print(ConstructionAssembliesShelf.overview(detailed=True))

    ConstructionAssembliesShelf.export_to_excel(str(db_path / 'construction_assemblies.ods'))


if __name__ == '__main__':
    main()
