"""
02.B. CONSTRUCTION ASSEMBLIES: ROOFS AND CEILINGS

Create construction assemblies and store them on the construction assemblies
shelf.
"""
import pandas as pd
from hvac import Quantity
from hvac.cooling_load_calc import (
    Geometry,
    HeatFlowDirection,
    SurfaceLayer,
    BuildingComponent,
    AirSpace,
    ConstructionAssembly
)

from hvac.cooling_load_calc.wtcb_catalog import (
    MaterialsShelf,
    ConstructionAssembliesShelf,
    db_path
)

Q_ = Quantity


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_ca_roof_wtcbF1(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection = HeatFlowDirection.UPWARDS,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    roof_tiles = BuildingComponent.create(
        ID='roof_tiles',
        geometry=Geometry(t=Q_(3.0, 'cm')),
        material=MaterialsShelf.load('kleidakpannen')
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(45.0, 'deg')
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.88, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.12, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    plywood_board = BuildingComponent.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('multiplexplaat')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcbF1 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            roof_tiles,
            air_space,
            insulation,
            plywood_board,
            int_surf_film
        ]
    )
    roof = roof.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1499, 'mm ** 2'),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['roof_tiles'].slices = 2
    roof.layers['insulation'].slices = 5
    roof.layers['plywood_board'].slices = 2
    return roof


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F2
# ------------------------------------------------------------------------------

def create_ca_roof_wtcbF2(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    roof_tiles = BuildingComponent.create(
        ID='roof_tiles',
        geometry=Geometry(t=Q_(3.0, 'cm')),
        material=MaterialsShelf.load('kleidakpannen')
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(45.0, 'deg')
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.88, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.12, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcbF2 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            roof_tiles,
            air_space,
            insulation,
            gypsum_board,
            int_surf_film
        ]
    )
    roof = roof.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1499, 'mm ** 2'),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['roof_tiles'].slices = 2
    roof.layers['insulation'].slices = 5
    roof.layers['gypsum_board'].slices = 2
    return roof


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F3
# ------------------------------------------------------------------------------

def create_ca_roof_wtcbF3(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    roof_tiles = BuildingComponent.create(
        ID='roof_tiles',
        geometry=Geometry(t=Q_(3.0, 'cm')),
        material=MaterialsShelf.load('kleidakpannen')
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(45.0, 'deg')
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.88, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.12, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcbF3 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            roof_tiles,
            air_space,
            insulation,
            gypsum_board,
            int_surf_film
        ]
    )
    roof = roof.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(501, 'mm ** 2'),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['roof_tiles'].slices = 2
    roof.layers['insulation'].slices = 5
    roof.layers['gypsum_board'].slices = 2
    return roof


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F4
# ------------------------------------------------------------------------------

def create_ca_roof_wtcbF4(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    roofing = BuildingComponent.create(
        ID='roofing',
        geometry=Geometry(t=Q_(2.5, 'mm')),
        material=MaterialsShelf.load('EPDM')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('geëxpandeerd polystyreen (EPS), plaat')
    )
    plywood_board = BuildingComponent.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(2, 'cm')),
        material=MaterialsShelf.load('multiplexplaat')
    )
    airspace_ = AirSpace.create(
        ID='airspace_',
        geometry=Geometry(t=Q_(15, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        dT=dT_asp,
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    # airspace with width less than 10 times the thickness = air void with
    # thermal resistance calculated acc. to NBN EN ISO 6946, annex D.4.
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=Q_(15, 'cm'), A=Q_(0.11, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    # per unit of area (1 m²) 11 % is wooden framework and 89 % is airspace.
    airspace = airspace_ // wood_
    airspace.ID = 'airspace'
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcbF4 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            roofing,
            insulation,
            plywood_board,
            airspace,
            gypsum_board,
            int_surf_film
        ]
    )
    roof = roof.apply_ventilated_layer_correction(
        ventilated_layer_ID=airspace.ID,
        area_of_openings=Q_(501, 'mm ** 2'),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['insulation'].slices = 5
    roof.layers['plywood_board'].slices = 2
    roof.layers['gypsum_board'].slices = 2
    return roof


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F5
# ------------------------------------------------------------------------------

def create_ca_roof_wtcbF5(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    roofing = BuildingComponent.create(
        ID='roofing',
        geometry=Geometry(t=Q_(2.5, 'mm')),
        material=MaterialsShelf.load('EPDM')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('geëxpandeerd polystyreen (EPS), plaat')
    )
    slope_concrete = BuildingComponent.create(
        ID='slope_concrete',
        geometry=Geometry(t=Q_(7, 'cm')),
        material=MaterialsShelf.load('licht beton, 1900 kg/m3')
    )
    floor_slabs = BuildingComponent.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('welfsel, zwaar beton, t=12 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcbF5 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            roofing,
            insulation,
            slope_concrete,
            floor_slabs,
            gypsum_layer,
            int_surf_film
        ]
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['insulation'].slices = 5
    roof.layers['slope_concrete'].slices = 5
    roof.layers['floor_slabs'].slices = 10
    roof.layers['gypsum_layer'].slices = 2
    return roof


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F6
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF6(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    multiplex_board = BuildingComponent.create(
        ID='multiplex_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('multiplexplaat')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF6 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            insulation,
            multiplex_board,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F7
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF7(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF7 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            insulation,
            gypsum_board,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F8
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF8(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    airspace = AirSpace.create(
        ID='airspace',
        geometry=Geometry(t=Q_(5, 'cm')),
        dT=dT_asp,
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF8 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            insulation,
            airspace,
            gypsum_board,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F9
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF9(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    plywood_board = BuildingComponent.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(2, 'cm')),
        material=MaterialsShelf.load('multiplexplaat')
    )
    airspace_ = AirSpace.create(
        ID='airspace_',
        geometry=Geometry(t=Q_(5, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        dT=dT_asp,
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    airspace = airspace_ // wood_
    airspace.ID = 'airspace'
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF9 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            plywood_board,
            airspace,
            insulation,
            gypsum_board,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F10
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF10(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    plywood_board = BuildingComponent.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(2, 'cm')),
        material=MaterialsShelf.load('multiplexplaat')
    )
    airspace_ = AirSpace.create(
        ID='airspace_',
        geometry=Geometry(t=Q_(5, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        dT=dT_asp,
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    airspace_01 = airspace_ // wood_
    airspace_01.ID = 'airspace_01'
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    airspace_02 = AirSpace.create(
        ID='airspace_02',
        geometry=Geometry(t=Q_(5, 'cm')),
        dT=dT_asp,
        heat_flow_direction=heat_flow_direction,
        Tmn=T_asp
    )
    gypsum_board = BuildingComponent.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF10 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            plywood_board,
            airspace_01,
            insulation,
            airspace_02,
            gypsum_board,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F11
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF11(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    concrete_slab = BuildingComponent.create(
        ID='concrete_slab',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('beton, gewapend, 2400 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF11 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            insulation,
            concrete_slab,
            gypsum_layer,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F12
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF12(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, onbeklede plaat, 135 kg/m3')
    )
    floor_slabs = BuildingComponent.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('welfsel, zwaar beton, t=12 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF12 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            insulation,
            floor_slabs,
            gypsum_layer,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F13
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF13(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    screed_concrete = BuildingComponent.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(6, 'cm')),
        material=MaterialsShelf.load('licht beton, 1600 kg/m3')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('geëxtrudeerd polystyreen (XPS), plaat')
    )
    floor_slabs = BuildingComponent.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('welfsel, zwaar beton, t=12 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF13 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            screed_concrete,
            insulation,
            floor_slabs,
            gypsum_layer,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F14
# ------------------------------------------------------------------------------

def create_ca_ceiling_wtcbF14(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext
    )
    screed_concrete = BuildingComponent.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(6, 'cm')),
        material=MaterialsShelf.load('licht beton, 1600 kg/m3')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('polyurethaan (PUR), plaat, 150 kg/m3')
    )
    floor_slabs = BuildingComponent.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('welfsel, zwaar beton, t=12 cm')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcbF14 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            screed_concrete,
            insulation,
            floor_slabs,
            gypsum_layer,
            int_surf_film
        ]
    )
    ceiling = ceiling.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return ceiling


# ------------------------------------------------------------------------------

def main():
    t_ins = Q_(12, 'cm')
    hfd = HeatFlowDirection.DOWNWARDS

    ca_roof_wtcbF1 = create_ca_roof_wtcbF1(t_ins, hfd)
    ca_roof_wtcbF2 = create_ca_roof_wtcbF2(t_ins, hfd)
    ca_roof_wtcbF3 = create_ca_roof_wtcbF3(t_ins, hfd)
    ca_roof_wtcbF4 = create_ca_roof_wtcbF4(t_ins, hfd)
    ca_roof_wtcbF5 = create_ca_roof_wtcbF5(t_ins, hfd)
    ca_ceiling_wtcbF6 = create_ca_ceiling_wtcbF6(t_ins, hfd)
    ca_ceiling_wtcbF7 = create_ca_ceiling_wtcbF7(t_ins, hfd)
    ca_ceiling_wtcbF8 = create_ca_ceiling_wtcbF8(t_ins, hfd)
    ca_ceiling_wtcbF9 = create_ca_ceiling_wtcbF9(t_ins, hfd)
    ca_ceiling_wtcbF10 = create_ca_ceiling_wtcbF10(t_ins, hfd)
    ca_ceiling_wtcbF11 = create_ca_ceiling_wtcbF11(t_ins, hfd)
    ca_ceiling_wtcbF12 = create_ca_ceiling_wtcbF12(t_ins, hfd)
    ca_ceiling_wtcbF13 = create_ca_ceiling_wtcbF13(t_ins, hfd)
    ca_ceiling_wtcbF14 = create_ca_ceiling_wtcbF14(t_ins, hfd)

    ConstructionAssembliesShelf.add(
        ca_roof_wtcbF1,
        ca_roof_wtcbF2,
        ca_roof_wtcbF3,
        ca_roof_wtcbF4,
        ca_roof_wtcbF5,
        ca_ceiling_wtcbF6,
        ca_ceiling_wtcbF7,
        ca_ceiling_wtcbF8,
        ca_ceiling_wtcbF9,
        ca_ceiling_wtcbF10,
        ca_ceiling_wtcbF11,
        ca_ceiling_wtcbF12,
        ca_ceiling_wtcbF13,
        ca_ceiling_wtcbF14
    )

    # print(ca_roof_wtcbF1, end='\n\n')
    # print(ca_roof_wtcbF2, end='\n\n')
    # print(ca_roof_wtcbF3, end='\n\n')
    # print(ca_roof_wtcbF4, end='\n\n')
    # print(ca_roof_wtcbF5, end='\n\n')
    # print(ca_ceiling_wtcbF6, end='\n\n')
    # print(ca_ceiling_wtcbF7, end='\n\n')
    # print(ca_ceiling_wtcbF8, end='\n\n')
    # print(ca_ceiling_wtcbF9, end='\n\n')
    # print(ca_ceiling_wtcbF10, end='\n\n')
    # print(ca_ceiling_wtcbF11, end='\n\n')
    # print(ca_ceiling_wtcbF12, end='\n\n')
    # print(ca_ceiling_wtcbF13, end='\n\n')
    # print(ca_ceiling_wtcbF14, end='\n\n')

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
