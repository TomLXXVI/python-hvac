"""
02.B. CONSTRUCTION ASSEMBLIES: ROOFS AND CEILINGS
Creates the construction assemblies and saves them on the construction assembly
shelf.
"""
import pandas as pd
from hvac import Quantity
from hvac.cooling_load_calc.core import (
    Geometry,
    HeatFlowDirection,
    SurfaceFilm,
    SolidLayer,
    AirLayer,
    ConstructionAssembly
)

from hvac.cooling_load_calc.wtcb.setup import (
    MaterialShelf,
    ConstructionAssemblyShelf,
    db_path
)

Q_ = Quantity


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_roof_wtcb_F1(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection = HeatFlowDirection.UPWARDS,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    roof_tiles = SolidLayer.create(
        ID='roof_tiles',
        geometry=Geometry(t=Q_(3.0, 'cm')),
        material=MaterialShelf.load('roof-tile-clay')
    )
    air_space = AirLayer.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp,
        surf_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        angle=Q_(45.0, 'deg')
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.88, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.12, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    plywood_board = SolidLayer.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('plywood')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F1_t_ins={t_ins.to("cm"):~P.0f}',
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
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['roof_tiles'].num_slices = 2
    roof.layers['insulation'].num_slices = 5
    roof.layers['plywood_board'].num_slices = 2
    return roof


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F2
# ------------------------------------------------------------------------------

def create_roof_wtcb_F2(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    roof_tiles = SolidLayer.create(
        ID='roof_tiles',
        geometry=Geometry(t=Q_(3.0, 'cm')),
        material=MaterialShelf.load('roof-tile-clay')
    )
    air_space = AirLayer.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp,
        surf_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        angle=Q_(45.0, 'deg')
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.88, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.12, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F2_t_ins={t_ins.to("cm"):~P.0f}',
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
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['roof_tiles'].num_slices = 2
    roof.layers['insulation'].num_slices = 5
    roof.layers['gypsum_board'].num_slices = 2
    return roof


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F3
# ------------------------------------------------------------------------------

def create_roof_wtcb_F3(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    roof_tiles = SolidLayer.create(
        ID='roof_tiles',
        geometry=Geometry(t=Q_(3.0, 'cm')),
        material=MaterialShelf.load('roof-tile-clay')
    )
    air_space = AirLayer.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp,
        surf_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        angle=Q_(45.0, 'deg')
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.88, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.12, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F3_t_ins={t_ins.to("cm"):~P.0f}',
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
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['roof_tiles'].num_slices = 2
    roof.layers['insulation'].num_slices = 5
    roof.layers['gypsum_board'].num_slices = 2
    return roof


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F4
# ------------------------------------------------------------------------------

def create_roof_wtcb_F4(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    roofing = SolidLayer.create(
        ID='roofing',
        geometry=Geometry(t=Q_(2.5, 'mm')),
        material=MaterialShelf.load('epdm')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polystyrene-expanded-sheet')
    )
    plywood_board = SolidLayer.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(2, 'cm')),
        material=MaterialShelf.load('plywood')
    )
    air_layer_ = AirLayer.create(
        ID='air_layer_',
        geometry=Geometry(t=Q_(15, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    # AirLayer with width less than 10 times the thickness = air void with
    # thermal resistance calculated acc. to NBN EN ISO 6946, annex D.4.
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=Q_(15, 'cm'), A=Q_(0.11, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    # per unit of area (1 m²) 11 % is wooden framework and 89 % is AirLayer.
    air_layer = air_layer_ // wood_
    air_layer.ID = 'air_layer'
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F4_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            ext_surf_film,
            roofing,
            insulation,
            plywood_board,
            air_layer,
            gypsum_board,
            int_surf_film
        ]
    )
    roof = roof.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_layer.ID,
        area_of_openings=Q_(501, 'mm ** 2'),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    roof = roof.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    roof.layers['insulation'].num_slices = 5
    roof.layers['plywood_board'].num_slices = 2
    roof.layers['gypsum_board'].num_slices = 2
    return roof


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F5
# ------------------------------------------------------------------------------

def create_roof_wtcb_F5(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    roofing = SolidLayer.create(
        ID='roofing',
        geometry=Geometry(t=Q_(2.5, 'mm')),
        material=MaterialShelf.load('epdm')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polystyrene-expanded-sheet')
    )
    slope_concrete = SolidLayer.create(
        ID='slope_concrete',
        geometry=Geometry(t=Q_(7, 'cm')),
        material=MaterialShelf.load('concrete-light-1900kg/m3')
    )
    floor_slabs = SolidLayer.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('precast-slab-heavy-concrete-t=12cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F5_t_ins={t_ins.to("cm"):~P.0f}',
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
    roof.layers['insulation'].num_slices = 5
    roof.layers['slope_concrete'].num_slices = 5
    roof.layers['floor_slabs'].num_slices = 10
    roof.layers['gypsum_layer'].num_slices = 2
    return roof


# ------------------------------------------------------------------------------
# CEILING CONSTRUCTION ASSEMBLY WTCB F6
# ------------------------------------------------------------------------------

def create_ceiling_wtcb_F6(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    multiplex_board = SolidLayer.create(
        ID='multiplex_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('plywood')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F6_t_ins={t_ins.to("cm"):~P.0f}',
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

def create_ceiling_wtcb_F7(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F7_t_ins={t_ins.to("cm"):~P.0f}',
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

def create_ceiling_wtcb_F8(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    air_layer = AirLayer.create(
        ID='air_layer',
        geometry=Geometry(t=Q_(5, 'cm')),
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F8_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            insulation,
            air_layer,
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

def create_ceiling_wtcb_F9(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    plywood_board = SolidLayer.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(2, 'cm')),
        material=MaterialShelf.load('plywood')
    )
    air_layer_ = AirLayer.create(
        ID='air_layer_',
        geometry=Geometry(t=Q_(5, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    air_layer = air_layer_ // wood_
    air_layer.ID = 'air_layer'
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F9_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            plywood_board,
            air_layer,
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

def create_ceiling_wtcb_F10(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    plywood_board = SolidLayer.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(2, 'cm')),
        material=MaterialShelf.load('plywood')
    )
    AirLayer_ = AirLayer.create(
        ID='AirLayer_',
        geometry=Geometry(t=Q_(5, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.89, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.11, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    AirLayer_01 = AirLayer_ // wood_
    AirLayer_01.ID = 'AirLayer_01'
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    AirLayer_02 = AirLayer.create(
        ID='AirLayer_02',
        geometry=Geometry(t=Q_(5, 'cm')),
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp
    )
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F10_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            plywood_board,
            AirLayer_01,
            insulation,
            AirLayer_02,
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

def create_ceiling_wtcb_F11(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    concrete_slab = SolidLayer.create(
        ID='concrete_slab',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('concrete-reinforced-2400kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F11_t_ins={t_ins.to("cm"):~P.0f}',
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

def create_ceiling_wtcb_F12(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-sheet-uncovered-135kg/m3')
    )
    floor_slabs = SolidLayer.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('precast-slab-heavy-concrete-t=12cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F12_t_ins={t_ins.to("cm"):~P.0f}',
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

def create_ceiling_wtcb_F13(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    screed_concrete = SolidLayer.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(6, 'cm')),
        material=MaterialShelf.load('concrete-light-1600kg/m3')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polystyrene-extruded-sheet')
    )
    floor_slabs = SolidLayer.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('precast-slab-heavy-concrete-t=12cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F13_t_ins={t_ins.to("cm"):~P.0f}',
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

def create_ceiling_wtcb_F14(
    t_ins: Quantity,
    heat_flow_dir: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_ext
    )
    screed_concrete = SolidLayer.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(6, 'cm')),
        material=MaterialShelf.load('concrete-light-1600kg/m3')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polyurethane-sheet-150kg/m3')
    )
    floor_slabs = SolidLayer.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('precast-slab-heavy-concrete-t=12cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F14_t_ins={t_ins.to("cm"):~P.0f}',
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

    ca_roof_wtcb_F1 = create_roof_wtcb_F1(t_ins, hfd)
    ca_roof_wtcb_F2 = create_roof_wtcb_F2(t_ins, hfd)
    ca_roof_wtcb_F3 = create_roof_wtcb_F3(t_ins, hfd)
    ca_roof_wtcb_F4 = create_roof_wtcb_F4(t_ins, hfd)
    ca_roof_wtcb_F5 = create_roof_wtcb_F5(t_ins, hfd)
    ca_ceiling_wtcb_F6 = create_ceiling_wtcb_F6(t_ins, hfd)
    ca_ceiling_wtcb_F7 = create_ceiling_wtcb_F7(t_ins, hfd)
    ca_ceiling_wtcb_F8 = create_ceiling_wtcb_F8(t_ins, hfd)
    ca_ceiling_wtcb_F9 = create_ceiling_wtcb_F9(t_ins, hfd)
    ca_ceiling_wtcb_F10 = create_ceiling_wtcb_F10(t_ins, hfd)
    ca_ceiling_wtcb_F11 = create_ceiling_wtcb_F11(t_ins, hfd)
    ca_ceiling_wtcb_F12 = create_ceiling_wtcb_F12(t_ins, hfd)
    ca_ceiling_wtcb_F13 = create_ceiling_wtcb_F13(t_ins, hfd)
    ca_ceiling_wtcb_F14 = create_ceiling_wtcb_F14(t_ins, hfd)

    ConstructionAssemblyShelf.add(
        ca_roof_wtcb_F1,
        ca_roof_wtcb_F2,
        ca_roof_wtcb_F3,
        ca_roof_wtcb_F4,
        ca_roof_wtcb_F5,
        ca_ceiling_wtcb_F6,
        ca_ceiling_wtcb_F7,
        ca_ceiling_wtcb_F8,
        ca_ceiling_wtcb_F9,
        ca_ceiling_wtcb_F10,
        ca_ceiling_wtcb_F11,
        ca_ceiling_wtcb_F12,
        ca_ceiling_wtcb_F13,
        ca_ceiling_wtcb_F14
    )

    with pd.option_context(
            'display.max_rows', None,
            'display.max_columns', None,
            'display.width', None,
            'display.colheader_justify', 'center'
    ):
        print(ConstructionAssemblyShelf.overview(detailed=True))

    ConstructionAssemblyShelf.export_to_excel(str(db_path / 'construction_assemblies.ods'))


if __name__ == '__main__':
    main()
