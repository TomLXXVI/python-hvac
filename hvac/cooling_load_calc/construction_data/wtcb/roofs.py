"""02.B. CONSTRUCTION ASSEMBLIES: ROOFS AND CEILINGS

Collection of functions that create the construction assemblies of roofs and
ceilings contained in the WTCB catalog (see /docs/wtcb_catalog/wtcb_catalog.pdf).
"""
from hvac import Quantity
from ....cooling_load_calc import (
    Geometry,
    HeatFlowDirection,
    SurfaceFilm,
    SolidLayer,
    AirLayer,
    ConstructionAssembly,
    AirLayerTemperatureSolver
)
from .setup import MaterialShelf


Q_ = Quantity


def _create_air_layer(
    T_ext: Quantity,
    T_int: Quantity,
    R_ea: Quantity,
    R_ai: Quantity,
    geometry: Geometry,
    surf_emissivities: tuple[Quantity, Quantity] = (Q_(0.9, 'frac'), Q_(0.9, 'frac')),
    angle: Quantity = Q_(0, 'deg')
) -> tuple[AirLayer, Quantity]:
    """Internal helper function to create an `AirLayer` object.

    Parameters
    ----------
    T_ext:
        Exterior temperature.
    T_int:
        Interior temperature.
    R_ea:
        Unit thermal resistance between exterior and outer air layer side.
    R_ai:
        Unit thermal resistance between inner air layer side and interior.
    geometry:
        The geometry (size) of the air gap.
    surf_emissivities:
        A tuple with the emissivities of the exterior side surface and the
        interior side surface.
    angle:
        Slope angle of the roof.

    Returns
    -------
    `AirLayer` object and average air gap temperature.
    """
    ats = AirLayerTemperatureSolver(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=R_ea,
        R_ai=R_ai
    )
    if T_ext < T_int:
        dT = T_int.to('K') - T_ext.to('K')
        T_ae_ini = T_ext.to('K') + dT / 2 - Q_(2.5, 'K')
        T_ai_ini = T_int.to('K') - dT / 2 + Q_(2.5, 'K')
        heat_flow_dir = HeatFlowDirection.UPWARDS
    else:
        dT = T_ext.to('K') - T_int.to('K')
        T_ae_ini = T_ext.to('K') - dT / 2 + Q_(2.5, 'K')
        T_ai_ini = T_int.to('K') + dT / 2 - Q_(2.5, 'K')
        heat_flow_dir = HeatFlowDirection.DOWNWARDS
    *_, dT_asp, T_asp = ats.solve(T_ae_ini, T_ai_ini)

    air_space = AirLayer.create(
        ID='air_space',
        geometry=geometry,
        dT=dT_asp,
        heat_flow_dir=heat_flow_dir,
        T_mn=T_asp,
        surf_emissivities=surf_emissivities,
        angle=angle
    )
    return air_space, T_asp


# ------------------------------------------------------------------------------
# ROOF CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------
def create_roof_F1(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_ext < T_int
        else HeatFlowDirection.DOWNWARDS
    )
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
    # air_space: see below
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
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + roof_tiles.R,
        R_ai=insulation.R + plywood_board.R + int_surf_film.R,
        geometry=Geometry(t=Q_(5, 'cm')),
        angle=Q_(45, 'deg')
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F1',
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

def create_roof_F2(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_ext < T_int
        else HeatFlowDirection.DOWNWARDS
    )
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
    # air_space: see below
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
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + roof_tiles.R,
        R_ai=insulation.R + gypsum_board.R + int_surf_film.R,
        geometry=Geometry(t=Q_(5, 'cm')),
        angle=Q_(45, 'deg')
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F2',
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

def create_roof_F3(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_ext < T_int
        else HeatFlowDirection.DOWNWARDS
    )
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
    # air_space: see below
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
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + roof_tiles.R,
        R_ai=insulation.R + gypsum_board.R + int_surf_film.R,
        geometry=Geometry(t=Q_(5, 'cm')),
        angle=Q_(45, 'deg')
    )
    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F3',
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

def create_roof_F4(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_ext < T_int
        else HeatFlowDirection.DOWNWARDS
    )
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
    # Air layer having width less than 10 times the thickness = air void with
    # thermal resistance calculated acc. to NBN EN ISO 6946, annex D.4.
    # Per unit area 11 % is wooden framework and 89 % is air.
    air_layer_, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + roofing.R + insulation.R + plywood_board.R,
        R_ai=gypsum_board.R + int_surf_film.R,
        geometry=Geometry(t=Q_(15, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        angle=Q_(0, 'deg')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=Q_(15, 'cm'), A=Q_(0.11, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    air_layer = air_layer_ // wood_
    air_layer.ID = 'air_layer'

    roof = ConstructionAssembly.create(
        ID=f'roof_wtcb_F4',
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

def create_roof_F5(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_ext < T_int
        else HeatFlowDirection.DOWNWARDS
    )
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
        ID=f'roof_wtcb_F5',
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

def create_ceiling_F6(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
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
        ID=f'ceiling_wtcb_F6',
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

def create_ceiling_F7(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
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
        ID=f'ceiling_wtcb_F7',
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

def create_ceiling_F8(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
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
    # air_layer: see below
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
    air_layer, T_asp = _create_air_layer(
        T_ext=T_adj,
        T_int=T_int,
        R_ea=adj_surf_film.R + insulation.R,
        R_ai=gypsum_board.R + int_surf_film.R,
        geometry=Geometry(t=Q_(5, 'cm')),
        angle=Q_(0, 'deg')
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F8',
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

def create_ceiling_F9(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
    )
    plywood_board = SolidLayer.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(2, 'cm')),
        material=MaterialShelf.load('plywood')
    )
    # air layer: see below
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
    air_layer_, T_asp = _create_air_layer(
        T_ext=T_adj,
        T_int=T_int,
        R_ea=adj_surf_film.R + plywood_board.R,
        R_ai=insulation.R + gypsum_board.R + int_surf_film.R,
        geometry=Geometry(t=Q_(5, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        angle=Q_(0, 'deg')
    )
    air_layer = air_layer_ // wood_
    air_layer.ID = 'air_layer'
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F9',
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

def create_ceiling_F10(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
    )
    plywood_board = SolidLayer.create(
        ID='plywood_board',
        geometry=Geometry(t=Q_(2, 'cm')),
        material=MaterialShelf.load('plywood')
    )
    # air layer 1: see below
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
    # air layer 2: see below
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
    air_layer_01_, T_asp = _create_air_layer(
        T_ext=T_adj,
        T_int=T_int,
        R_ea=adj_surf_film.R + plywood_board.R,
        R_ai=insulation.R + gypsum_board.R + int_surf_film.R,
        geometry=Geometry(t=Q_(5, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
        angle=Q_(0, 'deg')
    )
    air_layer_01 = air_layer_01_ // wood_
    air_layer_01.ID = 'air_layer_01'
    air_layer_02, _ = _create_air_layer(
        T_ext=T_adj,
        T_int=T_int,
        R_ea=adj_surf_film.R + plywood_board.R + insulation.R,
        R_ai=gypsum_board.R + int_surf_film.R,
        geometry=Geometry(t=Q_(5, 'cm'))
    )
    ceiling = ConstructionAssembly.create(
        ID=f'ceiling_wtcb_F10',
        layers=[
            adj_surf_film,
            plywood_board,
            air_layer_01,
            insulation,
            air_layer_02,
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

def create_ceiling_F11(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
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
        ID=f'ceiling_wtcb_F11',
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

def create_ceiling_F12(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
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
        ID=f'ceiling_wtcb_F12',
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

def create_ceiling_F13(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
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
        ID=f'ceiling_wtcb_F13',
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

def create_ceiling_F14(
    t_ins: Quantity,
    T_adj: Quantity = Q_(15.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.UPWARDS if T_adj < T_int
        else HeatFlowDirection.DOWNWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
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
    floor_slab = SolidLayer.create(
        ID='floor_slab',
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
        ID=f'ceiling_wtcb_F14',
        layers=[
            adj_surf_film,
            screed_concrete,
            insulation,
            floor_slab,
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
class RoofCatalog:
    """Class that bundles the functions to create construction assemblies of
    roofs.
    """
    def __init__(
        self,
        t_ins: Quantity = Q_(10, 'cm'),
        T_ext: Quantity = Q_(0, 'degC'),
        T_int: Quantity = Q_(20, 'degC'),
        v_wind: Quantity = Q_(4, 'm / s')
    ) -> None:
        """Creates an instance of `RoofCatalog`.

        Parameters
        ----------
        t_ins:
            Thickness of the insulation layer.
        T_ext:
            The outdoor air design temperature.
        T_int:
            The indoor air design temperature.
        v_wind:
            Design value for the wind speed.
        """
        self._d = {
            'F1': create_roof_F1,
            'F2': create_roof_F2,
            'F3': create_roof_F3,
            'F4': create_roof_F4,
            'F5': create_roof_F5
        }
        self.t_ins = t_ins
        self.T_ext = T_ext
        self.T_int = T_int
        self.v_wind = v_wind

    def __call__(
        self,
        ID: str,
        t_ins: Quantity | None = None,
        T_ext: Quantity | None = None,
        T_int: Quantity | None = None,
        v_wind: Quantity | None = None
    ) -> ConstructionAssembly:
        """Creates the construction assembly of the roof indicated by `ID`.
        `ID` refers to the sheet number in the WTCB catalog
        (see /docs/wtcb_catalog/wtcb_catalog.pdf).

        Parameters
        ----------
        ID:
            Sheet number of the roof in the WTCB catalog.
        t_ins:
            Thickness of the insulation layer. Overrides the value assigned on
            instantiation of the `RoofCatalog` class.
        T_ext:
            The outdoor air design temperature. Overrides the value assigned on
            instantiation of the `RoofCatalog` class.
        T_int:
            The indoor air design temperature. Overrides the value assigned on
            instantiation of the `RoofCatalog` class.
        v_wind:
            Design value for the wind speed. Overrides the value assigned on
            instantiation of the `RoofCatalog` class.
        """
        if ID in self._d.keys():
            t_ins = t_ins if t_ins is not None else self.t_ins
            T_ext = T_ext if T_ext is not None else self.T_ext
            T_int = T_int if T_int is not None else self.T_int
            v_wind = v_wind if v_wind is not None else self.v_wind
            return self._d[ID](t_ins, T_ext, T_int, v_wind)
        else:
            raise KeyError('ID unknown')


class CeilingCatalog:
    """Class that bundles the functions to create construction assemblies of
    ceilings.
    """
    def __init__(
        self,
        t_ins: Quantity = Q_(10, 'cm'),
        T_adj: Quantity = Q_(15, 'degC'),
        T_int: Quantity = Q_(20, 'degC')
    ) -> None:
        """Creates an instance of `RoofCatalog`.

        Parameters
        ----------
        t_ins:
            Thickness of the insulation layer.
        T_adj:
            The air design temperature in the adjacent space.
        T_int:
            The indoor air design temperature.
        """
        self._d = {
            'F6': create_ceiling_F6,
            'F7': create_ceiling_F7,
            'F8': create_ceiling_F8,
            'F9': create_ceiling_F9,
            'F10': create_ceiling_F10,
            'F11': create_ceiling_F11,
            'F12': create_ceiling_F12,
            'F13': create_ceiling_F13,
            'F14': create_ceiling_F14
        }
        self.t_ins = t_ins
        self.T_adj = T_adj
        self.T_int = T_int

    def __call__(
        self,
        ID: str,
        t_ins: Quantity | None = None,
        T_adj: Quantity | None = None,
        T_int: Quantity | None = None
    ) -> ConstructionAssembly:
        """Creates the construction assembly of the roof indicated by `ID`.
        `ID` refers to the sheet number in the WTCB catalog
        (see /docs/wtcb_catalog/wtcb_catalog.pdf).

        Parameters
        ----------
        ID:
            Sheet number of the roof in the WTCB catalog.
        t_ins:
            Thickness of the insulation layer. Overrides the value assigned on
            instantiation of the `CeilingCatalog` class.
        T_adj:
            The air design temperature in the adjacent space. Overrides the
            value assigned on instantiation of the `CeilingCatalog` class.
        T_int:
            The indoor air design temperature. Overrides the value assigned on
            instantiation of the `CeilingCatalog` class.
        """
        if ID in self._d.keys():
            t_ins = t_ins if t_ins is not None else self.t_ins
            T_adj = T_adj if T_adj is not None else self.T_adj
            T_int = T_int if T_int is not None else self.T_int
            return self._d[ID](t_ins, T_adj, T_int)
        else:
            raise KeyError('ID unknown')


if __name__ == '__main__':
    roof_catalog = RoofCatalog(T_ext=Q_(0, 'degC'), T_int=Q_(20, 'degC'))
    roof_F1 = roof_catalog('F1', t_ins=Q_(12, 'cm'))
    print(roof_F1)

    ceiling_catalog = CeilingCatalog(t_ins=Q_(18, 'cm'))
    ceiling_F13 = ceiling_catalog('F13')
    print(ceiling_F13)
