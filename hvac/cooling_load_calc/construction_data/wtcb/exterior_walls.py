"""02.A. CONSTRUCTION ASSEMBLIES: EXTERIOR WALLS

Collection of functions that create the construction assemblies of exterior
walls described in the WTCB catalog (see /docs/wtcb_catalog/wtcb_catalog.pdf).

Notes
-----
The function names all end with the wall ID number used in the WTCB catalog.
"""
from hvac import Quantity
from ....cooling_load_calc import (
    Geometry,
    Material,
    HeatFlowDirection,
    SurfaceFilm,
    SolidLayer,
    AirLayer,
    MechanicalFastening,
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
    t: Quantity,
    surf_emissivities: tuple[Quantity, Quantity] = (Q_(0.9, 'frac'), Q_(0.9, 'frac'))
) -> tuple[AirLayer, Quantity]:
    """Internal helper function for creating an `AirLayer` object.

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
    t:
        The thickness of the air gap.
    surf_emissivities:
        A tuple with the emissivities of the exterior side surface and the
        interior side surface of the air layer.

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
    else:
        dT = T_ext.to('K') - T_int.to('K')
        T_ae_ini = T_ext.to('K') - dT / 2 + Q_(2.5, 'K')
        T_ai_ini = T_int.to('K') + dT / 2 - Q_(2.5, 'K')
    *_, dT_asp, T_asp = ats.solve(T_ae_ini, T_ai_ini)
    air_space = AirLayer.create(
        ID='air_space',
        geometry=Geometry(t=t),
        dT=dT_asp,
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp,
        surf_emissivities=surf_emissivities,
        angle=Q_(0.0, 'deg')
    )
    return air_space, T_asp


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_ext_wall_F1(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext-surf-film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # airspace: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(15.0, 'cm')),
        material=MaterialShelf.load('concrete-aerated-block-glued-600kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F1',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F2
# ------------------------------------------------------------------------------

def create_ext_wall_F2(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('terracotta-block-1200kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F2',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F3
# ------------------------------------------------------------------------------

def create_ext_wall_F3(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(19.0, 'cm')),
        material=MaterialShelf.load('terracotta-block-1200kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F3',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F4
# ------------------------------------------------------------------------------

def create_ext_wall_F4(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('expanded-clay-solid-block')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F4',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F5
# ------------------------------------------------------------------------------

def create_ext_wall_F5(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('concrete-medium-solid-block-1700kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext-wall-wtcb-F5',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F6
# ------------------------------------------------------------------------------

def create_ext_wall_F6(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('expanded-clay-hollow-block-1200kg/m3-t=14cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F6',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F7
# ------------------------------------------------------------------------------

def create_ext_wall_F7(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polyurethane-sheet-150kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('expanded-clay-hollow-block-1200kg/m3-t=14cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F7',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F8
# ------------------------------------------------------------------------------

def create_ext_wall_F8(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(19.0, 'cm')),
        material=MaterialShelf.load('expanded-clay-hollow-block-1200kg/m3-t=19cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F8',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F9
# ------------------------------------------------------------------------------

def create_ext_wall_F9(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('concrete-heavy-hollow-block-1400kg/m3-t=14cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F9',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F10
# ------------------------------------------------------------------------------

def create_ext_wall_F10(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polyurethane-sheet-150kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('concrete-heavy-hollow-block-1400kg/m3-t=14cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F10',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F11
# ------------------------------------------------------------------------------

def create_ext_wall_F11(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialShelf.load('concrete-medium-solid-block-1800kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(15.0, 'cm')),
        material=MaterialShelf.load('concrete-aerated-block-glued-600kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F11',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F12
# ------------------------------------------------------------------------------

def create_ext_wall_F12(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialShelf.load('concrete-medium-solid-block-1800kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('expanded-clay-hollow-block-1200kg/m3-t=14cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F12',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F13
# ------------------------------------------------------------------------------

def create_ext_wall_F13(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialShelf.load('concrete-medium-solid-block-1800kg/m3'),
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('concrete-heavy-hollow-block-1400kg/m3-t=14cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F13',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F14
# ------------------------------------------------------------------------------

def create_ext_wall_F14(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    wood_siding = SolidLayer.create(
        ID='wood_siding',
        geometry=Geometry(t=Q_(2.5, 'cm')),
        material=Material(
            k=Q_(0.18, 'W / (m * K)'),
            rho=Q_(540, 'kg / m ** 3'),
            c=Q_(1.17, 'kJ / (kg * K)')
        )
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(15.0, 'cm')),
        material=MaterialShelf.load('concrete-aerated-block-glued-600kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + wood_siding.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F14',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    # ext_wall.layers['wood_siding'].num_slices = 10
    # --> removed from construction assembly as AirLayer is well ventilated
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F15
# ------------------------------------------------------------------------------

def create_ext_wall_F15(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    wood_siding = SolidLayer.create(
        ID='wood_siding',
        geometry=Geometry(t=Q_(2.5, 'cm')),
        material=Material(
            k=Q_(0.18, 'W / (m * K)'),
            rho=Q_(540, 'kg / m ** 3'),
            c=Q_(1.17, 'kJ / (kg * K)')
        )
    )
    # air_space: see below
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=Q_(6, 'cm')),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('expanded-clay-hollow-block-1200kg/m3-t=14cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + wood_siding.R,
        R_ai=insulation.R + inner_leaf.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F15',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    # ext_wall.layers['wood_siding'].num_slices = 10
    # --> removed from construction assembly as AirLayer is well ventilated
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F16
# ------------------------------------------------------------------------------

def create_ext_wall_F16(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3')
    )
    # air_space: see below
    OSB_board = SolidLayer.create(
        ID='OSB_board',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=MaterialShelf.load('osb-board')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + outer_leaf.R,
        R_ai=OSB_board.R + insulation.R + gypsum_board.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F16',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F17
# ------------------------------------------------------------------------------

def create_ext_wall_F17(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    wood_siding = SolidLayer.create(
        ID='wood_siding',
        geometry=Geometry(t=Q_(2.5, 'cm')),
        material=Material(
            k=Q_(0.18, 'W / (m * K)'),
            rho=Q_(540, 'kg / m ** 3'),
            c=Q_(1.17, 'kJ / (kg * K)')
        )
    )
    # air_space: see below
    OSB_board = SolidLayer.create(
        ID='OSB_board',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=MaterialShelf.load('osb-board')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    gypsum_board = SolidLayer.create(
        ID='gypsum_board',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + wood_siding.R,
        R_ai=OSB_board.R + insulation.R + gypsum_board.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F17',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    # ext_wall.layers['wood_siding'].num_slices = 10
    # --> removed from construction assembly as AirLayer is well ventilated
    ext_wall.layers['insulation'].num_slices = 5
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F18
# ------------------------------------------------------------------------------

def create_ext_wall_F18(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('terracotta-block-1200kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F18',
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F19
# ------------------------------------------------------------------------------

def create_ext_wall_F19(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(9.0, 'cm')),
        material=MaterialShelf.load('terracotta-brick-1500kg/m3')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('expanded-clay-hollow-block-1200kg/m3-t=14cm')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F19',
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
    ext_wall.layers['outer_leaf'].num_slices = 10
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F20
# ------------------------------------------------------------------------------

def create_ext_wall_F20(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    cement_plastering = SolidLayer.create(
        ID='cement_plastering',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=MaterialShelf.load('cement-plaster')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('concrete-aerated-block-glued-600kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F20',
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
    ext_wall.layers['cement_plastering'].num_slices = 2
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F21
# ------------------------------------------------------------------------------

def create_ext_wall_F21(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    cement_plastering = SolidLayer.create(
        ID='cement_plastering',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=MaterialShelf.load('cement-plaster')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    inner_leaf = SolidLayer.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(14.0, 'cm')),
        material=MaterialShelf.load('expanded-clay-solid-block')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F21',
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
    ext_wall.layers['cement_plastering'].num_slices = 2
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['inner_leaf'].num_slices = 15
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F22
# ------------------------------------------------------------------------------

def create_ext_wall_F22(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    vinyl_siding = SolidLayer.create(
        ID='vinyl_siding',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=Material(
            k=Q_(0.16, 'W / (m * K)'),
            rho=Q_(1375.0, 'kg / m ** 3'),
            c=Q_(250.0, 'J / (kg * K)')
        )
    )
    # air_space: see below
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(29.0, 'cm')),
        material=MaterialShelf.load('terracotta-block-1500kg/m3')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polystyrene-extruded-sheet')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + vinyl_siding.R,
        R_ai=outer_leaf.R + insulation.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F22',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 15
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['gypsum_layer'].num_slices = 2
    return ext_wall


# ------------------------------------------------------------------------------
# EXTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F23
# ------------------------------------------------------------------------------

def create_ext_wall_F23(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    vinyl_siding = SolidLayer.create(
        ID='vinyl_siding',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=Material(
            k=Q_(0.16, 'W / (m * K)'),
            rho=Q_(1375.0, 'kg / m ** 3'),
            c=Q_(250.0, 'J / (kg * K)')
        )
    )
    # air_space: see below
    outer_leaf = SolidLayer.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(29.0, 'cm')),
        material=MaterialShelf.load('terracotta-block-1500kg/m3')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polyurethane-sheet-150kg/m3')
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-plaster')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    air_space, T_asp = _create_air_layer(
        T_ext=T_ext,
        T_int=T_int,
        R_ea=ext_surf_film.R + vinyl_siding.R,
        R_ai=outer_leaf.R + insulation.R + gypsum_layer.R + int_surf_film.R,
        t=Q_(5.0, 'cm')
    )
    ext_wall = ConstructionAssembly.create(
        ID='ext-wall-wtcb-F23',
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
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_asp
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
    ext_wall.layers['outer_leaf'].num_slices = 15
    ext_wall.layers['insulation'].num_slices = 5
    ext_wall.layers['gypsum_layer'].num_slices = 2
    return ext_wall


# ------------------------------------------------------------------------------
class ExteriorWallCatalog:
    """Class that bundles all the functions to create construction assemblies of
    exterior walls.
    """
    def __init__(
        self,
        t_ins: Quantity = Q_(10, 'cm'),
        T_ext: Quantity = Q_(0, 'degC'),
        T_int: Quantity = Q_(20, 'degC'),
        v_wind: Quantity = Q_(4, 'm / s')
    ) -> None:
        """Creates an instance of `ExteriorWallCatalog`.

        Parameters
        ----------
        t_ins:
            Default thickness of the insulation layer.
        T_ext:
            Default design value for outdoor air temperature.
        T_int:
            Default design value for indoor air temperature.
        v_wind:
            Default design value for wind speed.
        """
        self._d = {
            'F1': create_ext_wall_F1,
            'F2': create_ext_wall_F2,
            'F3': create_ext_wall_F3,
            'F4': create_ext_wall_F4,
            'F5': create_ext_wall_F5,
            'F6': create_ext_wall_F6,
            'F7': create_ext_wall_F7,
            'F8': create_ext_wall_F8,
            'F9': create_ext_wall_F9,
            'F10': create_ext_wall_F10,
            'F11': create_ext_wall_F11,
            'F12': create_ext_wall_F12,
            'F13': create_ext_wall_F13,
            'F14': create_ext_wall_F14,
            'F15': create_ext_wall_F15,
            'F16': create_ext_wall_F16,
            'F17': create_ext_wall_F17,
            'F18': create_ext_wall_F18,
            'F19': create_ext_wall_F19,
            'F20': create_ext_wall_F20,
            'F21': create_ext_wall_F21,
            'F22': create_ext_wall_F22,
            'F23': create_ext_wall_F23
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
        """Creates the construction assembly of the exterior wall indicated by
        `ID`. `ID` refers to the wall ID number in the WTCB catalog
        (see /docs/wtcb_catalog/wtcb_catalog.pdf).

        Parameters
        ----------
        ID:
            Sheet number of the exterior wall in the WTCB catalog.
        t_ins:
            Thickness of the insulation layer. Overrides the value assigned on
            instantiation of the `ExteriorWallCatalog` class.
        T_ext:
            The outdoor air design temperature. Overrides the value assigned on
            instantiation of the `ExteriorWallCatalog` class.
        T_int:
            The indoor air design temperature. Overrides the value assigned on
            instantiation of the `ExteriorWallCatalog` class.
        v_wind:
            Design value for the wind speed. Overrides the value assigned on
            instantiation of the `ExteriorWallCatalog` class.
        """
        if ID in self._d.keys():
            t_ins = t_ins if t_ins is not None else self.t_ins
            T_ext = T_ext if T_ext is not None else self.T_ext
            T_int = T_int if T_int is not None else self.T_int
            v_wind = v_wind if v_wind is not None else self.v_wind
            return self._d[ID](t_ins, T_ext, T_int, v_wind)
        else:
            raise KeyError('ID unknown')


if __name__ == '__main__':
    catalog = ExteriorWallCatalog(T_ext=Q_(32, 'degC'), T_int=Q_(26, 'degC'))
    exterior_wall = catalog('F1', t_ins=Q_(12, 'cm'))
    print(exterior_wall)
