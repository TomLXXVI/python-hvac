"""
02.D. CONSTRUCTION ASSEMBLIES: INTERIOR WALLS
Creates construction assemblies and stores them on the construction assembly
shelf.
"""
import pandas as pd
from hvac import Quantity
from hvac.cooling_load_calc.core import (
    Geometry,
    HeatFlowDirection,
    SurfaceFilm,
    SolidLayer,
    MechanicalFastening,
    ConstructionAssembly
)

from hvac.cooling_load_calc.wtcb.setup import (
    MaterialShelf,
    ConstructionAssemblyShelf,
    db_path
)

Q_ = Quantity


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F1(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj,
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    wall = SolidLayer.create(
        ID='wall',
        geometry=Geometry(t=Q_(9, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F1_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            insulation,
            wall,
            gypsum_layer,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return int_wall


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F2
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F2(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj,
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    wall = SolidLayer.create(
        ID='wall',
        geometry=Geometry(t=Q_(14, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F2_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            insulation,
            wall,
            gypsum_layer,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return int_wall


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F3
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F3(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    wall = SolidLayer.create(
        ID='wall',
        geometry=Geometry(t=Q_(10, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F3_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            insulation,
            wall,
            gypsum_layer,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return int_wall


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F4
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F4(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj,
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    wall = SolidLayer.create(
        ID='wall',
        geometry=Geometry(t=Q_(20, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F4_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            insulation,
            wall,
            gypsum_layer,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return int_wall


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F5
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F5(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    wall = SolidLayer.create(
        ID='wall',
        geometry=Geometry(t=Q_(9.0, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F5_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            insulation,
            wall,
            gypsum_layer,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return int_wall


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F6
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F6(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    wall = SolidLayer.create(
        ID='wall',
        geometry=Geometry(t=Q_(14, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F6_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            insulation,
            wall,
            gypsum_layer,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return int_wall


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F7
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F7(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    wall = SolidLayer.create(
        ID='wall',
        geometry=Geometry(t=Q_(14, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F7_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            insulation,
            wall,
            gypsum_layer,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return int_wall


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F8
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F8(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj
    )
    wall = SolidLayer.create(
        ID='wall',
        geometry=Geometry(t=Q_(14, 'cm')),
        material=MaterialShelf.load('expanded-clay-solid-block')
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F8_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            wall,
            insulation,
            gypsum_board,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return int_wall


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F9
# ------------------------------------------------------------------------------

def create_int_wall_wtcb_F9(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_adj
    )
    gypsum_board_01 = SolidLayer.create(
        ID='gypsum_board_01',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    insulation_ = SolidLayer.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.85, 'm ** 2')),
        material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
    )
    wood_ = SolidLayer.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.15, 'm ** 2')),
        material=MaterialShelf.load('wood-pine')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board_02 = SolidLayer.create(
        ID='gypsum_board_02',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialShelf.load('gypsum-cardboard')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_int
    )
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F9_t_ins={t_ins.to("cm"):~P.0f}',
        layers=[
            adj_surf_film,
            gypsum_board_01,
            insulation,
            gypsum_board_02,
            int_surf_film
        ]
    )
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return int_wall


# ------------------------------------------------------------------------------

def main():
    t_ins = Q_(6, 'cm')

    ca_int_wall_wtcb_F1 = create_int_wall_wtcb_F1(t_ins)
    ca_int_wall_wtcb_F2 = create_int_wall_wtcb_F2(t_ins)
    ca_int_wall_wtcb_F3 = create_int_wall_wtcb_F3(t_ins)
    ca_int_wall_wtcb_F4 = create_int_wall_wtcb_F4(t_ins)
    ca_int_wall_wtcb_F5 = create_int_wall_wtcb_F5(t_ins)
    ca_int_wall_wtcb_F6 = create_int_wall_wtcb_F6(t_ins)
    ca_int_wall_wtcb_F7 = create_int_wall_wtcb_F7(t_ins)
    ca_int_wall_wtcb_F8 = create_int_wall_wtcb_F8(t_ins)
    ca_int_wall_wtcb_F9 = create_int_wall_wtcb_F9(t_ins)

    ConstructionAssemblyShelf.add(
        ca_int_wall_wtcb_F1,
        ca_int_wall_wtcb_F2,
        ca_int_wall_wtcb_F3,
        ca_int_wall_wtcb_F4,
        ca_int_wall_wtcb_F5,
        ca_int_wall_wtcb_F6,
        ca_int_wall_wtcb_F7,
        ca_int_wall_wtcb_F8,
        ca_int_wall_wtcb_F9
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
    