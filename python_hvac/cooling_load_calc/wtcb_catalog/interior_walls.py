"""
02.D. CONSTRUCTION ASSEMBLIES: INTERIOR WALLS

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
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_ca_int_wall_wtcbF1(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj,
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    wall = BuildingComponent.create(
        ID='wall',
        geometry=Geometry(t=Q_(9, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF1 (t_ins={t_ins.to("cm"):~P.0f})',
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

def create_ca_int_wall_wtcbF2(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj,
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    wall = BuildingComponent.create(
        ID='wall',
        geometry=Geometry(t=Q_(14, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF2 (t_ins={t_ins.to("cm"):~P.0f})',
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

def create_ca_int_wall_wtcbF3(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    wall = BuildingComponent.create(
        ID='wall',
        geometry=Geometry(t=Q_(10, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF3 (t_ins={t_ins.to("cm"):~P.0f})',
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

def create_ca_int_wall_wtcbF4(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj,
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    wall = BuildingComponent.create(
        ID='wall',
        geometry=Geometry(t=Q_(20, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF4 (t_ins={t_ins.to("cm"):~P.0f})',
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

def create_ca_int_wall_wtcbF5(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    wall = BuildingComponent.create(
        ID='wall',
        geometry=Geometry(t=Q_(9.0, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF5 (t_ins={t_ins.to("cm"):~P.0f})',
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

def create_ca_int_wall_wtcbF6(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    wall = BuildingComponent.create(
        ID='wall',
        geometry=Geometry(t=Q_(14, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF6 (t_ins={t_ins.to("cm"):~P.0f})',
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

def create_ca_int_wall_wtcbF7(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    wall = BuildingComponent.create(
        ID='wall',
        geometry=Geometry(t=Q_(14, 'cm')),
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF7 (t_ins={t_ins.to("cm"):~P.0f})',
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

def create_ca_int_wall_wtcbF8(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj
    )
    wall = BuildingComponent.create(
        ID='wall',
        geometry=Geometry(t=Q_(14, 'cm')),
        material=MaterialsShelf.load('volle betonblokken, geëxpandeerde klei')
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
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF8 (t_ins={t_ins.to("cm"):~P.0f})',
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

def create_ca_int_wall_wtcbF9(
    t_ins: Quantity,
    T_adj: Quantity = Q_(10, 'degC'),
    T_int: Quantity = Q_(20, 'degC')
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_adj
    )
    gypsum_board_01 = BuildingComponent.create(
        ID='gypsum_board_01',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    insulation_ = BuildingComponent.create(
        ID='insulation_',
        geometry=Geometry(t=t_ins, A=Q_(0.85, 'm ** 2')),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    wood_ = BuildingComponent.create(
        ID='wood_',
        geometry=Geometry(t=t_ins, A=Q_(0.15, 'm ** 2')),
        material=MaterialsShelf.load('naaldhout')
    )
    insulation = insulation_ // wood_
    insulation.ID = 'insulation'
    gypsum_board_02 = BuildingComponent.create(
        ID='gypsum_board_02',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipsplaat, tussen 2 lagen karton')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcbF9 (t_ins={t_ins.to("cm"):~P.0f})',
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

    ca_int_wall_wtcbF1 = create_ca_int_wall_wtcbF1(t_ins)
    ca_int_wall_wtcbF2 = create_ca_int_wall_wtcbF2(t_ins)
    ca_int_wall_wtcbF3 = create_ca_int_wall_wtcbF3(t_ins)
    ca_int_wall_wtcbF4 = create_ca_int_wall_wtcbF4(t_ins)
    ca_int_wall_wtcbF5 = create_ca_int_wall_wtcbF5(t_ins)
    ca_int_wall_wtcbF6 = create_ca_int_wall_wtcbF6(t_ins)
    ca_int_wall_wtcbF7 = create_ca_int_wall_wtcbF7(t_ins)
    ca_int_wall_wtcbF8 = create_ca_int_wall_wtcbF8(t_ins)
    ca_int_wall_wtcbF9 = create_ca_int_wall_wtcbF9(t_ins)

    # print(ca_int_wall_wtcbF1, end='\n\n')
    # print(ca_int_wall_wtcbF2, end='\n\n')
    # print(ca_int_wall_wtcbF3, end='\n\n')
    # print(ca_int_wall_wtcbF4, end='\n\n')
    # print(ca_int_wall_wtcbF5, end='\n\n')
    # print(ca_int_wall_wtcbF6, end='\n\n')
    # print(ca_int_wall_wtcbF7, end='\n\n')
    # print(ca_int_wall_wtcbF8, end='\n\n')
    # print(ca_int_wall_wtcbF9, end='\n\n')

    ConstructionAssembliesShelf.add(
        ca_int_wall_wtcbF1,
        ca_int_wall_wtcbF2,
        ca_int_wall_wtcbF3,
        ca_int_wall_wtcbF4,
        ca_int_wall_wtcbF5,
        ca_int_wall_wtcbF6,
        ca_int_wall_wtcbF7,
        ca_int_wall_wtcbF8,
        ca_int_wall_wtcbF9
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
