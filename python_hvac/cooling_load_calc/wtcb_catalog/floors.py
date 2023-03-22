"""
02.C. CONSTRUCTION ASSEMBLIES: FLOORS

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
    ConstructionAssembly
)

from hvac.cooling_load_calc.wtcb_catalog import (
    MaterialsShelf,
    ConstructionAssembliesShelf,
    db_path
)

Q_ = Quantity


# ------------------------------------------------------------------------------
# FLOOR CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_ca_floor_wtcbF1(
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
    floor_slabs = BuildingComponent.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('welfsel, zwaar beton, t=12 cm')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('geëxtrudeerd polystyreen (XPS), plaat')
    )
    screed_concrete = BuildingComponent.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=MaterialsShelf.load('licht beton, 1600 kg/m3')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    floor = ConstructionAssembly.create(
        ID=f'floor_wtcbF1 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            floor_slabs,
            insulation,
            screed_concrete,
            int_surf_film
        ]
    )
    floor = floor.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    floor.layers['floor_slabs'].slices = 10
    floor.layers['insulation'].slices = 5
    floor.layers['screed_concrete'].slices = 10
    return floor


# ------------------------------------------------------------------------------
# FLOOR CONSTRUCTION ASSEMBLY WTCB F2
# ------------------------------------------------------------------------------

def create_ca_floor_wtcbF2(
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
    concrete_slab = BuildingComponent.create(
        ID='concrete_slab',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('beton, gewapend (2 % staal)')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('geëxtrudeerd polystyreen (XPS), plaat')
    )
    screed_concrete = BuildingComponent.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=MaterialsShelf.load('licht beton, 1600 kg/m3')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    floor = ConstructionAssembly.create(
        ID=f'floor_wtcbF2 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            concrete_slab,
            insulation,
            screed_concrete,
            int_surf_film
        ]
    )
    floor = floor.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    floor.layers['concrete_slab'].slices = 10
    floor.layers['insulation'].slices = 5
    floor.layers['screed_concrete'].slices = 10
    return floor


# ------------------------------------------------------------------------------
# FLOOR CONSTRUCTION ASSEMBLY WTCB F3
# ------------------------------------------------------------------------------

def create_ca_floor_wtcbF3(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
) -> ConstructionAssembly:
    adj_surf_film = SurfaceLayer.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_ext,
    )
    floor_slab = BuildingComponent.create(
        ID='floor_slab',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('welfsel, zwaar beton, t=12 cm')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('geëxtrudeerd polystyreen (XPS), plaat')
    )
    screed_concrete = BuildingComponent.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=MaterialsShelf.load('licht beton, 1600 kg/m3')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    floor = ConstructionAssembly.create(
        ID=f'floor_wtcbF3 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            adj_surf_film,
            floor_slab,
            insulation,
            screed_concrete,
            int_surf_film
        ]
    )
    floor = floor.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    floor.layers['floor_slab'].slices = 10
    floor.layers['insulation'].slices = 5
    floor.layers['screed_concrete'].slices = 10
    return floor


# ------------------------------------------------------------------------------
# FLOOR CONSTRUCTION ASSEMBLY WTCB F4
# ------------------------------------------------------------------------------

def create_ca_floor_wtcbF4(
    t_ins: Quantity,
    heat_flow_direction: HeatFlowDirection,
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    floor_slab = BuildingComponent.create(
        ID='floor_slab',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialsShelf.load('beton, gewapend (2 % staal)')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('geëxtrudeerd polystyreen (XPS), plaat')
    )
    screed_concrete = BuildingComponent.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=MaterialsShelf.load('licht beton, 1600 kg/m3')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=heat_flow_direction,
        Tmn=T_int
    )
    floor = ConstructionAssembly.create(
        ID=f'floor_wtcbF4 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            floor_slab,
            insulation,
            screed_concrete,
            int_surf_film
        ]
    )
    floor = floor.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return floor


# ------------------------------------------------------------------------------

def main():
    t_ins = Q_(12, 'cm')
    hfd = HeatFlowDirection.UPWARDS

    ca_floor_wtcbF1 = create_ca_floor_wtcbF1(t_ins, hfd)
    ca_floor_wtcbF2 = create_ca_floor_wtcbF2(t_ins, hfd)
    ca_floor_wtcbF3 = create_ca_floor_wtcbF3(t_ins, hfd)
    ca_floor_wtcbF4 = create_ca_floor_wtcbF4(t_ins, hfd)

    ConstructionAssembliesShelf.add(
        ca_floor_wtcbF1,
        ca_floor_wtcbF2,
        ca_floor_wtcbF3,
        ca_floor_wtcbF4
    )

    # print(ca_floor_wtcbF1, end='\n\n')
    # print(ca_floor_wtcbF2, end='\n\n')
    # print(ca_floor_wtcbF3, end='\n\n')
    # print(ca_floor_wtcbF4, end='\n\n')

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
