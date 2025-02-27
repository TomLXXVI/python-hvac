"""02.C. CONSTRUCTION ASSEMBLIES: FLOORS

Collection of functions that create the construction assemblies of floors
contained in the WTCB catalog (see /docs/wtcb_catalog/wtcb_catalog.pdf).

Notes
-----
The function names all end with the floor ID number used in the WTCB catalog.
"""
from hvac import Quantity
from ....cooling_load_calc import (
    Geometry,
    HeatFlowDirection,
    SurfaceFilm,
    SolidLayer,
    ConstructionAssembly
)
from .setup import MaterialShelf


Q_ = Quantity


# ------------------------------------------------------------------------------
# FLOOR CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_floor_F1(
    t_ins: Quantity,
    T_adj: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.DOWNWARDS
        if T_adj < T_int
        else HeatFlowDirection.UPWARDS
    )
    ext_surf_film = SurfaceFilm.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj,
        is_internal_surf=False,
        wind_speed=v_wind
    )
    floor_slabs = SolidLayer.create(
        ID='floor_slabs',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('precast-slab-heavy-concrete-t=12cm')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polystyrene-extruded-sheet')
    )
    screed_concrete = SolidLayer.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=MaterialShelf.load('concrete-light-1600kg/m3')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    floor = ConstructionAssembly.create(
        ID=f'floor_wtcb_F1',
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
    floor.layers['floor_slabs'].num_slices = 10
    floor.layers['insulation'].num_slices = 5
    floor.layers['screed_concrete'].num_slices = 10
    return floor


# ------------------------------------------------------------------------------
# FLOOR CONSTRUCTION ASSEMBLY WTCB F2
# ------------------------------------------------------------------------------

def create_floor_F2(
    t_ins: Quantity,
    T_adj: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.DOWNWARDS
        if T_adj < T_int
        else HeatFlowDirection.UPWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj
    )
    concrete_slab = SolidLayer.create(
        ID='concrete_slab',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('concrete-reinforced-2%-steel')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polystyrene-extruded-sheet')
    )
    screed_concrete = SolidLayer.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=MaterialShelf.load('concrete-light-1600kg/m3')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    floor = ConstructionAssembly.create(
        ID=f'floor_wtcb_F2',
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
    floor.layers['concrete_slab'].num_slices = 10
    floor.layers['insulation'].num_slices = 5
    floor.layers['screed_concrete'].num_slices = 10
    return floor


# ------------------------------------------------------------------------------
# FLOOR CONSTRUCTION ASSEMBLY WTCB F3
# ------------------------------------------------------------------------------

def create_floor_F3(
    t_ins: Quantity,
    T_adj: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.DOWNWARDS
        if T_adj < T_int
        else HeatFlowDirection.UPWARDS
    )
    adj_surf_film = SurfaceFilm.create(
        ID='adj_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_adj,
    )
    floor_slab = SolidLayer.create(
        ID='floor_slab',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('precast-slab-heavy-concrete-t=12cm')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polystyrene-extruded-sheet')
    )
    screed_concrete = SolidLayer.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=MaterialShelf.load('concrete-light-1600kg/m3')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    floor = ConstructionAssembly.create(
        ID=f'floor_wtcb_F3',
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
    floor.layers['floor_slab'].num_slices = 10
    floor.layers['insulation'].num_slices = 5
    floor.layers['screed_concrete'].num_slices = 10
    return floor


# ------------------------------------------------------------------------------
# FLOOR CONSTRUCTION ASSEMBLY WTCB F4
# ------------------------------------------------------------------------------

def create_floor_F4(
    t_ins: Quantity,
    T_grd: Quantity,
    T_int: Quantity = Q_(20.0, 'degC')
) -> ConstructionAssembly:
    heat_flow_dir = (
        HeatFlowDirection.DOWNWARDS
        if T_grd < T_int
        else HeatFlowDirection.UPWARDS
    )
    floor_slab = SolidLayer.create(
        ID='floor_slab',
        geometry=Geometry(t=Q_(12, 'cm')),
        material=MaterialShelf.load('concrete-reinforced-2%-steel')
    )
    insulation = SolidLayer.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialShelf.load('polystyrene-extruded-sheet')
    )
    screed_concrete = SolidLayer.create(
        ID='screed_concrete',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=MaterialShelf.load('concrete-light-1600kg/m3')
    )
    int_surf_film = SurfaceFilm.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_dir=heat_flow_dir,
        T_mn=T_int
    )
    floor = ConstructionAssembly.create(
        ID=f'floor_wtcb_F4',
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
class FloorCatalog:
    """Class that bundles the functions to create construction assemblies of
    floors.
    """
    def __init__(
        self,
        t_ins: Quantity = Q_(10, 'cm'),
        T_adj: Quantity = Q_(0, 'degC'),
        T_int: Quantity = Q_(20, 'degC'),
        v_wind: Quantity = Q_(4, 'm / s'),
        T_grd: Quantity = Q_(20, 'degC')
    ) -> None:
        """Creates an instance of `FloorCatalog`.

        Parameters
        ----------
        t_ins:
            Thickness of the insulation layer.
        T_adj:
            The design air temperature in the adjacent space.
        T_int:
            The indoor air design temperature.
        v_wind:
            Design value for the wind speed.
        """
        self._d = {
            'F1': create_floor_F1,
            'F2': create_floor_F2,
            'F3': create_floor_F3,
            'F4': create_floor_F4,
        }
        self.t_ins = t_ins
        self.T_adj = T_adj
        self.T_int = T_int
        self.v_wind = v_wind
        self.T_grd = T_grd

    def __call__(
        self,
        ID: str,
        t_ins: Quantity | None = None,
        T_adj: Quantity | None = None,
        T_int: Quantity | None = None,
        v_wind: Quantity | None = None,
        T_grd: Quantity | None = None
    ) -> ConstructionAssembly:
        """Creates the construction assembly of the floor indicated by `ID`.
        `ID` refers to the sheet number in the WTCB catalog
        (see /docs/wtcb_catalog/wtcb_catalog.pdf).

        Parameters
        ----------
        ID:
            Sheet number of the floor in the WTCB catalog.
        t_ins:
            Thickness of the insulation layer. Overrides the value assigned on
            instantiation of the `FloorCatalog` class.
        T_adj:
            The outdoor air design temperature. Overrides the value assigned on
            instantiation of the `FloorCatalog` class.
        T_int:
            The indoor air design temperature. Overrides the value assigned on
            instantiation of the `FloorCatalog` class.
        v_wind:
            Design value for the wind speed. Overrides the value assigned on
            instantiation of the `FloorCatalog` class. Only for floor
            with ID 'F1'.
        T_grd:
            Design temperature of the ground under the floor. Overrides the
            value assigned on instantiation of the `FloorCatalog` class. Only
            for floor with ID 'F4'.
        """
        if ID in self._d.keys():
            t_ins = t_ins if t_ins is not None else self.t_ins
            T_adj = T_adj if T_adj is not None else self.T_adj
            T_int = T_int if T_int is not None else self.T_int
            v_wind = v_wind if v_wind is not None else self.v_wind
            T_grd = T_grd if T_grd is not None else self.T_grd
            if ID == 'F1':
                return self._d[ID](t_ins, T_adj, T_int, v_wind)
            elif ID == 'F4':
                return self._d[ID](t_ins, T_grd, T_int)
            else:
                return self._d[ID](t_ins, T_adj, T_int)
        else:
            raise KeyError('ID unknown')


if __name__ == '__main__':
    catalog = FloorCatalog(T_adj=Q_(32, 'degC'), T_int=Q_(26, 'degC'))
    floor_F4 = catalog('F4', t_ins=Q_(12, 'cm'))
    print(floor_F4)
