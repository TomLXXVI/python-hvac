"""02.D. CONSTRUCTION ASSEMBLIES: INTERIOR WALLS

Collection of functions that create the construction assemblies of interior
walls contained in the WTCB catalog (see /docs/wtcb_catalog/wtcb_catalog.pdf).

Notes
-----
The function names all end with the wall ID number used in the WTCB catalog.
"""
from hvac import Quantity
from ....cooling_load_calc import (
    Geometry,
    HeatFlowDirection,
    SurfaceFilm,
    SolidLayer,
    MechanicalFastening,
    ConstructionAssembly
)
from .setup import MaterialShelf

Q_ = Quantity


# ------------------------------------------------------------------------------
# INTERIOR WALL CONSTRUCTION ASSEMBLY WTCB F1
# ------------------------------------------------------------------------------

def create_int_wall_F1(
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
    if t_ins:
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
        )
    else:
        insulation = None
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
    if insulation is not None:
        layers = [adj_surf_film, insulation, wall, gypsum_layer, int_surf_film]
    else:
        layers = [adj_surf_film, wall, gypsum_layer, int_surf_film]
    int_wall = ConstructionAssembly.create('int_wall_wtcb_F1', layers)
    if insulation is not None:
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

def create_int_wall_F2(
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
    if t_ins:
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
        )
    else:
        insulation = None
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
    if insulation is not None:
        layers = [adj_surf_film, insulation, wall, gypsum_layer, int_surf_film]
    else:
        layers = [adj_surf_film, wall, gypsum_layer, int_surf_film]
    int_wall = ConstructionAssembly.create('int_wall_wtcb_F2', layers)
    if insulation is not None:
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

def create_int_wall_F3(
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
    if t_ins:
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
        )
    else:
        insulation = None
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
    if insulation is not None:
        layers = [adj_surf_film, insulation, wall, gypsum_layer, int_surf_film]
    else:
        layers = [adj_surf_film, wall, gypsum_layer, int_surf_film]
    int_wall = ConstructionAssembly.create('int_wall_wtcb_F3', layers)
    if insulation is not None:
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

def create_int_wall_F4(
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
    if t_ins:
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
        )
    else:
        insulation = None
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
    if insulation is not None:
        layers = [adj_surf_film, insulation, wall, gypsum_layer, int_surf_film]
    else:
        layers = [adj_surf_film, wall, gypsum_layer, int_surf_film]
    int_wall = ConstructionAssembly.create('int_wall_wtcb_F4', layers)
    if insulation is not None:
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

def create_int_wall_F5(
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
    if t_ins:
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
        )
    else:
        insulation = None
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
    if insulation is not None:
        layers = [adj_surf_film, insulation, wall, gypsum_layer, int_surf_film]
    else:
        layers = [adj_surf_film, wall, gypsum_layer, int_surf_film]
    int_wall = ConstructionAssembly.create('int_wall_wtcb_F5', layers)
    if insulation is not None:
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

def create_int_wall_F6(
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
    if t_ins:
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
        )
    else:
        insulation = None
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
    if insulation is not None:
        layers = [adj_surf_film, insulation, wall, gypsum_layer, int_surf_film]
    else:
        layers = [adj_surf_film, wall, gypsum_layer, int_surf_film]
    int_wall = ConstructionAssembly.create('int_wall_wtcb_F6', layers)
    if insulation is not None:
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

def create_int_wall_F7(
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
    if t_ins:
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
        )
    else:
        insulation = None
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
    if insulation is not None:
        layers = [adj_surf_film, insulation, wall, gypsum_layer, int_surf_film]
    else:
        layers = [adj_surf_film, wall, gypsum_layer, int_surf_film]
    int_wall = ConstructionAssembly.create(
        ID=f'int_wall_wtcb_F7',
        layers=layers
    )
    if insulation is not None:
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

def create_int_wall_F8(
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
    if t_ins:
        insulation = SolidLayer.create(
            ID='insulation',
            geometry=Geometry(t=t_ins),
            material=MaterialShelf.load('mineral-wool-glass-fibre-cover-sheet-22kg/m3')
        )
    else:
        insulation = None
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
    if insulation is not None:
        layers = [adj_surf_film, wall, insulation, gypsum_board, int_surf_film]
    else:
        layers = [adj_surf_film, wall, gypsum_board, int_surf_film]
    int_wall = ConstructionAssembly.create('int_wall_wtcb_F8', layers)
    if insulation is not None:
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

def create_int_wall_F9(
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
    layers = [adj_surf_film, gypsum_board_01, insulation, gypsum_board_02, int_surf_film]
    int_wall = ConstructionAssembly.create('int_wall_wtcb_F9', layers)
    int_wall = int_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=None
    )
    return int_wall


# ------------------------------------------------------------------------------
class InteriorWallCatalog:
    """Class that bundles the functions to create construction assemblies of
    interior walls.
    """
    def __init__(
        self,
        t_ins: Quantity = Q_(10, 'cm'),
        T_adj: Quantity = Q_(0, 'degC'),
        T_int: Quantity = Q_(20, 'degC'),
    ) -> None:
        """Creates an instance of `InteriorWallCatalog`.

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
            'F1': create_int_wall_F1,
            'F2': create_int_wall_F2,
            'F3': create_int_wall_F3,
            'F4': create_int_wall_F4,
            'F5': create_int_wall_F5,
            'F6': create_int_wall_F6,
            'F7': create_int_wall_F7,
            'F8': create_int_wall_F8,
            'F9': create_int_wall_F9
        }
        self.t_ins = t_ins
        self.T_adj = T_adj
        self.T_int = T_int

    def __call__(
        self,
        ID: str,
        t_ins: Quantity | None = None,
        T_adj: Quantity | None = None,
        T_int: Quantity | None = None,
    ) -> ConstructionAssembly:
        """Creates the construction assembly of the interior wall indicated by
        `ID`. `ID` refers to the sheet number in the WTCB catalog
        (see /docs/wtcb_catalog/wtcb_catalog.pdf).

        Parameters
        ----------
        ID:
            Sheet number of the exterior wall in the WTCB catalog.
        t_ins:
            Thickness of the insulation layer. Overrides the value assigned on
            instantiation of the `InteriorWallCatalog` class.
        T_adj:
            The air design temperature of the adjacent space. Overrides the
            value assigned on instantiation of the `InteriorWallCatalog` class.
        T_int:
            The indoor air design temperature. Overrides the value assigned on
            instantiation of the `InteriorWallCatalog` class.
        """
        if ID in self._d.keys():
            t_ins = t_ins if t_ins is not None else self.t_ins
            T_adj = T_adj if T_adj is not None else self.T_adj
            T_int = T_int if T_int is not None else self.T_int
            return self._d[ID](t_ins, T_adj, T_int)
        else:
            raise KeyError('ID unknown')


if __name__ == '__main__':
    catalog = InteriorWallCatalog(T_adj=Q_(26, 'degC'), T_int=Q_(24, 'degC'))
    iw = catalog('F9', t_ins=Q_(8, 'cm'))
    print(iw)
