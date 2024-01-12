"""
EXAMPLE 3
---------
CREATING A CONSTRUCTION ASSEMBLY
"""
from hvac import Quantity
from hvac.cooling_load_calc.core.construction_assembly import (
    Material,
    Geometry,
    HeatFlowDirection,
    SolidLayer,
    AirLayer,
    SurfaceFilm,
    ConstructionAssembly
)
from hvac.cooling_load_calc.core.utils import AirLayerTemperatureSolver


Q_ = Quantity


def create_construction_assembly(
    T_ext: Quantity,
    T_zone: Quantity
) -> ConstructionAssembly:
    """Creates the construction assembly of an exterior wall.

    Parameters
    ----------
    T_ext:
        Exterior temperature
    T_zone:
        Zone air temperature
    """
    ext_film = SurfaceFilm.create(
        ID='exterior_surface_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_ext,
        is_internal_surf=False,
        wind_speed=Q_(4, 'm / s')
    )
    outer_leaf = SolidLayer.create(
        ID='brick_wall_15cm',
        geometry=Geometry(t=Q_(15, 'cm')),
        material=Material(
            k=Q_(1.59, 'W / (m * K)'),
            rho=Q_(1500, 'kg / m**3'),
            c=Q_(840, 'J / (kg * K)')
        )
    )
    inner_leaf = SolidLayer.create(
        ID='brick_wall_8cm',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=Material(
            k=Q_(1.59, 'W / (m * K)'),
            rho=Q_(1500, 'kg / m**3'),
            c=Q_(840, 'J / (kg * K)')
        )
    )
    gypsum_layer = SolidLayer.create(
        ID='gypsum_1.5cm',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=Material(
            k=Q_(0.56, 'W / (m * K)'),
            rho=Q_(1300, 'kg / m**3'),
            c=Q_(840, 'J / (kg * K)')
        )
    )
    int_film = SurfaceFilm.create(
        ID='interior_surface_film',
        geometry=Geometry(),
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_zone
    )
    # Between the outer leaf and the inner leaf is an air layer. First we will
    # determine the temperature difference across the air gap and the average
    # temperature in the air gap, as we need these values to determine the
    # unit thermal resistance of the air layer.
    air_layer_temp_solver = AirLayerTemperatureSolver(
        T_ext=T_ext,
        T_zone=T_zone,
        R_ea=ext_film.R + outer_leaf.R,
        R_az=inner_leaf.R + gypsum_layer.R + int_film.R
    )
    *_, dT, T_mn = air_layer_temp_solver.solve(
        T_ae_guess=T_ext.to('K') - Q_(3, 'K'),
        T_ai_guess=T_zone.to('K') + Q_(5, 'K')
    )
    print(f"dT air layer = {dT.to('K'):~P.3f}")
    print(f"T_avg air layer = {T_mn.to('degC'):~P.3f}")

    air_layer = AirLayer.create(
        ID='air_layer',
        geometry=Geometry(t=Q_(5, 'cm'), w=Q_(float('inf'), 'm')),
        # We set the width of the air layer to infinity, to distinguish it from
        # an air void (an air void is an air layer with a width that is less
        # than 10x the air layer thickness).
        dT=dT,
        heat_flow_dir=HeatFlowDirection.HORIZONTAL,
        T_mn=T_mn,
        surf_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        angle=Q_(90, 'deg')
    )

    constr_assem = ConstructionAssembly.create(
        ID='exterior_wall',
        # The layers must always be ordered from the exterior towards the
        # interior:
        layers=[
            ext_film,
            outer_leaf,
            air_layer,
            inner_leaf,
            gypsum_layer,
            int_film
        ]
    )
    return constr_assem


def main():
    # Create the construction assembly of an exterior wall:
    ext_wall = create_construction_assembly(
        T_ext=Q_(60, 'degC'),
        T_zone=Q_(26, 'degC')
    )

    # Unit thermal resistance of exterior wall:
    print(f"exterior wall: R = {ext_wall.R.to('K * m**2 / W'):~P.3f}")

    # Unit thermal resistance of each layer in the construction assembly:
    for ID, layer in ext_wall.layers.items():
        print(f"{ID}: R = {layer.R.to('K * m**2 / W'):~P.3f}")


if __name__ == '__main__':
    main()
