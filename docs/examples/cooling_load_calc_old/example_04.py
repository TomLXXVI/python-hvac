"""
EXAMPLE 4
---------
CREATING AN EXTERIOR BUILDING ELEMENT - CONDUCTION HEAT GAIN
"""
import pandas as pd
from hvac import Quantity
from hvac.sun import Location, ClimateType, ReferenceDates
from hvac.cooling_load_calc_old.core import (
    Material,
    Geometry,
    HeatFlowDirection,
    SolidLayer,
    AirLayer,
    SurfaceFilm,
    ConstructionAssembly,
    WeatherData,
    ExteriorBuildingElement
)
from hvac.cooling_load_calc_old.core.utils import (
    AirLayerTemperatureSolver,
    convert_to_clock_time
)
from hvac.charts import LineChart

Q_ = Quantity


def create_construction_assembly(
    T_ext: Quantity,
    T_zone: Quantity
) -> ConstructionAssembly:
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
    outer_leaf.num_slices = 5
    # this will determine the number of temperature nodes of the outer leaf in
    # the linear thermal network of the exterior building element
    inner_leaf = SolidLayer.create(
        ID='brick_wall_8cm',
        geometry=Geometry(t=Q_(8, 'cm')),
        material=Material(
            k=Q_(1.59, 'W / (m * K)'),
            rho=Q_(1500, 'kg / m**3'),
            c=Q_(840, 'J / (kg * K)')
        )
    )
    inner_leaf.num_slices = 4
    gypsum_layer = SolidLayer.create(
        ID='gypsum_1.5cm',
        geometry=Geometry(t=Q_(2.0, 'cm')),
        material=Material(
            k=Q_(0.56, 'W / (m * K)'),
            rho=Q_(1300, 'kg / m**3'),
            c=Q_(840, 'J / (kg * K)')
        )
    )
    gypsum_layer.num_slices = 2
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
        T_int=T_zone,
        R_ea=ext_film.R + outer_leaf.R,
        R_ai=inner_leaf.R + gypsum_layer.R + int_film.R
    )
    *_, dT, T_mn = air_layer_temp_solver.solve(
        T_ae_guess=T_ext.to('K') - Q_(3, 'K'),
        T_ai_guess=T_zone.to('K') + Q_(5, 'K')
    )
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
    # CREATING AN EXTERIOR BUILDING ELEMENT
    # Set the geographic location:
    location = Location(
        fi=Q_(51.183, 'deg'),
        L_loc=Q_(3.8, 'deg'),
        altitude=Q_(8, 'm'),
        climate_type=ClimateType.MID_LATITUDE_SUMMER,
        timezone='Etc/GMT-1'
    )
    # Create a weather data model for the given location based on the clear-sky
    # model for solar radiation:
    weather_data = WeatherData.create_from_climatic_design_data(
        location=location,
        date=ReferenceDates.get_date_for('Jun'),  # this is the "design day"
        T_db_des=Q_(24.2, 'degC'),
        T_db_rng=Q_(11.8, 'K'),
        T_wb_mc=Q_(18.0, 'degC'),
        T_wb_rng=Q_(5.6, 'K')
    )
    # Create the construction assembly of the exterior building element:
    constr_assem_ext_wall = create_construction_assembly(
        T_ext=weather_data.T_db(12),
        T_zone=Q_(22, 'degC')
    )
    # Save the construction assembly to a shelf on disk for later reuse.
    # --> First we need to set the path to the shelf:
    ConstructionAssembly.db_path = "./example_04_constr_assem.db"
    # --> To save the construction assembly on the shelf, call method `save()` on
    #   the `ConstructionAssembly` object:
    constr_assem_ext_wall.save()
    # Create the exterior building element:
    ext_wall = ExteriorBuildingElement.create(
        ID='exterior-south-wall',
        T_zone=lambda t_sol_sec: Q_(22, 'degC'),
        constr_assem=constr_assem_ext_wall,
        gross_area=Q_(10, 'm') * Q_(3, 'm'),  # surface area of the wall
        weather_data=weather_data,
        gamma=Q_(0, 'deg'),  # oriented to the south
        beta=Q_(90, 'deg')   # vertical wall
    )

    # HEAT CONDUCTION THROUGH THE EXTERIOR BUILDING ELEMENT:
    # Solving for the node temperatures and the heat flows in the linear thermal
    # network model of the exterior building element:
    tnt, qdt = ext_wall.solve(num_cycles=5)
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print('NODE TEMPERATURE TABLE:')
        print(tnt)
        print()
        print('HEAT FLOW TABLE:')
        print(qdt)
        print()

    # INTERIOR CONDUCTION HEAT GAIN
    # The conduction heat gain at the interior side of the exterior building
    # element has a convective and a radiative component.
    Q_dot_cond, Q_dot_conv, Q_dot_rad = ext_wall.conductive_heat_gain(num_cycles=5)

    # Display a table with the interior conduction heat gains in the course of
    # time with their corresponding solar and local time at the location:
    n = len(Q_dot_cond)
    loc_time, sol_time = zip(*[
        convert_to_clock_time(
            time_index=k,
            dt_hr=1.0,
            date=weather_data.date,
            L_loc=location.L_loc,
            tz_loc=location.timezone
        ) for k in range(n)
    ])
    d = {
        'solar time': [st.time() for st in sol_time],
        f'{location.timezone}': [lt.time() for lt in loc_time],
        'Q_dot_cond': [q_dot_cond.to('W').m for q_dot_cond in Q_dot_cond]
    }
    df = pd.DataFrame(d)
    print('CONDUCTION HEAT GAIN THROUGH EXTERIOR WALL')
    print(df)

    # Plot the interior and exterior conduction heat gain as function of time:
    chart = LineChart()
    chart.add_xy_data(
        label="interior conduction heat gain",
        x1_values=[lt.hour for lt in loc_time],
        y1_values=Q_dot_cond.magnitude,
        style_props={'marker': 'o'}
    )
    chart.add_xy_data(
        label="exterior conduction heat gain",
        x1_values=[lt.hour for lt in loc_time],
        y1_values=[Q_dot for Q_dot in qdt['ESN'].values],
        style_props={'marker': 'o'}
    )
    chart.add_legend(anchor='upper right', position=(0.999, 0.999), columns=1)
    chart.x1.add_title('time index')
    chart.y1.add_title('conductive heat gain, W')
    chart.show()


if __name__ == '__main__':
    main()
