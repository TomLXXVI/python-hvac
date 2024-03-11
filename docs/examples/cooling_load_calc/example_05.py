"""
EXAMPLE 5
---------
CREATING AN EXTERIOR BUILDING ELEMENT WITH A WINDOW - CONDUCTION HEAT GAIN AND
SOLAR HEAT GAIN.
"""
from typing import Callable
import pandas as pd
from hvac import Quantity
from hvac.sun import Location, ClimateType, ReferenceDates
from hvac.cooling_load_calc import (
    WeatherData,
    ConstructionAssembly,
    ExteriorBuildingElement,
    Window
)
from hvac.cooling_load_calc.shelves import WindowPropertiesShelf
from hvac.cooling_load_calc.core.utils import convert_to_clock_time
from hvac.charts import LineChart

Q_ = Quantity


def create_exterior_wall(
    weather_data: WeatherData,
    T_zone: Callable[[float], Quantity]
) -> ExteriorBuildingElement:
    """Creates an exterior wall.

    Parameters
    ----------
    weather_data:
        `WeatherData` object that holds all the climatic design information
        needed to calculate the conduction heat gain through the exterior
        building element on the geographic location and on the design day that
        were specified on instantiation of the `WeatherData` object.
    T_zone:
        The zone air temperature, being a function with signature
        `f(t_sol_sec: float) -> Quantity` which takes the solar time in
        seconds and returns the temperature in the zone as a `Quantity`
        object. This may allow a time-variable temperature in the zone.

    Returns
    -------
    The `ExteriorBuildingElement` object that models the exterior wall.
    """
    # We will reuse the construction assembly we created in example 4:
    ConstructionAssembly.db_path = "./example_04_constr_assem.db"
    ca_ext_wall = ConstructionAssembly.load(ID='exterior_wall')
    # Create the exterior building element:
    ext_wall = ExteriorBuildingElement.create(
        ID='south-exterior-wall',
        T_zone=T_zone,
        constr_assem=ca_ext_wall,
        gross_area=Q_(10, 'm') * Q_(3, 'm'),  # surface area of the wall
        weather_data=weather_data,
        gamma=Q_(0, 'deg'),  # oriented to the south
        beta=Q_(90, 'deg')  # vertical wall
    )
    return ext_wall


def add_window(ebe: ExteriorBuildingElement) -> Window:
    """Adds a window to the `ExteriorBuildingElement` object passed to parameter
    `ebe` and returns this window.
    """
    # We will use the window properties of a window type stored on the shelf
    # located in the directory `wtcb-database` of the user's home directory
    # (see __init__.py-file in `cooling_load_calc.wtcb`).
    window_5a = WindowPropertiesShelf.load('window-5a-operable-wood/vinyl')
    ebe.add_window(
        ID='window_1',
        width=Q_(3, 'm'),
        height=Q_(2, 'm'),
        props=window_5a
    )
    return ebe.windows['window_1']


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

    # Create the exterior wall:
    ext_wall = create_exterior_wall(
        weather_data,
        T_zone=lambda t: Q_(22, 'degC')
    )

    # Add the window to the exterior wall:
    window = add_window(ext_wall)

    # CONDUCTION HEAT GAIN AT THE INTERIOR SIDE OF THE EXTERIOR WALL

    # Set the time step as a fraction of an hour between two time moments at
    # which the heat gains will be calculated. By default, this time step is 1
    # hr.
    dt_hr = 1 / 2

    # The conduction heat gain at the interior side of the exterior building
    # element has a convective and a radiative component.
    wall_cnd_Q_dot, wall_cnd_Q_dot_conv, wall_cnd_Q_dot_rad = ext_wall.conductive_heat_gain(
        dt_hr=dt_hr,
        num_cycles=5
    )

    # CONDUCTION HEAT GAIN THROUGH THE WINDOW
    # The conduction heat gain through the window also has a convective and a
    # radiative component.
    wnd_cnd_Q_dot, wnd_cnd_Q_dot_conv, wnd_cnd_Q_dot_rad = window.conductive_heat_gain(dt_hr)

    # SOLAR HEAT GAIN THROUGH THE WINDOW
    wnd_sol_Q_dot, wnd_sol_Q_dot_conv, wnd_sol_Q_dot_rad = window.solar_heat_gain(dt_hr)

    # Display a table with the interior conduction heat gains in the course of
    # time with their corresponding solar and local time at the location:
    n = len(wall_cnd_Q_dot)
    loc_time, sol_time = zip(*[
        convert_to_clock_time(
            time_index=k,
            dt_hr=dt_hr,
            date=weather_data.date,
            L_loc=location.L_loc,
            tz_loc=location.timezone
        ) for k in range(n)
    ])
    d = {
        'solar time': [st.time() for st in sol_time],
        f'{location.timezone}': [lt.time() for lt in loc_time],
        'HG cond. wall': wall_cnd_Q_dot.magnitude,
        'HG cond. window': wnd_cnd_Q_dot.magnitude,
        'HG sol. window': wnd_sol_Q_dot.magnitude
    }
    df = pd.DataFrame(d)
    print('CONDUCTION AND SOLAR HEAT GAIN THROUGH EXTERIOR WALL')
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print(df)

    # Show the heat gains as function of time in a line chart
    chart = LineChart()
    chart.add_xy_data(
        label='HG cond. wall',
        x1_values=df.index,
        y1_values=df['HG cond. wall'],
        style_props={'marker': 'o'}
    )
    chart.add_xy_data(
        label='HG cond. window',
        x1_values=df.index,
        y1_values=df['HG cond. window'],
        style_props={'marker': 'o'}
    )
    chart.add_xy_data(
        label='HG sol. window',
        x1_values=df.index,
        y1_values=df['HG sol. window'],
        style_props={'marker': 'o'}
    )
    chart.add_legend()
    chart.x1.add_title('time index')
    chart.y1.add_title('heat gain, W')
    chart.show()


if __name__ == '__main__':
    main()
