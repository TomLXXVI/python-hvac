"""
PART-LOAD OPERATION ANALYSIS OF A SINGLE-ZONE VAV AIR-COOLING SYSTEM WITH
ECONOMIZER.
"""
from datetime import date
import dill as pickle
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.charts import LineChart, PsychrometricChart, StatePoint
from hvac.air_conditioning.single_zone.vav_cooling_sim import (
    VAVSingleZoneAirCoolingSystem,
    DesignData,
    CoolingSimData,
    Output
)
from hvac.climate import (
    ClimateData,
    Location,
    TMY
)
from exposition_hall import ExpositionHall


Q_ = Quantity

# ------------------------------------------------------------------------------
# CONFIGURATION OF THE SINGLE-ZONE VAV AIR-COOLING SYSTEM

# To create the `VAVSingleZoneAirCoolingSystem`, we need the design load
# calculations for the summer peak design day:
# - the sensible cooling zone load
# - the latent cooling zone load
# - the minimum required volume flow rate of outdoor ventilation air
# - the volume flow rate of supply air
# - the state of outdoor air under summer peak design conditions
# - the desired state of zone air
# - the state of supply air

# To indicate that the air-cooling system also has a heating coil, we can set
# parameter `heating_coil_present` to True. This is the default value.
# It is also possible to omit the heating coil (False) and instead set parameter
# `variable_cooling_setpoint` to True. In that case, the setpoint temperature of
# the air leaving the cooling coil is the same as the required supply air
# temperature needed to keep the zone air temperature constant.
# If both parameters are set to False, the cooling coil operates with a constant
# setpoint temperature for the air leaving the cooling coil and being supplied
# to the zone.

# Parameter `rel_m_dot_supply_min` represents the minimum allowed supply air mass
# flow rate as a fraction of the design supply air mass flow rate. This minimum
# flow rate must still ensure the proper mixing of supply air with zone air.
# By default, it is set to 60 % of the design value.

# Finally, the optional parameter `units` can be used to set the measuring units
# of quantities when displaying the results.

airco_system = VAVSingleZoneAirCoolingSystem(
    design_data=DesignData(
        Q_zone_sen=Q_(31.126, 'kW'),
        Q_zone_lat=Q_(17.936, 'kW'),
        V_dot_vent_ntp=Q_(4400, 'm ** 3 / hr'),
        V_dot_supply_ntp=Q_(7599.946, 'm ** 3 / hr'),
        outdoor_air=HumidAir(Tdb=Q_(26.7, 'degC'), Twb=Q_(19.2, 'degC')),
        zone_air=HumidAir(Tdb=Q_(26.0, 'degC'), RH=Q_(50, 'pct')),
        supply_air=HumidAir(Tdb=Q_(14, 'degC'), W=Q_(7.783, 'g / kg')),
    ),
    heating_coil_present=True
)

# ------------------------------------------------------------------------------
# WEATHER DATA AND PART-LOAD ZONE LOADS

# We use hourly TMY-data for a certain day of the year to get outdoor air
# temperatures and to determine the part-load zone loads for each hour of the
# day. See class `TMY` in module `tmy` of subpackage `climate` for more
# information.

# We also need a thermal model of our single-zone building to calculate the zone
# loads at part-load. For this, we use the model `ExpositionHall` we created in
# example 01 of the cooling load calculation examples (see docs/examples/
# cooling_load_calc).

# The function `get_simulation_data` returns for a specified day (indicated by
# a month and day index) a `CoolingSimData` object. See module `vav_cooling_sim`
# in subpackage `air_conditioning.single_zone` for more information. We use this
# object to get the outdoor air temperature and the relative sensible and latent
# zone loads at part-load for each hour of the day.


def get_simulation_data(
    month: int, day: int,
    design_data: DesignData
) -> CoolingSimData:
    climate = ClimateData.create_from_TMY_data(
        day=date(2022, month, day),
        location=Location(
            name='Ghent',
            lat=Q_(51.183, 'deg'),
            lon=Q_(3.8, 'deg'),
            alt=Q_(8.0, 'm'),
            tz='Europe/Brussels'
        ),
        tmy=TMY('tmy_gent_2005_2020.csv')
    )

    hall = ExpositionHall.create(
        climate,
        max_num_people=50,
        min_num_people=5
    )

    sim_data = CoolingSimData(
        T_outdoor_db_rng=climate.Tdb_profile['T'],
        T_outdoor_wb_rng=climate.Twb_profile['T'],
        df_Q_zone=hall.get_heat_gains(unit='kW'),
        design_data=design_data
    )
    return sim_data


# ------------------------------------------------------------------------------
# The function `save_data` is just a helper function to pickle some results
# of the simulation for later use in other scripts. It is not really necessary
# to use this function in a simulation.

def save_data(outputs: list[Output], file_path: str) -> None:
    data = [
        (
            output.mixed_air.Tdb.to('degC'),
            output.mixed_air.W.to('g / kg'),
            output.m_dot_supply.to('kg / hr')
        )
        for output in outputs
    ]
    with open(file_path, 'wb') as fh:
        pickle.dump(data, fh)


# ------------------------------------------------------------------------------
# SIMULATION OF PART-LOAD OPERATION AT A GIVEN DAY

def main():
    # We will save the line charts generated further on this function to a
    # subfolder in the current working directory.
    # The results we will pickle with the function `save_data` will be saved in
    # another subfolder (this can be omitted).
    chart_folder = "./charts/"
    data_folder = "./data/"

    # Get the simulation working data for the given month and day:
    sim_data = get_simulation_data(7, 21, airco_system.design_data)

    # Run the simulation for each hour of the day we have selected:
    outputs = []
    for hour, (outdoor_air, rel_Q_zone_sen, rel_Q_zone_lat) in enumerate(sim_data()):
        print(f"hour: {hour}")
        print(
            "outdoor air: "
            f"{outdoor_air.Tdb.to('degC'):~P.3f} DB, "
            f"{outdoor_air.W.to('g / kg'):~P.3f} AH "
            f"({outdoor_air.RH.to('pct'):~P.3f} RH)",
            "relative sensible zone load: "
            f"{rel_Q_zone_sen.to('pct'):~P.3f}",
            "relative latent zone load: "
            f"{rel_Q_zone_lat.to('pct'):~P.3f}",
            sep='\n'
        )

        # Analyze part-load operation of the air-cooling system by passing the
        # state of outdoor air and the relative sensible and latent zone loads,
        # and the setpoint of the zone air temperature at the given hour of the
        # selected day:
        output = airco_system.analyze(
            outdoor_air=outdoor_air,
            rel_Q_zone_sen=rel_Q_zone_sen,
            rel_Q_zone_lat=rel_Q_zone_lat,
            T_set_zone=Q_(26, 'degC')
        )
        outputs.append(output)
        print(output)

    # Pickle the results in `outputs` with the function `save_data` (this can
    # be omitted).
    save_data(outputs, data_folder + "data.pickle")

    # DRAW LINE CHARTS WITH THE RESULTS

    # Daily profile of dry-bulb and wet-bulb outdoor air temperature:
    chart_00 = LineChart()
    chart_00.add_xy_data(
        label='outdoor air dry-bulb temperature',
        x1_values=[h for h in range(24)],
        y1_values=[Tdb.to('degC').m for Tdb in sim_data.T_outdoor_db_rng]
    )
    chart_00.add_xy_data(
        label='outdoor air wet-bulb temperature',
        x1_values=[h for h in range(24)],
        y1_values=[Twb.to('degC').m for Twb in sim_data.T_outdoor_wb_rng]
    )
    chart_00.x1.add_title('hour of the day')
    chart_00.x1.scale(0, 24, 1)
    chart_00.y1.add_title('temperature, °C')
    chart_00.add_legend()
    chart_00.show()
    chart_00.save('chart_00', location=chart_folder)

    # Relative zone loads:
    chart_01 = LineChart()
    chart_01.add_xy_data(
        label='relative sensible zone load',
        x1_values=[i for i in range(len(sim_data.rel_Q_zone_sen_rng))],
        y1_values=[
            rel_Q_zone_sen.to('pct').m
            for rel_Q_zone_sen in sim_data.rel_Q_zone_sen_rng
        ]
    )
    chart_01.add_xy_data(
        label='relative latent zone load',
        x1_values=[i for i in range(len(sim_data.rel_Q_zone_lat_rng))],
        y1_values=[
            rel_Q_zone_lat.to('pct').m
            for rel_Q_zone_lat in sim_data.rel_Q_zone_lat_rng
        ]
    )
    chart_01.x1.add_title('hour of the day')
    chart_01.x1.scale(0, 24, 1)
    chart_01.y1.add_title('relative zone load, %')
    chart_01.y1.scale(-20, 220, 20)
    chart_01.add_legend()
    chart_01.show()
    chart_01.save('chart_01', location=chart_folder)

    # Temperatures:
    chart_02 = LineChart()
    chart_02.add_xy_data(
        label='supply air temperature',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.supply_air.Tdb.to('degC').m for output in outputs],
        style_props={'color': 'tab:cyan'}
    )
    chart_02.add_xy_data(
        label='cooling air temperature',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.cooled_air.Tdb.to('degC').m for output in outputs],
        style_props={'color': 'tab:blue', 'linestyle': ':'}
    )
    chart_02.add_xy_data(
        label='mixed air temperature',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.mixed_air.Tdb.to('degC').m for output in outputs],
        style_props={'color': 'tab:red', 'linestyle': '--'}
    )
    chart_02.add_xy_data(
        label='return air temperature',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.return_air.Tdb.to('degC').m for output in outputs],
        style_props={'color': 'tab:orange', 'linestyle': ':'}
    )
    chart_02.add_xy_data(
        label='outdoor air temperature',
        x1_values=[h for h in range(24)],
        y1_values=[Tdb.to('degC').m for Tdb in sim_data.T_outdoor_db_rng],
        style_props={'color': 'tab:red', 'linestyle': ':'}
    )
    chart_02.x1.add_title('hour of the day')
    chart_02.x1.scale(0, 24, 1)
    chart_02.y1.add_title('air temperature, °C')
    chart_02.y1.scale(10, 45, 5)
    chart_02.add_legend()
    chart_02.show()
    chart_02.save('chart_02', location=chart_folder)

    # Humidity levels:
    chart_03 = LineChart()
    chart_03.add_xy_data(
        label='supply air humidity',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.supply_air.W.to('g / kg').m for output in outputs],
        style_props={'color': 'tab:cyan', 'linestyle': ':'}
    )
    chart_03.add_xy_data(
        label='cooling air humidity',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.cooled_air.W.to('g / kg').m for output in outputs],
        style_props={'color': 'tab:blue', 'linestyle': ':'}
    )
    chart_03.add_xy_data(
        label='mixed air humidity',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.mixed_air.W.to('g / kg').m for output in outputs],
        style_props={'color': 'tab:red', 'linestyle': '--'}
    )
    chart_03.add_xy_data(
        label='return air humidity',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.return_air.W.to('g / kg').m for output in outputs],
        style_props={'color': 'tab:orange', 'linestyle': ':'}
    )
    chart_03.add_xy_data(
        label='outdoor air humidity',
        x1_values=[i for i in range(len(sim_data.outdoor_air_rng))],
        y1_values=[oa.W.to('g / kg').m for oa in sim_data.outdoor_air_rng],
        style_props={'color': 'tab:red', 'linestyle': ':'}
    )
    chart_03.x1.add_title('hour of the day')
    chart_03.x1.scale(0, 24, 1)
    chart_03.y1.add_title('air humidity, g/kg')
    chart_03.add_legend()
    chart_03.show()
    chart_03.save('chart_03', location=chart_folder)

    # Volume flow rates:
    chart_04 = LineChart()
    chart_04.add_xy_data(
        label='supply air NTP volume flow rate',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.V_dot_supply_ntp.to('m ** 3 / hr').m for output in outputs],
        style_props={'color': 'tab:cyan'}
    )
    chart_04.add_xy_data(
        label='ventilation air NTP volume flow rate',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.V_dot_vent_ntp.to('m ** 3 / hr').m for output in outputs],
        style_props={'color': 'tab:red', 'linestyle': ':'}
    )
    chart_04.add_xy_data(
        label='recirculated air NTP volume flow rate',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.V_dot_recir_ntp.to('m ** 3 / hr').m for output in outputs],
        style_props={'color': 'tab:orange', 'linestyle': ':'}
    )
    chart_04.x1.add_title('hour of the day')
    chart_04.x1.scale(0, 24, 1)
    chart_04.y1.add_title('air volume flow rate, m³/h')
    chart_04.y1.scale(0, 11000, 1000)
    chart_04.add_legend()
    chart_04.show()
    chart_04.save('chart_04', location=chart_folder)

    # Cooling coil and heating coil loads:
    chart_05 = LineChart()
    chart_05.add_xy_data(
        label='cooling coil load',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.Q_dot_cc.to('kW').m for output in outputs],
        style_props={'color': 'tab:blue'}
    )
    chart_05.add_xy_data(
        label='heating coil load',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.Q_dot_hc.to('kW').m for output in outputs],
        style_props={'color': 'tab:red'}
    )
    chart_05.x1.add_title('hour of the day')
    chart_05.x1.scale(0, 24, 1)
    chart_05.y1.add_title('cooling/heating coil load, kW')
    chart_05.y1.scale(0, 80, 10)
    chart_05.add_legend()
    chart_05.show()
    chart_05.save('chart_05', location=chart_folder)

    # States of zone air:
    chart_06 = PsychrometricChart()
    state_points = [
        StatePoint(output.return_air.Tdb, output.return_air.W)
        for output in outputs
    ]
    for i, state_point in enumerate(state_points):
        chart_06.plot_point(f"Pnt. {i}", state_point)
    chart_06.show()


if __name__ == '__main__':
    main()
