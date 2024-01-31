"""
Steady-state analysis on an hourly basis of a single-zone VAV air heating and
humidification system during a selected winter day.
"""
import warnings
import pathlib

from hvac.fluids import CoolPropWarning
warnings.filterwarnings('ignore', category=CoolPropWarning)

from concurrent.futures import ProcessPoolExecutor

from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.sun.tmy import TMY
from hvac.air_conditioning.single_zone.vav_heating_sim import (
    DesignData, VAVSingleZoneAirHeatingSystem, Output
)
from hvac.charts import LineChart, PsychrometricChart, StatePoint
from hvac.logging import ModuleLogger, Logger


Q_ = Quantity

logger: Logger

# Winter peak design data:
design_data = DesignData(
    Q_dot_zone_sen=Q_(-41.487, 'kW'),
    Q_dot_zone_lat=Q_(-8.792, 'kW'),  # the zone loses moisture @ design
    V_dot_vent_ntp=Q_(1699.011, 'm ** 3 / hr'),
    V_dot_supply_ntp=Q_(7612.546, 'm ** 3 / hr'),
    outdoor_air=HumidAir(Tdb=Q_(0.51, 'degC'), RH=Q_(53, 'pct')),
    zone_air=HumidAir(Tdb=Q_(22.22, 'degC'), RH=Q_(50, 'pct')),
    supply_air=HumidAir(Tdb=Q_(38.188, 'degC'), RH=Q_(24, 'pct')),
    T_steam=Q_(93.33, 'degC'),
    Q_dot_zone_ihg=Q_(0, 'kW')
)


def task(outdoor_air: HumidAir) -> Output:
    """Performs the analysis of the AC system at the given outdoor air state.

     Notes
     -----
     Function `task` will run in separate child processes. Therefore, each call
     to `task` should create its own `VAVSingleZoneAirHeatingSystem` object.
     """
    global design_data, logger

    logger.info(
        "Analyze AC system with "
        f"outdoor_air = ({outdoor_air.Tdb.to('degC'):~P.2f} DB, "
        f"{outdoor_air.RH.to('pct'):~P.2f} RH), "
    )

    ac_system = VAVSingleZoneAirHeatingSystem(
        design_data=design_data,
        humidifier_present=True
    )

    output = ac_system.analyze(
        outdoor_air=outdoor_air,
        rel_Q_zone_sen=None,  # sensible load depends on outdoor air temperature and K_zone
        rel_Q_zone_lat=-Q_(10, 'pct')  # moisture is added to the zone
    )

    logger.info(
        "Zone air state: "
        f"{output.return_air.Tdb.to('degC'):~P.2f} DB, "
        f"{output.return_air.W.to('g / kg'):~P.2f} AH "
        f"({output.return_air.RH.to('pct'):~P.2f} RH)"
    )

    return output


def init_worker(file_path: pathlib.Path):
    global logger
    logger = ModuleLogger.get_logger(__name__, file_path=str(file_path))


def main():
    global design_data

    file_path = pathlib.Path('sim.log')
    if file_path.exists():
        file_path.unlink()

    # Get the simulation data for the selected day of the year (hourly dry-bulb
    # temperature and relative humidity) from a TMY-datafile:
    tmy = TMY(file_path='tmy_gent_2005_2020.csv')
    selected_day = tmy.get_day(month=12, day=31)

    T_outdoor_rng = Q_(selected_day['T2m'].values, 'degC')
    RH_outdoor_rng = Q_(selected_day['RH'].values, 'pct')
    outdoor_air_rng = [
        HumidAir(Tdb=Tdb, RH=RH)
        for Tdb, RH in zip(T_outdoor_rng, RH_outdoor_rng)
    ]

    # Plot the daily profile of hourly dry-bulb and wet-bulb outdoor air
    # temperature on the selected day:
    chart_00 = LineChart()
    chart_00.add_xy_data(
        label='outdoor air dry-bulb temperature',
        x1_values=[h for h in range(24)],
        y1_values=[oa.Tdb.to('degC').m for oa in outdoor_air_rng]
    )
    chart_00.add_xy_data(
        label='outdoor air wet-bulb temperature',
        x1_values=[h for h in range(24)],
        y1_values=[oa.Twb.to('degC').m for oa in outdoor_air_rng]
    )
    chart_00.x1.add_title('hour of the day')
    chart_00.x1.scale(0, 24, 1)
    chart_00.y1.add_title('temperature, °C')
    chart_00.add_legend()
    chart_00.show()

    # Run the simulation:
    outputs = []
    with ProcessPoolExecutor(initializer=init_worker, initargs=(file_path,)) as executor:
        for output in executor.map(task, outdoor_air_rng):
            outputs.append(output)

    # Show the results:
    chart_01 = LineChart()
    chart_01.add_xy_data(
        label='heating coil load',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.Q_dot_hc.to('kW').m for output in outputs]
    )
    chart_01.x1.add_title('hour of the day')
    chart_01.x1.scale(0, 24, 1)
    chart_01.y1.add_title('heating coil load, kW')
    chart_01.y1.scale(0, 80, 10)
    chart_01.show()

    chart_02 = LineChart()
    chart_02.add_y2_axis()
    chart_02.add_xy_data(
        label='supply air NTP volume flow rate',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.V_dot_supply_ntp.to('m ** 3 / hr').m for output in outputs]
    )
    chart_02.add_xy_data(
        label='ventilation air NTP volume flow rate',
        x1_values=[i for i in range(len(outputs))],
        y2_values=[output.V_dot_vent_ntp.to('m ** 3 / hr').m for output in outputs],
        style_props={'color': 'orange'}
    )
    chart_02.x1.add_title('hour of the day')
    chart_02.x1.scale(0, 24, 1)
    chart_02.y1.add_title('supply air volume flow rate, m³/h')
    chart_02.y2.add_title('ventilation air volume flow rate, m³/h')
    chart_02.y1.scale(1000, 10000, 1000)
    chart_02.y2.scale(1000, 10000, 1000)
    chart_02.add_legend()
    chart_02.show()

    chart_03 = LineChart()
    chart_03.add_xy_data(
        label='supply air temperature',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.supply_air.Tdb.to('degC').m for output in outputs],
    )
    chart_03.x1.add_title('hour of the day')
    chart_03.x1.scale(0, 24, 1)
    chart_03.y1.add_title('supply air temperature, °C')
    chart_03.y1.scale(0, 55, 5)
    chart_03.show()

    chart_04 = LineChart()
    chart_04.add_y2_axis()
    chart_04.add_xy_data(
        label='absolute humidity of zone air',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.return_air.W.to('g / kg').m for output in outputs],
    )
    chart_04.add_xy_data(
        label='relative humidity of zone air',
        x1_values=[i for i in range(len(outputs))],
        y2_values=[output.return_air.RH.to('pct').m for output in outputs],
        style_props={'color': 'orange'}
    )
    chart_04.x1.add_title('hour of the day')
    chart_04.x1.scale(0, 24, 1)
    chart_04.y1.add_title('absolute humidity of zone air, g/kg')
    chart_04.y2.add_title('relative humidity of zone air, %')
    chart_04.y1.scale(0, 16, 1)
    chart_04.y2.scale(0, 110, 10)
    chart_04.add_legend()
    chart_04.show()

    chart_05 = LineChart()
    chart_05.add_xy_data(
        label='zone air temperature',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.return_air.Tdb.to('degC').m for output in outputs]
    )
    chart_05.x1.add_title('hour of the day')
    chart_05.x1.scale(0, 24, 1)
    chart_05.y1.add_title('zone air temperature, °C')
    chart_05.y1.scale(0, 32, 2)
    chart_05.show()

    chart_06 = LineChart()
    chart_06.add_xy_data(
        label='humidification',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.m_dot_steam.to('kg / hr').m for output in outputs]
    )
    chart_06.x1.add_title('hour of the day')
    chart_06.x1.scale(0, 24, 1)
    chart_06.y1.add_title('steam mass flow rate, kg/h')
    chart_06.y1.scale(0, 32, 2)
    chart_06.show()

    chart_07 = PsychrometricChart()
    state_points = [
        StatePoint(output.return_air.Tdb, output.return_air.W)
        for output in outputs
    ]
    for i, state_point in enumerate(state_points):
        chart_07.plot_point(f"Pnt. {i}", state_point)
    chart_07.show()

    ModuleLogger.sort_log_by_process_id(file_path)


if __name__ == '__main__':
    main()
