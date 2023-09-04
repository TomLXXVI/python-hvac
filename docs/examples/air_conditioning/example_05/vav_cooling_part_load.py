"""
Steady-state analysis of a single-zone VAV air cooling system hour by hour
during the summer peak design day.
"""
import warnings
import pathlib
from hvac.fluids import CoolPropWarning
warnings.filterwarnings('ignore', category=CoolPropWarning)

from datetime import date
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid
from hvac.climate import ClimateData, Location
from hvac.heat_transfer.heat_exchanger.fin_tube.core.plain_fin_tube import PlainFinTubeHeatExchangerCore
from hvac.air_conditioning.single_zone.vav_cooling_sim import (
    DesignData, VAVSingleZoneAirCoolingSystem,
    CoolingSimData, Output
)
from hvac.charts import LineChart, PsychrometricChart, StatePoint
from hvac.logging import ModuleLogger, Logger


Q_ = Quantity

logger: Logger


# Summer peak design data:
design_data = DesignData(
    Q_zone_sen=Q_(36.592, 'kW'),
    Q_zone_lat=Q_(11.0, 'kW'),
    V_dot_vent_ntp=Q_(4160.611, 'm ** 3 / hr'),
    V_dot_supply_ntp=Q_(10_721.476, 'm ** 3 / hr'),
    outdoor_air=HumidAir(Tdb=Q_(32, 'degC'), Twb=Q_(21, 'degC')),
    zone_air=HumidAir(Tdb=Q_(26, 'degC'), RH=Q_(50, 'pct')),
    supply_air=HumidAir(Tdb=Q_(14.0, 'degC'), W=Q_(9.1, 'g / kg'))
)


# noinspection PyUnusedLocal
def _warn_handler(message, category, filename, lineno, file=None, line=None):
    global logger
    logger.error(message)


def task(args) -> Output:
    """Performs the analysis of the AC system at the outdoor air state,
     the relative sensible cooling load, and the relative latent cooling load
     contained in `args`.

     Notes
     -----
     Function `task` will run in separate child processes. Therefore, each call
     to `task` should create its own `VAVSingleZoneAirCoolingSystem` object.
     """
    global design_data, logger

    outdoor_air = args[0]
    rel_Q_zone_sen = args[1]
    rel_Q_zone_lat = args[2]

    logger.info(
        "Analyze AC system with "
        f"outdoor_air = ({outdoor_air.Tdb.to('degC'):~P.2f} DB, "
        f"{outdoor_air.RH.to('pct'):~P.2f} RH), "
        f"sensible cooling load fraction = {rel_Q_zone_sen.to('pct'):~P.2f}, "
        f"latent cooling load fraction = {rel_Q_zone_lat.to('pct'):~P.2f}"
    )

    # Define the heat exchanger core of the DX air-cooling coil:
    hex_core = PlainFinTubeHeatExchangerCore(
        L1=Q_(1297, 'mm'),
        L3=Q_(1000, 'mm'),
        S_t=Q_(22.42, 'mm'),
        S_l=Q_(22.27, 'mm'),
        D_i=Q_(8.422, 'mm'),
        D_o=Q_(10.2, 'mm'),
        t_f=Q_(0.3302, 'mm'),
        N_f=1 / Q_(2.8, 'mm'),
        k_f=Q_(237, 'W / (m * K)'),
        N_r=3
    )
    # Configure the AC system:
    ac_system = VAVSingleZoneAirCoolingSystem(
        design_data=design_data,
        dx_coil_hex_core=hex_core,
        dx_coil_Rfg=Fluid('R410a'),
        dx_coil_T_evp=Q_(9, 'degC'),
        T_zone=Q_(26, 'degC'),    # setpoint zone air temperature
        T_cool=Q_(14, 'degC'),    # setpoint cooled air temperature
        heating_coil_present=False
    )
    with warnings.catch_warnings(category=RuntimeWarning):
        warnings.showwarning = _warn_handler
        # Analyze part-load operation:
        output = ac_system.analyze(
            outdoor_air=outdoor_air,
            rel_Q_zone_sen=rel_Q_zone_sen,
            rel_Q_zone_lat=rel_Q_zone_lat
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

    # Get hourly dry-bulb and wet-bulb temperature profile on summer peak design
    # day:
    climate = ClimateData.create(
        location=Location(
            name='Ghent',
            lat=Q_(51.183, 'deg'),
            lon=Q_(3.8, 'deg'),
            alt=Q_(8.0, 'm'),
            tz='Europe/Brussels'
        ),
        design_day=date(2022, 7, 21),
        Tdb_avg=Q_(31.6, 'degC'),
        Tdb_range=Q_(7.3, 'K'),
        Twb_mc=Q_(20.6, 'degC'),
        tau_beam=0.426,
        tau_dif=2.247
    )
    T_outdoor_db_rng = climate.Tdb_profile['T'][:-1]
    T_outdoor_wb_rng = climate.Twb_profile['T'][:-1]
    outdoor_air_rng = [
        HumidAir(Tdb=Tdb, Twb=Twb)
        for Tdb, Twb in zip(T_outdoor_db_rng, T_outdoor_wb_rng)
    ]

    # Prepare simulation data:
    sim_data = CoolingSimData(
        T_outdoor_db_rng=[oa.Tdb for oa in outdoor_air_rng],
        T_outdoor_wb_rng=[oa.Twb for oa in outdoor_air_rng],
        df_Q_zone=pd.read_excel("cooling_load.ods"),
        design_data=design_data
    )

    # Plot the daily profile of hourly dry-bulb and wet-bulb outdoor air
    # temperature:
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
        for output in executor.map(task, sim_data()):
            outputs.append(output)
    
    # Show results:
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
    chart_01.y1.scale(0, 160, 10)
    chart_01.add_legend()
    chart_01.show()
    
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
        style_props={'color': 'tab:blue'}
    )
    chart_02.add_xy_data(
        label='return air temperature',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.return_air.Tdb.to('degC').m for output in outputs],
        style_props={'color': 'tab:orange'}
    )
    chart_02.add_xy_data(
        label='dry-bulb outdoor air temperature',
        x1_values=[h for h in range(24)],
        y1_values=[oa.Tdb.to('degC').m for oa in outdoor_air_rng],
        style_props={'color': 'tab:red'}
    )
    chart_02.x1.add_title('hour of the day')
    chart_02.x1.scale(0, 24, 1)
    chart_02.y1.add_title('air temperature, °C')
    chart_02.y1.scale(10, 45, 5)
    chart_02.add_legend()
    chart_02.show()

    chart_03 = LineChart()
    chart_03.add_xy_data(
        label='supply air humidity',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.supply_air.W.to('g / kg').m for output in outputs],
        style_props={'color': 'tab:cyan'}
    )
    chart_03.add_xy_data(
        label='cooling air humidity',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.cooled_air.W.to('g / kg').m for output in outputs],
        style_props={'color': 'tab:blue'}
    )
    chart_03.add_xy_data(
        label='return air humidity',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.return_air.W.to('g / kg').m for output in outputs],
        style_props={'color': 'tab:orange'}
    )
    chart_03.add_xy_data(
        label='outdoor air humidity',
        x1_values=[i for i in range(len(outdoor_air_rng))],
        y1_values=[oa.W.to('g / kg').m for oa in outdoor_air_rng],
        style_props={'color': 'tab:red'}
    )
    chart_03.x1.add_title('hour of the day')
    chart_03.x1.scale(0, 24, 1)
    chart_03.y1.add_title('air humidity, g/kg')
    chart_03.add_legend()
    chart_03.show()

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
        style_props={'color': 'tab:red'}
    )
    chart_04.add_xy_data(
        label='recirculated air NTP volume flow rate',
        x1_values=[i for i in range(len(outputs))],
        y1_values=[output.V_dot_recir_ntp.to('m ** 3 / hr').m for output in outputs],
        style_props={'color': 'tab:orange'}
    )
    chart_04.x1.add_title('hour of the day')
    chart_04.x1.scale(0, 24, 1)
    chart_04.y1.add_title('air volume flow rate, m³/h')
    chart_04.y1.scale(0, 11000, 1000)
    chart_04.add_legend()
    chart_04.show()
    
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

    chart_06 = PsychrometricChart()
    state_points = [
        StatePoint(output.return_air.Tdb, output.return_air.W)
        for output in outputs
    ]
    for i, state_point in enumerate(state_points):
        chart_06.plot_point(f"Pnt. {i}", state_point)
    chart_06.show()

    ModuleLogger.sort_log_by_process_id(file_path)


if __name__ == '__main__':
    main()
