"""Steady-state analysis of a single-stage vapor compression machine at
different compressor speeds using multiprocessing (with `ProcessPoolExecutor`)
and with logging to console and file (note: afterward, the log file can be
sorted by process ID using the method `sort_log_by_process_id` of class
`ModuleLogger`).
"""
import time
import warnings
import traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from hvac import Quantity
from hvac.logging import ModuleLogger, Logger
from hvac.fluids import Fluid, HumidAir, CoolPropWarning
from hvac.vapor_compression import VariableSpeedCompressor
from hvac.vapor_compression.machine_bis import SingleStageVaporCompressionMachine
from hvac.heat_transfer.heat_exchanger.fin_tube import air_evaporator, air_condenser

# Turn off any runtime warnings:
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Turn off warnings coming from module `hvac.fluids.fluid`:
warnings.filterwarnings('ignore', category=CoolPropWarning)

Q_ = Quantity

# Get evaporator and condenser model classes
Evaporator = air_evaporator.rating.PFT_CO_EVP  # plain fin-tube, counter-flow evaporator
Condenser = air_condenser.rating.PFT_CO_CND  # plain fin-tube, counter-flow condenser

# Define the global variables for each process:
machine: SingleStageVaporCompressionMachine = None
logger: Logger = None


def create_machine(logger) -> SingleStageVaporCompressionMachine:
    # Define refrigerant:
    R134a = Fluid('R134a')
    # Create compressor model
    compressor = VariableSpeedCompressor(
        coeff_file=Path("./compressor_data/VTZ038-G_R134a.csv"),
        refrigerant_type=R134a,
        units={'m_dot': 'kg / hr', 'speed': '1 / s'}
    )
    # Create evaporator model:
    evaporator = Evaporator(
        L1=Q_(0.731, 'm'),
        L3=Q_(0.244, 'm'),
        N_r=2,
        S_t=Q_(25.4, 'mm'),
        S_l=Q_(22.0, 'mm'),
        D_i=Q_(8.422, 'mm'),
        D_o=Q_(10.2, 'mm'),
        t_f=Q_(0.3302, 'mm'),
        N_f=1 / Q_(3.175, 'mm')
    )
    # Create condenser model:
    condenser = Condenser(
        L1=Q_(0.909, 'm'),
        L3=Q_(0.303, 'm'),
        N_r=4,
        S_t=Q_(25.4, 'mm'),
        S_l=Q_(22.0, 'mm'),
        D_i=Q_(8.422, 'mm'),
        D_o=Q_(10.2, 'mm'),
        t_f=Q_(0.3302, 'mm'),
        N_f=1 / Q_(3.175, 'mm')
    )
    # Configure single-stage vapor compression machine:
    machine = SingleStageVaporCompressionMachine(
        evaporator,
        condenser,
        compressor,
        R134a,
        n_cmp_min=Q_(2100, '1 / min'),
        n_cmp_max=Q_(5400, '1 / min'),
        logger=logger
    )
    return machine


def analyze_performance(n_cmp: Quantity):
    # This function (the actual task) will run in parallel in different processes.
    # Indicate that in each process the global variables `machine` and `logger`
    # of that process will be used inside this function:
    global machine, logger
    # Set operating conditions on machine:
    # - mass flow rate of air through evaporator
    # - mass flow rate of air through condenser
    # - state of air entering evaporator
    # - state of air entering condenser
    # - degree of superheat set on expansion device
    # - compressor speed
    machine.set_operating_conditions(
        evp_m_dot_air=Q_(1500.0, 'kg / hr'),
        cnd_m_dot_air=Q_(2233.861, 'kg / hr'),
        evp_air_in=HumidAir(Tdb=Q_(24.0, 'degC'), RH=Q_(50, 'pct')),
        cnd_air_in=HumidAir(Tdb=Q_(35.0, 'degC'), RH=Q_(30, 'pct')),
        dT_sh=Q_(10, 'K'),
        n_cmp=n_cmp
    )
    # Run rating routine to fully determine steady-state machine performance
    # It may happen that no solution is possible for the current set of
    # operating conditions. For this compressor speed each of the traced
    # performance parameters is assigned a NaN-value.
    try:
        logger.info(
            f"Rating @ n_cmp = {n_cmp.to('1 / min'):~P.0f}"
        )
        machine.rate_root(
            T_evp_ini=Q_(5.0, 'degC'),  # initial guess evaporator temperature
            T_cnd_ini=Q_(50.0, 'degC')  # initial guess condenser temperature
        )
    except Exception as err:
        # Add exception traceback to the log message:
        traceback_str = traceback.format_exc()
        logger.error(f"{err}, {traceback_str}")
        return {
            'n_cmp [rpm]': n_cmp.to('1 / min').m,
            'T_evp [degC]': Q_(float('nan'), 'degC').m,
            'P_evp [bar]': Q_(float('nan'), 'bar').m,
            'T_cnd [degC]': Q_(float('nan'), 'degC').m,
            'P_cnd [bar]': Q_(float('nan'), 'bar').m,
            'Q_evp [kW]': Q_(float('nan'), 'kW').m,
            'W_cmp [kW]': Q_(float('nan'), 'kW').m,
            'Q_cnd [kW]': Q_(float('nan'), 'kW').m,
            'COP [-]': Q_(float('nan'), 'frac').m,
            'm_dot_rfg [kg/h]': Q_(float('nan'), 'kg / hr').m,
            'dT_sc [K]': Q_(float('nan'), 'K').m,
            'eps_evp [-]': Q_(float('nan'), 'frac').m,
            'eps_cnd [-]': Q_(float('nan'), 'frac').m,
            'TDB evp_air_out [degC]': Q_(float('nan'), 'degC').m,
            'RH evp_air_out [%]': Q_(float('nan'), 'pct').m,
            'dP evp_air_out [Pa]': Q_(float('nan'), 'Pa').m,
            'TDB cnd_air_out [degC]': Q_(float('nan'), 'degC').m,
            'RH cnd_air_out [%]': Q_(float('nan'), 'pct').m,
            'dP cnd_air_out [Pa]': Q_(float('nan'), 'Pa').m,
            'mass balance error [%]': Q_(float('nan'), 'pct').m,
            'energy balance error [%]': Q_(float('nan'), 'pct').m
        }
    else:
        return {
            'n_cmp [rpm]': n_cmp.to('1 / min').m,
            'T_evp [degC]': machine.Te.to('degC').m,
            'P_evp [bar]': machine.Pe.to('bar').m,
            'T_cnd [degC]': machine.Tc.to('degC').m,
            'P_cnd [bar]': machine.Pc.to('bar').m,
            'Q_evp [kW]': machine.Qc_dot.to('kW').m,
            'W_cmp [kW]': machine.Wc_dot.to('kW').m,
            'Q_cnd [kW]': machine.Qh_dot.to('kW').m,
            'COP [-]': machine.COP.to('frac').m,
            'm_dot_rfg [kg/h]': machine.m_dot.to('kg / hr').m,
            'dT_sc [K]': machine.sub_cooling.to('K').m,
            'eps_evp [-]': machine.evaporator.eps.to('frac').m,
            'eps_cnd [-]': machine.condenser.eps.to('frac').m,
            'TDB evp_air_out [degC]': machine.evaporator.air_out.Tdb.to('degC').m,
            'RH evp_air_out [%]': machine.evaporator.air_out.RH.to('pct').m,
            'dP evp_air_out [Pa]': machine.evaporator.dP_air.to('Pa').m,
            'TDB cnd_air_out [degC]': machine.condenser.air_out.Tdb.to('degC').m,
            'RH cnd_air_out [%]': machine.condenser.air_out.RH.to('pct').m,
            'dP cnd_air_out [Pa]': machine.condenser.dP_air.to('Pa').m,
            'mass balance error [%]': machine.check_mass_balance()[1].m,
            'energy balance error [%]': machine.check_energy_balance()[1].m
        }


def _init_worker(log_file_path: Path):
    global logger, machine
    # Create logger:
    logger = ModuleLogger.get_logger(
        __name__,
        file_path=log_file_path
    )
    # Note: The logger of each process writes directly to the same file.
    # Create machine:
    machine = create_machine(logger)


def main():
    # Define path to log file:
    log_file_path = Path(
        "C:/Users/Tom/PycharmProjects/python-hvac/"
        "docs/examples/vapor_compression/log_files/vcm_ex07_1.log"
    )
    # Delete previous log file if it exists:
    if log_file_path.exists():
        log_file_path.unlink()
    # Run `analyze_performance` with different compressor speeds parallel in
    # multiple processes with `ProcessPoolExecutor`. Before a process is
    # started, it is first initialized by calling function `_init_worker`,
    # passing the `log_file_path` to it to write any log messages to this file:
    n_cmp_rng = [Q_(n, '1 / min') for n in range(2100, 5700, 300)]
    data = []
    with ProcessPoolExecutor(initializer=_init_worker, initargs=(log_file_path,)) as executor:
        start_time = time.perf_counter()
        for record in executor.map(analyze_performance, n_cmp_rng):
            data.append(record)
        finish_time = time.perf_counter()
        print(f"Run time = {finish_time - start_time} seconds")
    # Put the records in a Pandas DataFrame object:
    df = pd.DataFrame(data)
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print(df)
    # Save dataframe to spreadsheet file:
    df.to_excel('./analysis_results/07_multi_analysis_1.ods')


if __name__ == '__main__':
    main()
    # Function `main` is only called by the main process.
