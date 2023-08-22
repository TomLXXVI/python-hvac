"""Analyzing a single-stage vapor compression machine at different compressor
speeds and at different outdoor temperatures. Other operating conditions
(mass flow rate of air at condenser, entering air temperature and mass flow
rate at evaporator, set degree of superheat) remain fixed.

Part 1: Create the model, run the simulation, and save the simulation results
to disk.
"""
import time
import pathlib
import warnings
import traceback
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import Fluid, HumidAir, CoolPropWarning
from hvac.heat_transfer.heat_exchanger.fin_tube.air_condenser.rating import PlainFinTubeCounterFlowCondenser
from hvac.heat_transfer.heat_exchanger.fin_tube.air_evaporator.rating import PlainFinTubeCounterFlowEvaporator
from hvac.vapor_compression.real_compressor import VariableSpeedCompressor
from hvac.vapor_compression.machine_bis import SingleStageVaporCompressionMachine


Q_ = Quantity

# Hide warnings from module hvac.fluids.fluid and RuntimeWarnings:
warnings.filterwarnings('ignore', category=CoolPropWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)


# Get evaporator and condenser model classes:
Evaporator = PlainFinTubeCounterFlowEvaporator
Condenser = PlainFinTubeCounterFlowCondenser


def create_inputs() -> list[tuple[Quantity, Quantity]]:
    """Creates list of outdoor temperatures and compressor speeds at which the
    performance of the vapor compression machine is to be analyzed.
    """
    T_outdoor_rng = Q_(np.arange(26.0, 40.0, 2.0), 'degC')
    n_cmp_rng = Q_(np.arange(2100, 5700, 300), '1 / min')
    input_list = [
        (T_outdoor, n_cmp)
        for n_cmp in n_cmp_rng
        for T_outdoor in T_outdoor_rng
    ]
    return input_list


def create_machine(logger) -> SingleStageVaporCompressionMachine:
    """Creates the components of the single-stage vapor compression machine
    and configures the machine.
    """
    # Define refrigerant:
    R134a = Fluid('R134a')
    # Create compressor model
    compressor = VariableSpeedCompressor(
        coeff_file=pathlib.Path('./compressor_data/VTZ038-G_R134a.csv'),
        refrigerant_type=R134a,
        units={'m_dot': 'kg / hr', 'speed': '1 / s'}
    )
    # Create evaporator model
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
    # Create condenser model
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
    # Configure single-stage vapor compression machine model
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


def analyze_performance(inputs: tuple[Quantity, Quantity]) -> dict[str, float]:
    """Analyzes performance of vapor compression machine at given inputs.
    Parameter `inputs` is a tuple of outdoor temperature and compressor speed.
    Other input parameters that determine machine performance remain fixed:
    - mass flow rate of air through evaporator
    - mass flow rate of air through condenser
    - state of air entering the evaporator
    - degree of superheat set on expansion device
    - relative humidity of outdoor air entering the condenser
    Returns a dictionary with the calculated performance parameters:
    - COP
    - mass flow rate of refrigerant in order to maintain set degree of superheat
    - evaporator pressure
    - condenser pressure
    - heat transfer effectiveness of evaporator
    - heat transfer effectiveness of condenser
    - heat absorption rate at evaporator
    - heat rejection rate at condenser
    - compressor power
    - ...
    """
    global machine, logger
    T_outdoor, n_cmp = inputs[0], inputs[1]
    machine.set_operating_conditions(
        evp_m_dot_air=Q_(1500, 'kg / hr'),
        cnd_m_dot_air=Q_(2233.861, 'kg / hr'),
        evp_air_in=HumidAir(Tdb=Q_(24.0, 'degC'), RH=Q_(50, 'pct')),
        cnd_air_in=HumidAir(Tdb=T_outdoor, RH=Q_(50, 'pct')),
        dT_sh=Q_(10, 'K'),
        n_cmp=n_cmp
    )
    try:
        logger.info(
            f"Rating @ n_cmp = {n_cmp.to('1 / min'):~P.0f} "
            f"and T_out = {T_outdoor.to('degC'):~P.2f}"
        )
        machine.rate_by_minimization(
            T_evp_ini=Q_(5, 'degC'),
            T_cnd_ini=Q_(50, 'degC')
        )
    except Exception as err:
        traceback_str = traceback.format_exc()
        logger.error(f"{err}, {traceback_str}")
        return {
            'T_outdoor [degC]': T_outdoor.to('degC').m,
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
            'T_outdoor [degC]': T_outdoor.to('degC').m,
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


def _init_worker(log_file_path: pathlib.Path):
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
    log_file_path = pathlib.Path(
        "C:/Users/Tom/PycharmProjects/python-hvac/"
        "docs/examples/vapor_compression/log_files/vcm_ex07_2a.log"
    )
    # Delete previous log file if it exists:
    if log_file_path.exists():
        log_file_path.unlink()
    # Run analyze-performance-tasks in parallel:
    input_list = create_inputs()
    data = []
    with ProcessPoolExecutor(initializer=_init_worker, initargs=(log_file_path,)) as executor:
        start_time = time.perf_counter()
        for record in executor.map(analyze_performance, input_list):
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
    df.to_excel('./analysis_results/07_multi_analysis_2.ods')


if __name__ == '__main__':
    main()
