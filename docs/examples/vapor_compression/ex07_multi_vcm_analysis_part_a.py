"""
STEADY-STATE ANALYSIS OF A SINGLE-STAGE VAPOR COMPRESSION MACHINE

A model of a single-stage vapor compression machine is used to analyze its
steady-state operation at given operating conditions:
    - the mass flow rate of air through the evaporator
    - the state of air entering the evaporator
    - the mass flow rate of air through the condenser
    - the state of air entering the condenser
    - the speed of the compressor
    - the degree of refrigerant superheating set on the expansion device.

The degree of refrigerant superheating set on the expansion device is considered
to be a fixed property of the machine, and it is specified when creating the
machine model (see function `create_machine` below).

The analysis will be made at different compressor speeds using multiprocessing
(with `ProcessPoolExecutor`). The other operating conditions remain unchanged.
(Note: each child process created by the `ProcessPoolExecutor` will run this
entire script.)

The analysis or rating routine tries to find the evaporation and condensation
temperature for which the mass flow rate of the refrigerant displaced by the
compressor equals the mass flow rate of the refrigerant let trough by the
expansion device. Once the evaporation and condensation temperature are
determined, the full operating state of the machine is known.
The analysis starts with an initial guess for the evaporation and condensation
temperature. Using a least-squares algorithm, the evaporation and condensation
temperature are iteratively searched for until the deviation between the two
mass flow rates becomes minimal.
"""
import time
import warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from hvac import Quantity
from hvac.fluids import Fluid, HumidAir, CoolPropWarning
from hvac.vapor_compression import VariableSpeedCompressor
from hvac.vapor_compression.machine import (
    SingleStageVaporCompressionMachine,
    Output,
    logger
)
from hvac.heat_exchanger.recuperator.fintube.continuous_fin import (
    PlainFinTubeCounterFlowAirEvaporator,
    PlainFinTubeCounterFlowAirCondenser
)
from hvac.logging import ModuleLogger


# Turn off any runtime warnings:
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Turn off warnings coming from module `hvac.fluids.fluid`:
warnings.filterwarnings('ignore', category=CoolPropWarning)

# Create a shortcut for the `Quantity`-constructor:
Q_ = Quantity


# ==============================================================================
def create_machine() -> SingleStageVaporCompressionMachine:

    # Define refrigerant:
    R134a = Fluid('R134a')

    # Create the compressor model
    compressor = VariableSpeedCompressor(
        coeff_file=Path("./compressor_data/VTZ054-G_R134a.csv"),
        refrigerant=R134a,
        units={'n': 'rps'}
    )

    # Create the evaporator model:
    evaporator = PlainFinTubeCounterFlowAirEvaporator(
        W_fro=Q_(0.731, 'm'),
        H_fro=Q_(0.244, 'm'),
        N_rows=3,
        S_trv=Q_(25.4, 'mm'),
        S_lon=Q_(22.0, 'mm'),
        D_int=Q_(8.422, 'mm'),
        D_ext=Q_(10.2, 'mm'),
        t_fin=Q_(0.3302, 'mm'),
        N_fin=1 / Q_(3.175, 'mm')
    )

    # Create the condenser model:
    condenser = PlainFinTubeCounterFlowAirCondenser(
        W_fro=Q_(1.003, 'm'),
        H_fro=Q_(0.334, 'm'),
        N_rows=5,
        S_trv=Q_(25.4, 'mm'),
        S_lon=Q_(22.0, 'mm'),
        D_int=Q_(8.422, 'mm'),
        D_ext=Q_(10.2, 'mm'),
        t_fin=Q_(0.3302, 'mm'),
        N_fin=1 / Q_(3.175, 'mm')
    )

    # Configure the single-stage vapor compression machine:
    machine = SingleStageVaporCompressionMachine(
        evaporator,
        condenser,
        compressor,
        dT_sh=Q_(5, 'K'),
        n_cmp_min=Q_(2100, 'rpm'),
        n_cmp_max=Q_(5400, 'rpm')
    )
    return machine


# ==============================================================================
def analyze_performance(n_cmp: Quantity) -> Output:
    # This function (the task) will run in parallel in different processes.

    # Create machine:
    machine = create_machine()

    # Set the other, fixed operating conditions of the machine:
    # - mass flow rate of air through evaporator
    # - mass flow rate of air through condenser
    # - state of air entering evaporator
    # - state of air entering condenser

    evp_air_in = HumidAir(Tdb=Q_(24.0, 'degC'), RH=Q_(50, 'pct'))
    evp_m_dot_air = Q_(1500.0, 'kg / hr')

    cnd_air_in = HumidAir(Tdb=Q_(35.0, 'degC'), RH=Q_(30, 'pct'))
    cnd_m_dot_air = Q_(3216.315, 'kg / hr')

    # Run rating routine to determine steady-state machine performance:
    logger.info(
        f"Rating @ n_cmp = {n_cmp.to('rpm'):~P.0f}"
    )

    output = machine.rate(
        evp_air_in=evp_air_in,
        evp_air_m_dot=evp_m_dot_air,
        cnd_air_in=cnd_air_in,
        cnd_air_m_dot=cnd_m_dot_air,
        n_cmp=n_cmp,
        # T_evp_ini=Q_(5.0, 'degC'),  # optional initial guess by user
        # T_cnd_ini=Q_(50.0, 'degC')  # optional initial guess by user
    )
    return output


def init_worker(log_file_path: Path):
    # In each worker process, add a file handler to the imported logger from
    # module `machine`:
    file_handler = ModuleLogger.create_file_handler(log_file_path)
    logger.addHandler(file_handler)


# ==============================================================================
def main():
    # Set path to log file:
    log_file_path = Path(
        "C:/Users/Tom/PycharmProjects/python-hvac/"
        "docs/examples/vapor_compression/log_files/vcm_ex07_1.log"
    )

    # Delete the previous log file if it exists:
    if log_file_path.exists():
        log_file_path.unlink()

    # Run `analyze_performance` with different compressor speeds parallel in
    # multiple processes with `ProcessPoolExecutor`. Before a process is
    # started, it is first initialized by calling function `_init_worker`,
    # passing the `log_file_path` to it to write any log messages to this file:
    n_cmp_rng = [Q_(n, 'rpm') for n in range(2100, 5700, 300)]
    outputs = []
    with ProcessPoolExecutor(
        initializer=init_worker,
        initargs=(log_file_path,)
    ) as executor:
        start_time = time.perf_counter()
        for output in executor.map(analyze_performance, n_cmp_rng):
            outputs.append(output)
        finish_time = time.perf_counter()
        print(f"Run time = {finish_time - start_time} seconds")

    # Put the records in a Pandas DataFrame object:
    df = pd.DataFrame([output.to_dict() for output in outputs])
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print(df)

    # Save dataframe to spreadsheet file:
    df.to_excel('./analysis_results/07_multi_analysis_1.ods')

    # Sort the log file by process ID:
    ModuleLogger.sort_log_by_process_id(log_file_path)


if __name__ == '__main__':
    main()
    # Function `main` is only called by the main process.
    # Within the function `main` the `ProcessPoolExecutor` launches separate
    # child processes that will run this same script (and so will import
    # everything mentioned at the head of this script), but without calling
    # function `main`; only functions `init_worker` and `analyze_performance`
    # will be called from this script.
