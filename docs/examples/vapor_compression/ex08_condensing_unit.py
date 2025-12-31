"""ANALYSIS OF THE REFRIGERATION CAPACITY OF A CONDENSING UNIT.

A condensing unit is a combination of a compressor and condenser. In the
analysis of a condensing unit, the steady-state performance is determined for
a range of condensing temperatures, while the evaporation temperature remains
constant. Other operating conditions remain fixed (compressor speed, degree
of superheat set on the expansion device, state and mass flow rate of air
entering the condenser).

For the implementation of the analysis, multiprocessing is applied using the
`ProcessPoolExecutor` class from `concurrent.futures`. The analysis for each
pair of evaporating and condensing temperatures runs parallel in a separate
process.
"""
import warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import itertools
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.charts import LineChart
from hvac.fluids import Fluid, HumidAir, CoolPropWarning
from hvac.heat_exchanger.recuperator.fintube.continuous_fin import (
    PlainFinTubeCounterFlowAirCondenser
)
from hvac.vapor_compression import (
    VariableSpeedCompressor,
    FixedSpeedCompressor
)

warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=CoolPropWarning)

logger = ModuleLogger.get_logger(__name__)

Q_ = Quantity

Condenser = PlainFinTubeCounterFlowAirCondenser
Compressor = VariableSpeedCompressor | FixedSpeedCompressor
R134a = Fluid('R134a')


class CondensingUnit:
    """Model class of a condensing unit."""

    def __init__(
        self,
        compressor: Compressor,
        condenser: Condenser,
        n_cmp: Quantity,
        dT_sh: Quantity,
        air_in: HumidAir,
        m_dot_air: Quantity
    ) -> None:
        """Creates an instance of the `CondensingUnit` class.

        Parameters
        ----------
        compressor:
            Compressor model.
        condenser:
            Condenser model.
        n_cmp:
            Compressor speed at which the condensing unit will be analyzed.
        dT_sh:
            Degree of superheat set on the expansion device of the vapor
            compression machine.
        air_in:
            State of air entering the condenser.
        m_dot_air:
            Mass flow rate of air through the condenser.
        """
        self.compressor = compressor
        self.condenser = condenser
        self.n_cmp = n_cmp
        self.dT_sh = dT_sh
        self.air_in = air_in
        self.m_dot_air = m_dot_air

    def rate(self, T_evp: Quantity, T_cnd: Quantity) -> tuple[Quantity, Quantity]:
        """Determines the performance of the condensing unit for the given
        operating conditions and for the given evaporation and condensation
        temperature.

        Parameters
        ----------
        T_evp:
            Evaporation temperature.
        T_cnd:
            Condensing temperature.

        Returns
        -------
        The compressor power to the refrigerant and the heat rejection rate
        at the condenser for the given operating conditions.
        """
        # Set compressor parameters:
        self.compressor.T_evp = T_evp
        self.compressor.T_cnd = T_cnd
        self.compressor.speed = self.n_cmp
        self.compressor.dT_sh = self.dT_sh
        # Note: passing the degree of superheat to the compressor allows the
        # compressor to determine the state of the suction gas entering the
        # compressor and the state of the discharge gas leaving the compressor.
        # Determine performance of condenser:
        self.condenser.solve(
            air_m_dot=self.m_dot_air,
            rfg_m_dot=self.compressor.m_dot,
            air_in=self.air_in,
            rfg_in=self.compressor.discharge_gas
        )
        return self.compressor.W_dot, self.condenser.Q_dot


# Condenser model:
condenser = Condenser(
    W_fro=Q_(0.999, 'm'),
    H_fro=Q_(0.333, 'm'),
    N_rows=5,
    S_trv=Q_(25.4, 'mm'),
    S_lon=Q_(22.0, 'mm'),
    D_int=Q_(8.422, 'mm'),
    D_ext=Q_(10.2, 'mm'),
    t_fin=Q_(0.3302, 'mm'),
    N_fin=1 / Q_(3.175, 'mm')
)

# Compressor model:
compressor = VariableSpeedCompressor(
    coeff_file=Path("./compressor_data/VTZ054-G_R134a.csv"),
    refrigerant=R134a,
    units={'n': 'rps'}
)

# Condensing unit model:
condensing_unit = CondensingUnit(
    compressor=compressor,
    condenser=condenser,
    n_cmp=Q_(3000, 'rpm'),
    dT_sh=Q_(10, 'K'),
    air_in=HumidAir(Tdb=Q_(35, 'degC'), RH=Q_(30, 'pct')),
    m_dot_air=Q_(3000, 'kg / hr')
)


def task(T_evp: Quantity, T_cnd: Quantity) -> Quantity:
    """Runs the analysis of the condensing unit for the given evaporating
    temperature `T_evp` and condensing temperature `T_cnd` and returns the
    refrigeration capacity of the condensing unit, being the heat rejection
    rate in the condenser minus the mechanical power added to the refrigerant
    by the compressor.
    """
    global condensing_unit
    logger.info(f"Running task with T_evp {T_evp:~P.2f} and T_cnd {T_cnd:~P.2f}")
    try:
        W_cmp, Q_cnd = condensing_unit.rate(T_evp, T_cnd)
    except Exception:
        return Q_(float('nan'), 'kW')
    else:
        Q_evp = Q_cnd - W_cmp
        return Q_evp.to('kW')


if __name__ == '__main__':
    # Set evaporation temperature:
    T_evp = 5.0  # degC
    # Set range of condensation temperatures:
    T_cnd_rng = Q_(np.arange(54.0, 72.0, 2.0), 'degC')
    # Create an array of the same evaporation temperature having the same length
    # as the range of condensation temperatures:
    T_evp_rng = Q_(np.array(list(itertools.repeat(T_evp, len(T_cnd_rng)))), 'degC')
    # Create a list to store the refrigeration capacities of the condensing unit:
    Q_evp_rng = []

    # For each pair of evaporation and condensation temperature, run the analysis
    # of the condensing unit in a separate process:
    with ProcessPoolExecutor() as executor:
        for r in executor.map(task, T_evp_rng, T_cnd_rng):
            Q_evp_rng.append(r)  # collect the refrigeration capacities

    # Process the results --> We want to see the refrigeration capacity of the
    # condensing unit as a function of the condensation temperature:
    d = {
        'T_cnd': [T_cnd.to('degC').m for T_cnd in T_cnd_rng],
        'Q_evp': [Q_evp.to('kW').m for Q_evp in Q_evp_rng]
    }

    df = pd.DataFrame(d)
    print(df)

    fig = LineChart()
    fig.add_xy_data(
        label=f"Q_evp at T_evp = {T_evp} °C",
        x1_values=df['T_cnd'],
        y1_values=df['Q_evp'],
        style_props={'marker': 'o', 'linestyle': 'none'}
    )
    fig.x1.add_title('T_cnd, °C')
    fig.y1.add_title('Q_evp, kW')
    fig.show()
