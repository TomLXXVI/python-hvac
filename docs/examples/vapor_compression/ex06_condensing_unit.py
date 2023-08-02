"""Analysis of the refrigeration capacity of a condensing unit.

A condensing unit is a combination of a compressor and condenser. In the
analysis of a condensing unit the steady-state performance is determined for
a range of condensing temperatures, while the evaporation temperature remains
constant. Other operating conditions remain fixed (compressor speed, degree
of superheat set on expansion device, state and mass flow rate of air entering
condenser).

For the implementation of the analysis multiprocessing is applied using the
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
from hvac.heat_transfer.heat_exchanger.fin_tube import air_condenser, air_evaporator
from hvac.vapor_compression import VariableSpeedCompressor, FixedSpeedCompressor

warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=CoolPropWarning)
warnings.filterwarnings('ignore', category=air_condenser.CondenserWarning)

logger = ModuleLogger.get_logger(__name__)

Q_ = Quantity

Evaporator = air_evaporator.rating.PlainFinTubeCounterFlowEvaporator
Condenser = air_condenser.rating.PlainFinTubeCounterFlowCondenser
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
        """Creates instance of `CondensingUnit` class.

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
        self.compressor.Te = T_evp
        self.compressor.Tc = T_cnd
        self.compressor.speed = self.n_cmp
        self.compressor.dT_sh = self.dT_sh
        # Note: passing the degree of superheat to the compressor allows the
        # compressor to determine the state of the suction gas entering the
        # compressor and the state of the discharge gas leaving the compressor.
        # Determine performance of condenser:
        self.condenser(
            m_dot_air=self.m_dot_air,
            m_dot_rfg=self.compressor.m_dot,
            air_in=self.air_in,
            rfg_in=self.compressor.discharge_gas
        )
        return self.compressor.Wc_dot, self.condenser.Q_dot


# Evaporator model:
evaporator = Evaporator(
    L1=Q_(0.731, 'm'),
    L3=Q_(0.244, 'm'),
    N_r=3,
    S_t=Q_(25.4, 'mm'),
    S_l=Q_(22.0, 'mm'),
    D_i=Q_(8.422, 'mm'),
    D_o=Q_(10.2, 'mm'),
    t_f=Q_(0.3302, 'mm'),
    N_f=1 / Q_(3.175, 'mm')
)

# Condenser model:
condenser = Condenser(
    L1=Q_(0.999, 'm'),
    L3=Q_(0.333, 'm'),
    N_r=5,
    S_t=Q_(25.4, 'mm'),
    S_l=Q_(22.0, 'mm'),
    D_i=Q_(8.422, 'mm'),
    D_o=Q_(10.2, 'mm'),
    t_f=Q_(0.3302, 'mm'),
    N_f=1 / Q_(3.175, 'mm')
)

# Compressor model:
compressor = VariableSpeedCompressor(
    coeff_file=Path("./compressor_data/VTZ054-G_R134a.csv"),
    refrigerant_type=R134a,
    units={'m_dot': 'kg / hr', 'speed': '1 / s'}
)

# Condensing unit model:
condensing_unit = CondensingUnit(
    compressor=compressor,
    condenser=condenser,
    n_cmp=Q_(3000, '1 / min'),
    dT_sh=Q_(10, 'K'),
    air_in=HumidAir(Tdb=Q_(35, 'degC'), RH=Q_(30, 'pct')),
    m_dot_air=Q_(3000, 'kg / hr')
)


def task(T_evp: Quantity, T_cnd: Quantity) -> Quantity:
    """Runs the analysis of the condensing unit for the given evaporating
    temperature `T_evp` and condensing temperature `T_cnd` and returns the
    refrigeration capacity of the condensing unit.
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
    T_evp = 3.0  # degC
    # Set range of condensation temperatures:
    T_cnd_rng = Q_(np.arange(40.0, 62.0, 2.0), 'degC')
    # Create an array of the same evaporation temperature having the same length
    # as the range of condensation temperatures:
    T_evp_rng = Q_(np.array(list(itertools.repeat(T_evp, len(T_cnd_rng)))), 'degC')
    # Create a list to store the refrigeration capacities of the condensing unit:
    Q_evp_rng = []

    # For each pair of evaporation and condensation temperature, run the analysis
    # of the condensing unit in a separate process:
    with ProcessPoolExecutor() as executor:
        for r in executor.map(task, T_evp_rng, T_cnd_rng):
            Q_evp_rng.append(r) # collect the refrigeration capacities

    # Process the results --> We want to see the refrigeration capacity as a
    # function of the condensation temperature:
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
