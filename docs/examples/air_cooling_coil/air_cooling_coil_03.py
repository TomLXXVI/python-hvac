"""
Given:
- An air-cooling coil in a space is cooling the space air.
- The mass flow rate of water through the air-cooling coil is known.
- The entering water temperature is known.
- The mass flow rate of air through the air-cooling coil is known.
- The sensible and latent space cooling load.

Calculate the space air state when steady-state operation (thermal equilibrium)
has been established (the thermal power absorbed by the water in the air-cooling
coil balances the total cooling load of the space).

Do this calculation for a range of space-cooling loads having different sensible
heat ratios. It's assumed that the space cooling load remains constant as the
space air temperature and humidity ratio evolve towards their steady-state
values.
"""
import warnings
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid, STANDARD_TEMPERATURE, STANDARD_PRESSURE
from hvac.heat_exchanger.recuperator.fintube.continuous_fin import (
    PlainFinTubeAirToWaterCounterFlowHeatExchanger as AirCoil
)
from hvac.air_conditioning import AirConditioningProcess
from hvac.charts import LineChart


warnings.filterwarnings('ignore', category=RuntimeWarning)


Q_ = Quantity
Water = Fluid('Water')
Air = Fluid('Air')
standard_air = Air(T=STANDARD_TEMPERATURE, P=STANDARD_PRESSURE)


air_coil = AirCoil(
    width=Q_(900, 'mm'),
    height=Q_(180, 'mm'),
    num_rows=3,
    pitch_trv=Q_(25.4, 'mm'),
    pitch_lon=Q_(22.0, 'mm'),
    d_i=Q_(8.422, 'mm'),
    d_o=Q_(10.2, 'mm'),
    t_fin=Q_(0.3302, 'mm'),
    fin_density=1 / Q_(3.175, 'mm'),
    num_circuits=2
)

# Fixed input air-cooling coil parameters:
air_V_dot = Q_(850.0, 'm ** 3 / hr')
air_m_dot = standard_air.rho * air_V_dot

water_in = Water(T=Q_(7, 'degC'), P=Q_(2, 'bar'))

water_V_dot = Q_(853.02, 'L / hr')
water_m_dot = water_in.rho * water_V_dot


def air_coil_fun(air_in: HumidAir) -> dict:
    """
    Given the state of air at the entry of the air cooling coil, returns the
    performance of the air cooling coil, while the other coil performance
    input parameters (the mass flow rate of air through the coil, the mass flow
    rate of water, and the entering water temperature) remain fixed.
    """
    air_coil.set_operating_conditions(
        air_in=air_in,
        water_in=water_in,
        air_m_dot=air_m_dot,
        water_m_dot=water_m_dot
    )
    result = air_coil.rate(eps_ini=0.5)
    return result


def space_fun(air_in: HumidAir, space_Q_dot: Quantity, space_SHR: Quantity) -> HumidAir:
    """
    Given the state of air entering the space and the space load, returns
    the state of air leaving the space and entering the air-cooling coil.
    """
    space = AirConditioningProcess(
        air_in=air_in,
        m_da=air_m_dot,
        Q=space_Q_dot,
        SHR=space_SHR
    )
    return space.air_out


def system_fun(args) -> HumidAir:
    """
    Combines the air-cooling coil and the space into a single thermal system,
    and determines the steady-state operation of this system for a given space
    load (i.e., without any intervention of a room thermostat).
    At steady-state operation, the space load and the cooling coil capacity have
    reached thermal equilibrium.
    In steady-state conditions, the state of space air entering the cooling coil
    becomes constant in time. Also, the state of air leaving the cooling coil,
    and the state of water leaving the cooling coil become constant.

    Returns
    -------
    The state of the space air when steady-state operation has been reached.
    """
    space_Q_dot, space_SHR = args

    solutions = {}
    i_max = 50
    i = 0
    dev_prev = None

    # Initial guess to start the loop:
    space_air = HumidAir(Tdb=Q_(26, 'degC'), RH=Q_(50, 'pct'))

    while i < i_max:
        # Calculate the coil performance for the current space air state:
        coil_perf = air_coil_fun(space_air)
        # Get the cooling capacity of the air-cooling coil:
        coil_Q_dot = coil_perf['Q_dot']
        # Get the state of air leaving the cooling coil:
        coil_air_out = coil_perf['air_out']
        # Check the deviation between the cooling capacity of the cooling coil
        # and the space cooling load:
        dev = abs(coil_Q_dot - space_Q_dot).to('W')
        if dev_prev is not None and abs(dev - dev_prev) < Q_(1.e-4, 'W'):
            break
        # Calculate the state of space air using the state of air leaving the
        # air-cooling coil:
        space_air = space_fun(coil_air_out, space_Q_dot, space_SHR)
        solutions[dev] = space_air
        dev_prev = dev
        i += 1
    else:
        print(f"No solution within tolerance after {i_max} iterations.")

    # Return the space air state for which the deviation is minimal:
    dev_min = min(dev for dev in solutions.keys())
    space_air = solutions[dev_min]
    return space_air


def main():
    """
    Calculate the space air state at thermal equilibrium for space-cooling loads
    in the range 1 to 7 kW having at a sensible heat ratio of 0.5 and 0.75.
    The steady-state space air temperature and space air humidity ratio are
    shown in two separate diagrams as a function of space load.

    We will use parallel programming to distribute the calculations over
    multiple processes.
    """
    space_Q_dot_list = [Q_(val, 'kW') for val in np.arange(1.0, 7.5, 0.5)]
    space_SHR_list = [0.5, 0.75]
    # To each total space cooling load in `space_Q_dot_list` we connect the SHRs
    # in `space_SHR_list`. Each tuple `(Q_dot, SHR)` in the two resulting lists
    # in `args_list` will serve as an input argument to a separate child process.
    # Thus, list `args_list` contains two lists. The first list contains the
    # tuples of all the space cooling loads in `space_Q_dot_list` combined with
    # the first SHR-value in `space_SHR_list`. The second list then contains the
    # tuples of all the space cooling loads combined with the second element in
    # `space_SHR_list`.
    args_list = [
        list(zip(space_Q_dot_list, [space_SHR] * len(space_Q_dot_list)))
        for space_SHR in space_SHR_list
    ]
    solutions = {}
    for space_SHR, args in zip(space_SHR_list, args_list):
        # For each list in `args_list` we run the function `system_fun` in
        # parallel processes with the tuples in that list. Each process returns
        # the state of the space air that corresponds with the given space
        # cooling load and SHR that was passed in a tuple to this process.
        print(f"Calculating space air state when SHR = {space_SHR}")
        space_air_list = []
        with ProcessPoolExecutor() as executor:
            start_time = time.perf_counter()
            for space_air in executor.map(system_fun, args):
                space_air_list.append(space_air)
            finish_time = time.perf_counter()
            print(f"Run time = {finish_time - start_time} seconds")
        space_T_list = [
            space_air.Tdb.to('degC').m
            for space_air in space_air_list
        ]
        space_W_list = [
            space_air.W.to('g / kg').m
            for space_air in space_air_list
        ]
        solutions[space_SHR] = (space_T_list, space_W_list)

    space_Q_dot_list = [
        space_Q_dot.to('kW').m
        for space_Q_dot in space_Q_dot_list
    ]

    T_diagram = LineChart()
    for space_SHR, sols in solutions.items():
        T_diagram.add_xy_data(
            label=f'space SHR = {space_SHR}',
            x1_values=space_Q_dot_list,
            y1_values=sols[0],
            style_props={'marker': 'o'}
        )
    T_diagram.add_legend()
    T_diagram.x1.add_title('total space load, kW')
    T_diagram.y1.add_title('space air temperature, °C')
    T_diagram.show()

    W_diagram = LineChart()
    for space_SHR, sols in solutions.items():
        W_diagram.add_xy_data(
            label=f'space SHR = {space_SHR}',
            x1_values=space_Q_dot_list,
            y1_values=sols[1],
            style_props={'marker': 'o'}
        )
    W_diagram.add_legend()
    W_diagram.x1.add_title('total space load, kW')
    W_diagram.y1.add_title('space humidity ratio, g/kg')
    W_diagram.show()


if __name__ == '__main__':
    main()
