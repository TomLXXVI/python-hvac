"""
EXAMPLE 1
---------
DYNAMIC HEAT TRANSFER THROUGH A CONCRETE WALL
Example of the dynamic conduction heat transfer through a concrete wall using a
linear thermal network. Comparison with the assumption of steady-state
conduction heat transfer (which neglects the effects of thermal capacity).

Uses the implementation of a linear thermal network in module:
`hvac.cooling_load_calc.core.thermal_model.linear_thermal_network.py`.
"""
import numpy as np
import pandas as pd

from hvac import Quantity
from hvac.cooling_load_calc.core.thermal_models.linear_thermal_network import (
    TemperatureNode,
    LinearThermalNetwork
)
from hvac.charts.matplotlibwrapper import LineChart


Q_ = Quantity


# Thermal properties of the concrete wall.
rho = Q_(2400, 'kg / m**3')
c = Q_(840, 'J / (kg * K)')
k = Q_(1.87, 'W / (m * K)')


# Surface area and thickness of the concrete wall.
A_wall = Q_(30, 'm**2')
wall_thick = Q_(10, 'cm')


# Average outdoor temperature:
T_out_avg = Q_(25.0, 'degC')


# Fixed indoor air temperature:
T_ind = Q_(22.0, 'degC')


def T_fun1(t_sec: float) -> Quantity:
    # harmonic fluctuating temperature at the left side (exterior side) of the
    # concrete wall
    period = 24 * 3600  # one day = 24 hours * 3600 sec/hr
    ampl = Q_(5.0, 'K')
    T = T_out_avg.to('K') + ampl * np.sin(2 * np.pi * t_sec / period)
    return T


def T_fun2(_) -> Quantity:
    # constant temperature at the right side of the concrete wall
    return T_ind


def create_linear_thermal_network(num_nodes: int) -> LinearThermalNetwork:
    """Creates a linear thermal network of the concrete wall with the number
    of nodes given by `num_nodes`.
    """
    # We will divide the concrete wall in a number of layers. Each layer will
    # be represented by a temperature node in the linear thermal network.
    layer_thick = wall_thick / num_nodes

    C = rho * layer_thick.to('m') * c  # unit capacitance of 1 layer
    R = layer_thick.to('m') / k        # unit thermal resistance of 1 layer

    nodes = []
    # Create the first node adjacent to the outdoor environment:
    first_node = TemperatureNode.create(
        ID='NODE 1',
        R1=R/2,
        R2=R,
        C=C,
        A=A_wall,
        T_fun1=T_fun1
    )
    nodes.append(first_node)
    # Create the intermediate nodes:
    i = 2
    i_max = num_nodes - 1
    while i <= i_max:
        n = TemperatureNode.create(
            ID=f'NODE {i}',
            R1=R,
            R2=R,
            C=C,
            A=A_wall
        )
        nodes.append(n)
        i += 1
    # Create the last node adjacent to the indoor environment:
    last_node = TemperatureNode.create(
        ID=f'NODE {num_nodes}',
        R1=R,
        R2=R/2 + Q_(0.13, 'K * m**2 / W'),
        C=C,
        A=A_wall,
        T_fun2=T_fun2
    )
    nodes.append(last_node)
    # Create the linear thermal network, being a list of temperature nodes:
    nw = LinearThermalNetwork(nodes)
    return nw


def steady_state_conduction_heat_transfer(t_sec: float):
    """Calculates the conduction heat transfer through the concrete wall at
    time moment `t_sec` from time zero based on the equation for steady-state
    conduction heat transfer.
    """
    R_wall_unit = wall_thick.to('m') / k.to('W / (K * m)')
    R_wall = R_wall_unit.to('K * m**2 / W') / A_wall.to('m**2')
    T_out = T_fun1(t_sec)
    Q_dot = (T_out.to('K') - T_ind.to('K')) / R_wall.to('K / W')
    return Q_dot


def main():
    # Set the number of nodes (or layers) the concrete wall is composed of:
    num_nodes = 10

    # Set the time span for the calculation of the node temperatures in the
    # course of time and the time span for which we want to get the results
    # back:
    num_days = 5
    eval_between_days = (num_days - 1, num_days)
    time_step = 3600  # seconds

    # Create and solve the linear thermal network:
    nw = create_linear_thermal_network(num_nodes)
    nw.solve(
        time_span=[0.0, num_days * 24.0 * 3600.0],
        T_node_init_seq=[T_ind for _ in range(num_nodes)],
        time_eval=np.arange(
            eval_between_days[0] * 24.0 * 3600.0,
            eval_between_days[1] * 24.0 * 3600.0,
            time_step
        ),
        continuous=True
    )
    T_node_table = nw.get_T_node_table()
    Q_dot_table = nw.get_Q_dot_table()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print('TEMPERATURE NODE TABLE:')
        print(T_node_table)
        print('HEAT FLOW TABLE:')
        print(Q_dot_table)

    # Get the conduction heat transfer through the concrete wall based on the
    # equation for steady-state conduction heat transfer:
    Q_dot_ss = [steady_state_conduction_heat_transfer(t) for t in nw.time_axis]

    # Plot the heat flow at the exterior and interior side of the wall:
    chart = LineChart()
    chart.add_xy_data(
        label='exterior heat flow',
        x1_values=Q_dot_table.index,
        y1_values=Q_dot_table['NODE 1']
    )
    chart.add_xy_data(
        label='interior heat flow',
        x1_values=Q_dot_table.index,
        y1_values=Q_dot_table['OUT']
    )
    chart.add_xy_data(
        label='steady-state heat flow',
        x1_values=Q_dot_table.index,
        y1_values=[Q.to('W').m for Q in Q_dot_ss]
    )
    chart.add_legend()
    chart.x1.add_title('time, hr')
    chart.y1.add_title('heat flow, W')
    chart.show()


if __name__ == '__main__':
    main()
