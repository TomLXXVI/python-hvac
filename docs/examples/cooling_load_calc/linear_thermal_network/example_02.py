"""
EXAMPLE 2
---------
DYNAMIC HEAT TRANSFER THROUGH A CONCRETE WALL

Uses the implementation of a linear thermal network in module:
`hvac.cooling_load_calc.core.thermal_model.exterior_building_element.py`.

This is the implementation which is further used in subpackage
`hvac.cooling_load_calc`.
"""
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.cooling_load_calc.core.thermal_models import (
    ExteriorSurfaceNode, BuildingMassNode, InteriorSurfaceNode,
    ExteriorBuildingElementLTN
)
from hvac.charts import LineChart

Q_ = Quantity


# Thermal properties of the concrete wall:
rho = Q_(2300, 'kg / m**3')
c = Q_(840, 'J / (kg * K)')
k = Q_(1.87, 'W / (m * K)')


# Surface area of the concrete wall:
A_wall = Q_(30, 'm**2')
layer_thick = Q_(1, 'cm')


def T_ext(t_sec: float) -> Quantity:
    # Harmonic fluctuating temperature at the left side (exterior side) of the
    # concrete wall:
    period = 24 * 3600  # one day = 24 hours * 3600 s/hr
    T_out_avg = Q_(25.0, 'degC')
    ampl = Q_(5.0, 'K')
    T = T_out_avg.to('K') + ampl * np.sin(2 * np.pi * t_sec / period)
    return T


def T_zone(_) -> Quantity:
    # Constant temperature at the right side of the concrete wall:
    T_ind = Q_(22.0, 'degC')
    return T_ind


def create_linear_thermal_network(num_nodes: int) -> ExteriorBuildingElementLTN:
    esn = ExteriorSurfaceNode.create(
        ID='ESN',
        C=c * rho * layer_thick,
        A=A_wall,
        R1=layer_thick / k / 2,
        R2=layer_thick / k
    )
    bmn_lst = [
        BuildingMassNode.create(
            ID=f'BMN {i + 1}',
            C=c * rho * layer_thick,
            A=A_wall,
            R1=layer_thick / k,
            R2=layer_thick / k
        ) for i in range(num_nodes - 2)
    ]
    isn = InteriorSurfaceNode.create(
        ID='ISN',
        A=A_wall,
        R1=layer_thick / k / 2,
        R2=Q_(0.13, 'K * m**2 / W')
    )
    nw = ExteriorBuildingElementLTN([esn])
    nw.extend(bmn_lst)
    nw.append(isn)
    return nw


def steady_state_conduction_heat_transfer(t_sec: float, num_nodes: int):
    """Calculates the conduction heat transfer through the concrete wall at
    time moment `t_sec` based on the equation for steady-state conduction heat
    transfer.
    """
    wall_thick = num_nodes * layer_thick
    R_wall_unit = wall_thick.to('m') / k.to('W / (K * m)')
    R_wall = R_wall_unit.to('K * m**2 / W') / A_wall.to('m**2')
    T_out = T_ext(t_sec)
    Q_dot = (T_out.to('K') - T_zone(t_sec).to('K')) / R_wall.to('K / W')
    return Q_dot


def main():
    num_nodes = 10
    nw = create_linear_thermal_network(num_nodes)
    init_values = [[Q_(22.0, 'degC')] * len(nw) for _ in range(2)]
    dt_hr = 1 / 4
    num_steps = int(24 / dt_hr)
    nw.solve(
        num_steps=num_steps,
        dt_hr=dt_hr,
        T_ext=T_ext,
        T_zone=T_zone,
        init_values=init_values,
        num_cycles=5
    )
    T_node_table = nw.get_node_temperatures()
    Q_dot_table = nw.get_heat_flows()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print('NODE TEMPERATURE TABLE:')
        print(T_node_table)
        print('HEAT FLOW TABLE:')
        print(Q_dot_table)

    # Get the conduction heat transfer through the concrete wall based on the
    # equation for steady-state conduction heat transfer:
    Q_dot_ss = [
        steady_state_conduction_heat_transfer(k * dt_hr * 3600, num_nodes)
        for k in range(num_steps)
    ]

    # Plot the heat flow at the exterior and interior side of the wall:
    chart = LineChart()
    chart.add_xy_data(
        label='exterior heat flow',
        x1_values=Q_dot_table.index,
        y1_values=Q_dot_table['ESN']
    )
    chart.add_xy_data(
        label='interior heat flow',
        x1_values=Q_dot_table.index,
        y1_values=Q_dot_table['ISN']
    )
    chart.add_xy_data(
        label='steady-state heat flow',
        x1_values=Q_dot_table.index,
        y1_values=[Q_dot.to('W').m for Q_dot in Q_dot_ss]
    )
    chart.add_legend()
    chart.x1.add_title('time, hr')
    chart.y1.add_title('heat flow, W')
    chart.show()


if __name__ == '__main__':
    main()
