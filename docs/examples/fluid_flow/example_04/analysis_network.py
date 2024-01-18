"""ANALYSIS OF A PIPE NETWORK
-----------------------------
The pipe network is part of a HVAC system equipped with fan coil units that are
operated in cooling mode. In the PDF-document `network_scheme.pdf` a schematic
diagram of the network is shown with all the details about the pipe network.
The circulation pump is a GRUNDFOS MAGNA3 32-120 F. From the GRUNDFOS website
the pump curves at different operating speeds can be retrieved and downloaded
in an Excel file (see file `pump-curve-data.xlsx`). Using this data, a second-
order polynomial of the pump curve can be derived. This is done in the script
`pump_curve.py`. The coefficients a0, a1, and a2 of this polynomial
(`dP_pump = a0 + a1 * V_dot + a2 * V_dot**2`) are needed to add a pump to a pipe
in the network. These coefficients are entered in the CSV-file
`input-analysis-csv` which is used to configure the pipe network. More
information about this configuration file can be found in the docstring of
the method `load_from_csv()` of the `Network` class in module
`hvac.fluid_flow.network.py`.
"""
import warnings
import pandas as pd
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning

warnings.filterwarnings('ignore', category=CoolPropWarning)

from hvac.fluid_flow import (
    PipeSchedule,
    PipeNetwork,
    Pipe,
    Circular,
    SystemCurve,
    PumpCurve,
    plot_curves
)
from hvac.fluid_flow.fittings import pipe as fittings
from hvac.charts import LineChart

# ------------------------------------------------------------------------------
Q_ = Quantity

# Set the average state of the water (a cooling water regime in the fan coil
# units of 7 °C in / 12 °C out and a system pressure of 1.5 bar are assumed):
Water = Fluid('Water')
water = Water(T=Q_((12 + 7) / 2, 'degC'), P=Q_(1.5, 'bar'))
# ------------------------------------------------------------------------------


def main():
    # Create the pipe network:
    network = create_network()

    # Now we can change the valve openings of the control valves to see the
    # effect of this on the working point of the pump and system. The control
    # valves are installed upstream of the fan coil units in the cross-over
    # pipes of the network. The IDs of these cross-overs are:
    cross_over_lst = [
        'D8R8', 'D7R7', 'D6R6', 'D5R5',
        'D4R4', 'D3R3', 'D9R9', 'D10R10'
    ]
    # Set the valve opening of each control valve:
    for ID in cross_over_lst:
        network.set_control_valve_opening(ID, percent_open=100)

    # Now, we can analyze the pipe network with the method of Hardy Cross. This
    # will determine the volume flow rate in each pipe of the network:
    network.analyze(tolerance=Q_(0.01, 'kPa'), i_max=100)

    # Print the working point of the pipe network:
    print(
        "working point: ",
        f"{network.volume_flow_rate.to('L/s'):~P.3f} ",
        f"{network.total_pressure_difference.to('bar'):~P.1f}"
    )

    # Draw the system curve of the network and the pump curve:
    diagram = create_curves(network)
    diagram.show()


def create_network() -> PipeNetwork:
    # Create the pipe schedule used for the pipes in the network:
    pipe_schedule, pipe_roughness = create_pipe_schedule()

    # Create a `PipeNetwork` object:
    network = PipeNetwork.create(
        ID='first-floor',
        fluid=water,
        wall_roughness=pipe_roughness,
        schedule=pipe_schedule,
        start_node_ID='D2',  # where fluid enters the network
        end_node_ID='R2',    # where fluid leaves the network
        units={
            'length': 'cm',
            'volume_flow_rate': 'L / s',
            'pressure': 'kPa'
        }
    )
    # Load the pipe network configuration file into the `PipeNetwork` object:
    network.load_from_csv('input-analysis.csv')

    # Add fittings to the pipes (see network diagram):
    # Supply pipes:
    add_fittings_D2D3(network)
    add_fittings_D3D4(network)
    add_fittings_D4D5(network)
    add_fittings_D5D6(network)
    add_fittings_D6D7(network)
    add_fittings_D2D9(network)
    add_fittings_D9D10(network)
    # Return pipes:
    add_fittings_R7R6(network)
    add_fittings_R6R5(network)
    add_fittings_R5R4(network)
    add_fittings_R4R3(network)
    add_fittings_R3R2(network)
    add_fittings_R10R9(network)
    add_fittings_R9R2(network)
    # Cross-overs:
    add_fittings_to_cross_over(
        network=network,
        pipe=network.conduits['D8R8'],
        supply_pipe=network.conduits['D7D8'],
        return_pipe=network.conduits['R8R7'],
        dP_fcu=Q_(11.2, 'kPa'),
        num_elbows=4,
        has_tees=False
    )
    add_fittings_to_cross_over(
        network=network,
        pipe=network.conduits['D7R7'],
        supply_pipe=network.conduits['D6D7'],
        return_pipe=network.conduits['R7R6'],
        dP_fcu=Q_(9.7, 'kPa')
    )
    add_fittings_to_cross_over(
        network=network,
        pipe=network.conduits['D6R6'],
        supply_pipe=network.conduits['D5D6'],
        return_pipe=network.conduits['R6R5'],
        dP_fcu=Q_(9.7, 'kPa')
    )
    add_fittings_to_cross_over(
        network=network,
        pipe=network.conduits['D5R5'],
        supply_pipe=network.conduits['D4D5'],
        return_pipe=network.conduits['R5R4'],
        dP_fcu=Q_(9.7, 'kPa')
    )
    add_fittings_to_cross_over(
        network=network,
        pipe=network.conduits['D4R4'],
        supply_pipe=network.conduits['D3D4'],
        return_pipe=network.conduits['R4R3'],
        dP_fcu=Q_(9.7, 'kPa'),
        num_elbows=0
    )
    add_fittings_to_cross_over(
        network=network,
        pipe=network.conduits['D3R3'],
        supply_pipe=network.conduits['D2D3'],
        return_pipe=network.conduits['R3R2'],
        dP_fcu=Q_(14.0, 'kPa'),
        num_elbows=0
    )
    add_fittings_to_cross_over(
        network=network,
        pipe=network.conduits['D9R9'],
        supply_pipe=network.conduits['D2D9'],
        return_pipe=network.conduits['R9R2'],
        dP_fcu=Q_(15.3, 'kPa'),
        num_elbows=0
    )
    add_fittings_to_cross_over(
        network=network,
        pipe=network.conduits['D10R10'],
        supply_pipe=network.conduits['D9D10'],
        return_pipe=network.conduits['R10R9'],
        dP_fcu=Q_(9.1, 'kPa'),
        num_elbows=0,
        has_tees=False
    )
    # Add balancing and control valves to the cross-overs:
    add_valves(network)
    return network


def create_pipe_schedule() -> tuple[PipeSchedule, Quantity]:
    """Returns the pipe schedule for the pipe network and the pipe wall
    roughness. The data is based on GEBERIT MAPRESS CARBON STEEL pipe.
    """
    pipe_records = [
        (10, 12, 1.2),  # (nominal diameter, outside diameter, wall thickness)
        (12, 15, 1.2),
        (15, 18, 1.2),
        (20, 22, 1.5),
        (25, 28, 1.5),
        (32, 35, 1.5),
        (40, 42, 1.5),
        (50, 54, 1.5),
        (65, 76.1, 2),
        (80, 88.9, 2),
        (100, 108, 2)
    ]
    D_nom, D_ext, t_wall = zip(*pipe_records)
    D_nom = Q_(D_nom, 'mm')
    D_ext = Q_(D_ext, 'mm')
    t_wall = Q_(t_wall, 'mm')
    D_int = D_ext - 2 * t_wall
    # Create a `PipeSchedule` object and set the unit in which the diameters
    # and pipe wall thickness will be entered to the lookup-table:
    pipe_schedule = PipeSchedule(unit='mm')
    # Create the lookup-table, which is a Pandas DataFrame (note: mind to use the
    # correct keys for the column names):
    pipe_schedule.lookup_table = pd.DataFrame({
        'D_nom': D_nom.to('mm').m,
        'D_ext': D_ext.to('mm').m,
        't': t_wall.to('mm').m,
        'D_int': D_int.to('mm').m
    })
    # Column 'D_nom' must be set as the index of the dataframe:
    pipe_schedule.lookup_table.set_index(['D_nom'], inplace=True)
    # Set the absolute pipe wall roughness:
    pipe_roughness = Q_(0.01, 'mm')
    return pipe_schedule, pipe_roughness


def add_fittings_D2D3(network: PipeNetwork) -> None:
    """Adds the fittings to the pipe with ID 'D2D3'."""
    tee_1 = fittings.Tee(
        flow_pattern='diverging',
        combined_pipe=Pipe.create(
            length=Q_(50, 'cm'),
            wall_roughness=network.wall_roughness,
            fluid=network.fluid,
            cross_section=Circular.create(
                nominal_diameter=Q_(50, 'mm'),
                schedule=network.schedule
            ),
            volume_flow_rate=Q_(1.873, 'L / s')
        ),
        branch_pipe=network.conduits['D2D3']
    )
    network.conduits['D2D3'].add_fitting(zeta=tee_1.zeta_b, ID='TEE-D2')

    elbow = fittings.PipeFitting(
        pipe=network.conduits['D2D3'],
        ELR=30.0
    )
    network.conduits['D2D3'].add_fitting(
        elbow.zeta,
        ID=f"ELB90-{network.conduits['D2D3'].ID}"
    )

    tee_2 = fittings.Tee(
        flow_pattern='diverging',
        combined_pipe=network.conduits['D2D3'],
        branch_pipe=network.conduits['D3R3']
    )
    network.conduits['D2D3'].add_fitting(zeta=tee_2.zeta_c, ID='TEE-D3')


def add_fittings_D3D4(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='diverging',
        combined_pipe=network.conduits['D3D4'],
        branch_pipe=network.conduits['D4R4']
    )
    network.conduits['D3D4'].add_fitting(zeta=tee.zeta_c, ID='TEE-D4')


def add_fittings_D4D5(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='diverging',
        combined_pipe=network.conduits['D4D5'],
        branch_pipe=network.conduits['D5R5']
    )
    network.conduits['D4D5'].add_fitting(zeta=tee.zeta_c, ID='TEE-D5')


def add_fittings_D5D6(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='diverging',
        combined_pipe=network.conduits['D5D6'],
        branch_pipe=network.conduits['D6R6']
    )
    network.conduits['D5D6'].add_fitting(zeta=tee.zeta_c, ID='TEE-D6')


def add_fittings_D6D7(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='diverging',
        combined_pipe=network.conduits['D6D7'],
        branch_pipe=network.conduits['D7R7']
    )
    network.conduits['D6D7'].add_fitting(zeta=tee.zeta_c, ID='TEE-D7')


def add_fittings_D2D9(network: PipeNetwork) -> None:
    tee_1 = fittings.Tee(
        flow_pattern='diverging',
        combined_pipe=Pipe.create(
            length=Q_(50, 'cm'),
            wall_roughness=network.wall_roughness,
            fluid=network.fluid,
            cross_section=Circular.create(
                nominal_diameter=Q_(50, 'mm'),
                schedule=network.schedule
            ),
            volume_flow_rate=Q_(1.873, 'L / s')
        ),
        branch_pipe=network.conduits['D2D9']
    )
    network.conduits['D2D9'].add_fitting(zeta=tee_1.zeta_b, ID='TEE-D2')

    for i in range(2):
        elbow = fittings.PipeFitting(
            pipe=network.conduits['D2D9'],
            ELR=30.0
        )
        network.conduits['D2D9'].add_fitting(
            elbow.zeta,
            ID=f"ELB90-{network.conduits['D2D9'].ID}-{i+1}"
        )

    tee_2 = fittings.Tee(
        flow_pattern='diverging',
        combined_pipe=network.conduits['D2D9'],
        branch_pipe=network.conduits['D9R9']
    )
    network.conduits['D2D9'].add_fitting(zeta=tee_2.zeta_c, ID='TEE-D9')


def add_fittings_D9D10(network: PipeNetwork) -> None:
    for i in range(3):
        elbow = fittings.PipeFitting(
            pipe=network.conduits['D9D10'],
            ELR=30.0
        )
        network.conduits['D9D10'].add_fitting(
            elbow.zeta,
            ID=f"ELB90-{network.conduits['D9D10'].ID}-{i+1}"
        )


def add_fittings_R7R6(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='converging',
        combined_pipe=network.conduits['R7R6'],
        branch_pipe=network.conduits['D7R7']
    )
    network.conduits['R7R6'].add_fitting(zeta=tee.zeta_c, ID='TEE-R7')


def add_fittings_R6R5(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='converging',
        combined_pipe=network.conduits['R6R5'],
        branch_pipe=network.conduits['D6R6']
    )
    network.conduits['R6R5'].add_fitting(zeta=tee.zeta_c, ID='TEE-R6')


def add_fittings_R5R4(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='converging',
        combined_pipe=network.conduits['R5R4'],
        branch_pipe=network.conduits['D5R5']
    )
    network.conduits['R5R4'].add_fitting(zeta=tee.zeta_c, ID='TEE-R5')


def add_fittings_R4R3(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='converging',
        combined_pipe=network.conduits['R4R3'],
        branch_pipe=network.conduits['D4R4']
    )
    network.conduits['R4R3'].add_fitting(zeta=tee.zeta_c, ID='TEE-R4')


def add_fittings_R3R2(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='converging',
        combined_pipe=network.conduits['R3R2'],
        branch_pipe=network.conduits['D3R3']
    )
    network.conduits['R3R2'].add_fitting(zeta=tee.zeta_c, ID='TEE-R3')

    elbow = fittings.PipeFitting(
        pipe=network.conduits['R3R2'],
        ELR=30.0
    )
    network.conduits['R3R2'].add_fitting(
        elbow.zeta,
        ID=f"ELB90-{network.conduits['R3R2'].ID}"
    )


def add_fittings_R10R9(network: PipeNetwork) -> None:
    for i in range(3):
        elbow = fittings.PipeFitting(
            pipe=network.conduits['R10R9'],
            ELR=30.0
        )
        network.conduits['R10R9'].add_fitting(
            elbow.zeta,
            ID=f"ELB90-{network.conduits['R10R9'].ID}-{i+1}"
        )


def add_fittings_R9R2(network: PipeNetwork) -> None:
    tee = fittings.Tee(
        flow_pattern='converging',
        combined_pipe=network.conduits['R9R2'],
        branch_pipe=network.conduits['D9R9']
    )
    network.conduits['R9R2'].add_fitting(zeta=tee.zeta_c, ID='TEE-R9')

    for i in range(2):
        elbow = fittings.PipeFitting(
            pipe=network.conduits['R9R2'],
            ELR=30.0
        )
        network.conduits['R9R2'].add_fitting(
            elbow.zeta,
            ID=f"ELB90-{network.conduits['R9R2'].ID}-{i+1}"
        )


def add_fittings_to_cross_over(
    network: PipeNetwork,
    pipe: Pipe,
    supply_pipe: Pipe,
    return_pipe: Pipe,
    dP_fcu: Quantity,
    num_elbows: int = 2,
    has_tees: bool = True
) -> None:
    """General function to place the fittings in the cross-over pipes."""
    if has_tees:
        tee_1 = fittings.Tee(
            flow_pattern='diverging',
            combined_pipe=supply_pipe,
            branch_pipe=pipe
        )
        pipe.add_fitting(
            zeta=tee_1.zeta_b,
            ID=f"TEE-{pipe.start_node.ID}"
        )

    reducer = fittings.Reducer(
        large_pipe=supply_pipe,
        small_pipe=pipe,
        L=Q_(1, 'cm')
    )
    pipe.add_fitting(
        zeta=reducer.zeta_small,
        ID=f"RED-{pipe.ID}"
    )

    for i in range(num_elbows):
        elbow = fittings.PipeFitting(
            pipe=pipe,
            ELR=30.0
        )
        pipe.add_fitting(
            zeta=elbow.zeta,
            ID=f"ELB90-{pipe.ID}-{i+1}"
        )

    ball_valve = fittings.PipeFitting(
        pipe=pipe,
        ELR=3.0
    )
    pipe.add_fitting(
        zeta=ball_valve.zeta,
        ID=f"BALL-VALVE-{pipe.ID}"
    )

    Kv_fcu = fittings.FlowCoefficient.get_Kv(
        volume_flow_rate=pipe.volume_flow_rate,
        pressure_drop=dP_fcu,
        fluid=network.fluid
    )
    zeta_fcu = fittings.ResistanceCoefficient.from_Kv(
        Kv=Kv_fcu,
        diameter=pipe.cross_section.equivalent_diameter
    )
    pipe.add_fitting(
        zeta=zeta_fcu,
        ID=f"FCU-{pipe.ID}"
    )

    enlarger = fittings.Enlarger(
        large_pipe=return_pipe,
        small_pipe=pipe,
        L=Q_(1, 'cm')
    )
    pipe.add_fitting(
        zeta=enlarger.zeta_small,
        ID=f"ENL-{pipe.ID}"
    )

    if has_tees:
        tee_2 = fittings.Tee(
            flow_pattern='converging',
            combined_pipe=return_pipe,
            branch_pipe=pipe
        )
        pipe.add_fitting(
            zeta=tee_2.zeta_b,
            ID=f"TEE-{pipe.end_node.ID}"
        )


def add_valves(network: PipeNetwork) -> None:
    """Adds the balancing valves and the control valves to the cross-overs."""
    cross_over_lst = [
        'D8R8', 'D7R7', 'D6R6', 'D5R5',
        'D4R4', 'D3R3', 'D9R9', 'D10R10'
    ]
    for ID in cross_over_lst:
        network.add_balancing_valve(ID)
        network.set_balancing_valve_Kvs(ID, 4.46)
        # balancing valve DN20, based on Caleffi Serie 130 DN20 (art. code 130500)

        network.add_control_valve(ID)
        network.set_control_valve_Kvs(ID, 6.3)
        # control valve DN20, based on Caleffi Serie 6452 DN20 (1/2")

    # Statically balance the pipe network for the design volume flow rates
    # through the cross-overs; this sets the valve openings of the balancing
    # valves, while the control valves remain fully open:
    network.balance_network_at_design()


def create_curves(network: PipeNetwork) -> LineChart:
    """Creates the system curve of the pipe network and the pump curve of the
    circulation pump.
    """
    system_curve = SystemCurve.create(
        R_hyd=network.hydraulic_resistance,
        name='system curve'
    )
    pump_curve = PumpCurve.create(
        coeffs=network.conduits['R2D2'].machine_coefficients,
        # the polynomial coefficients of the second-order polynomial that models
        # the pump curve can be retrieved from the `Pipe` object containing the
        # pump as indicated in the network configuration file.
        name='pump curve'
    )
    diagram = plot_curves(
        pump_curves=[pump_curve],
        system_curves=[system_curve],
        working_point=(
            network.volume_flow_rate,
            network.total_pressure_difference
        ),
        fig_size=(8, 6),
        V_unit='L / s',
        dP_unit='bar',
        V_max=Q_(6.5, 'L / s'),
        dP_max=Q_(1.6, 'bar'),
        V_step=Q_(0.5, 'L / s'),
        dP_step=Q_(0.2, 'bar')
    )
    diagram.add_legend(
        anchor='upper right',
        position=(0.99, 0.99),
        columns=1
    )
    return diagram


if __name__ == '__main__':
    main()
