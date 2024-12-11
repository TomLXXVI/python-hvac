"""
DESIGN OF SIDE-WALL AIR SUPPLY

For a room with given dimensions, room air temperature, and room load, determine
the required supply air volume flow rate, supply air temperature, and size and
aspect ratio of the supply opening, such that the critical room Archimedes
number and the allowable mean room air speed are met.

This script also compares the trajectory of the confined jet to the trajectory
under the assumption of a free jet.
"""
import numpy as np
from hvac import Quantity
from hvac.charts import LineChart
from hvac.air_diffusion.design.side_wall_supply import RoomInfo, SideWallSupply
from hvac.air_diffusion.air_jets.archimedes_number import room_archimedes_number

Q_ = Quantity


def find_target(target: Quantity, values: Quantity) -> int:
    diff_arr = np.absolute(values - target)
    index = diff_arr.argmin()
    return index


def main():
    room = RoomInfo(
        L=Q_(23 / 4, 'm'),
        B=Q_(11, 'm'),
        H=Q_(6, 'm'),
        Z=Q_(1.8, 'm'),
        T_r=Q_(26, 'degC'),
        Q_dot=Q_(25.452 / 4, 'kW')
    )
    print(f"room load per unit area: {room.q_dot.to('W / m**2'):~P.3f}")

    # Throw constant of supply outlet:
    K1 = 0.5
    # Height between upper edge of supply opening and ceiling:
    d = Q_(0.5, 'm')

    sws = SideWallSupply(room, d, K1, Ar_r_crit=11_000)
    # The critical room Archimedes number depends on the ratio room length to
    # room height: see e.g. Awbi, H. B. (2003). Ventilation of Buildings.
    # Taylor & Francis. p. 246, table 6.2.

    output = sws.design(
        h_o=Q_(0.2, 'm'),
        L_th_frac=0.25,
        v_r_minmax=(Q_(0.1, 'm / s'), Q_(0.3, 'm / s'))
    )
    print(output)

    # Check room Archimedes number:
    Ar_r = room_archimedes_number(
        B=room.B,
        H=room.H,
        T_r=room.T_r,
        V_dot=output.V_dot,
        dT_o=output.dT_o
    )
    print(round(Ar_r))
    
    # Confined jet:
    x = Q_(np.linspace(0.1, 23.0, 500, endpoint=True), 'm')
    y_cj, U_x_cj, T_x_cj = sws.trajectory(x)

    # Get the index of the y-value nearest to the target value (the y-value that
    # corresponds with the height of the occupied zone):
    i = find_target(target=-room.H + room.Z, values=y_cj)
    x_i_cj, y_i_cj, U_x_i_cj, T_x_i_cj = x[i], y_cj[i], U_x_cj[i], T_x_cj[i]
    
    # Free jet:
    y_fj = sws.compact_jet.trajectory(x)
    U_x_fj = sws.compact_jet.centerline_velocity_zone3(x)
    T_x_fj = sws.compact_jet.centerline_temperature_zone3(x)
    i = find_target(target=-room.H + room.Z, values=y_fj)
    x_i_fj, y_i_fj, U_x_i_fj, T_x_i_fj = x[i], y_fj[i], U_x_fj[i], T_x_fj[i]

    # Draw trajectories of confined and free jet:
    chart1 = LineChart(size=(8, 6))

    chart1.add_xy_data(
        label='confined jet',
        x1_values=x,
        y1_values=y_cj
    )
    chart1.add_xy_data(
        label='occupied zone (cj)',
        x1_values=[x_i_cj.to('m').m],
        y1_values=[y_i_cj.to('m').m],
        style_props={'marker': 'o', 'color': 'red'}
    )
    chart1.add_note(
        text=f"({U_x_i_cj.to('m / s'):~P.2f}, {T_x_i_cj.to('degC'):~P.1f})",
        x_pos=x_i_cj.to('m').m + 0.3,
        y_pos=y_i_cj.to('m').m,
        font_size=10,
        vert_align='bottom',
        hor_align='left',
        use_normalized_coordinates=False
    )

    chart1.add_xy_data(
        label='free jet',
        x1_values=x,
        y1_values=y_fj
    )
    chart1.add_xy_data(
        label='occupied zone (fj)',
        x1_values=[x_i_fj.to('m').m],
        y1_values=[y_i_fj.to('m').m],
        style_props={'marker': 'o', 'color': 'red'}
    )
    chart1.add_note(
        text=f"({U_x_i_fj.to('m / s'):~P.2f}, {T_x_i_fj.to('degC'):~P.1f})",
        x_pos=x_i_fj.to('m').m - 0.3,
        y_pos=y_i_fj.to('m').m,
        font_size=10,
        vert_align='top',
        hor_align='right',
        use_normalized_coordinates=False
    )

    chart1.add_xy_data(
        label='room length / drop height',
        x1_values=[0.0, room.L.m, room.L.m],
        y1_values=[-output.dy_max.m, -output.dy_max.m, -room.H.m],
        style_props={'linestyle': '--', 'color': 'red'}
    )

    chart1.x1.add_title('horizontal distance, m')
    chart1.y1.add_title('vertical distance, m')
    chart1.add_legend()
    chart1.x1.scale(0.0, 24.0, 1.0)
    chart1.y1.scale(-6.0, 1.0, 0.5)
    chart1.show()

    # Draw centerline velocity of confined and free jet:
    chart2 = LineChart(size=(8, 6))
    chart2.add_xy_data(
        label='confined jet',
        x1_values=U_x_cj,
        y1_values=room.H + y_cj
    )
    chart2.add_xy_data(
        label='free jet',
        x1_values=U_x_fj,
        y1_values=room.H + y_fj,
        style_props={'linestyle': '--'}
    )
    chart2.add_xy_data(
        label='0.25 m/s',
        x1_values=[0.25, 0.25],
        y1_values=[0, 6.5],
        style_props={'color': 'red', 'linestyle': ':'}
    )
    chart2.x1.add_title('centerline velocity, m/s')
    chart2.y1.add_title('vertical distance, m')
    chart2.x1.scale(0, 2.1, 0.1)
    chart2.y1.scale(0, 7.0, 0.5)
    chart2.show()

    # Draw centerline temperature of confined and free jet:
    chart3 = LineChart(size=(8, 6))
    chart3.add_xy_data(
        label='confined jet',
        x1_values=T_x_cj.to('degC'),
        y1_values=room.H + y_cj
    )
    chart3.add_xy_data(
        label='free jet',
        x1_values=T_x_fj.to('degC'),
        y1_values=room.H + y_fj,
        style_props={'linestyle': '--'}
    )
    chart3.x1.add_title('centerline temperature, °C')
    chart3.y1.add_title('vertical distance, m')
    chart3.x1.scale(round(output.T_o.to('degC').m), round(room.T_r.to('degC').m) + 0.5, 0.5)
    chart3.y1.scale(0, 7.0, 0.5)
    chart3.show()


if __name__ == '__main__':
    main()
