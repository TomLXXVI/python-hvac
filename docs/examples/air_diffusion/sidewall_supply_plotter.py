import numpy as np
from collections import namedtuple
from hvac import Quantity
from hvac.charts import LineChart
from hvac.air_diffusion import SideWallSupply, RoomInfo
from hvac.air_diffusion.air_jets.archimedes_number import room_archimedes_number

Q_ = Quantity

ConfinedJet = namedtuple(
    'ConfinedJet',
    ('x', 'y', 'u', 'T', 'x_oz', 'y_oz', 'u_oz', 'T_oz')
)

FreeJet = namedtuple(
    'FreeJet',
    ('x', 'y', 'u', 'T', 'x_oz', 'y_oz', 'u_oz', 'T_oz')
)


class JetPlotter:

    def __init__(self, sws: SideWallSupply, room: RoomInfo, L: Quantity):
        self.sws = sws
        self.room = room
        self.x = Q_(np.linspace(0.1, L.to('m').m, 500, endpoint=True), 'm')
        self.fj = self.free_jet()
        self.cj = self.confined_jet()

    @staticmethod
    def find_target(target: Quantity, values: Quantity) -> int:
        diff_arr = np.absolute(values - target)
        index = diff_arr.argmin()
        return index

    def confined_jet(self) -> ConfinedJet:
        # Confined jet: trajectory, centerline velocity and centerline temperature.
        y_cj, u_cj, T_cj = self.sws.trajectory(self.x)
        # Get the index of the y-value nearest to the target value (the y-value that
        # corresponds with the height of the occupied zone):
        i = self.find_target(target=-self.room.H + self.room.Z, values=y_cj)
        cj = ConfinedJet(
            x=self.x, y=y_cj, u=u_cj, T=T_cj,
            x_oz=self.x[i], y_oz=y_cj[i], u_oz=u_cj[i], T_oz=T_cj[i]
        )
        return cj

    def free_jet(self) -> FreeJet:
        # Free jet: trajectory, centerline velocity and centerline temperature.
        y_fj = self.sws.compact_jet.trajectory(self.x)
        u_fj = self.sws.compact_jet.centerline_velocity_zone3(self.x)
        T_fj = self.sws.compact_jet.centerline_temperature_zone3(self.x)
        # Get the index of the y-value nearest to the target value (the y-value that
        # corresponds with the height of the occupied zone):
        i = self.find_target(target=-self.room.H + self.room.Z, values=y_fj)
        fj = FreeJet(
            x=self.x, y=y_fj, u=u_fj, T=T_fj,
            x_oz=self.x[i], y_oz=y_fj[i], u_oz=u_fj[i], T_oz=T_fj[i]
        )
        return fj

    def plot_jet_trajectory(self) -> LineChart:
        # Draw trajectories of confined and free jet:
        chart1 = LineChart(size=(8, 6))
        chart1.add_xy_data(
            label='confined jet',
            x1_values=self.cj.x,
            y1_values=self.cj.y,
            style_props={'color': 'tab:blue'}
        )
        chart1.add_xy_data(
            label=f"{self.cj.u_oz.to('m / s'):~P.2f}, {self.cj.T_oz.to('degC'):~P.1f}",
            x1_values=[self.cj.x_oz.to('m').m],
            y1_values=[self.cj.y_oz.to('m').m],
            style_props={'marker': 'o', 'color': 'tab:blue'}
        )
        # chart1.add_note(
        #     text=f"({self.cj.u_oz.to('m / s'):~P.2f}, {self.cj.T_oz.to('degC'):~P.1f})",
        #     x_pos=self.cj.x_oz.to('m').m + 0.3,
        #     y_pos=self.cj.y_oz.to('m').m,
        #     font_size=10,
        #     vert_align='bottom',
        #     hor_align='left',
        #     use_normalized_coordinates=False
        # )
        chart1.add_xy_data(
            label='free jet',
            x1_values=self.fj.x,
            y1_values=self.fj.y,
            style_props={'color': 'tab:orange'}
        )
        chart1.add_xy_data(
            label=f"{self.fj.u_oz.to('m / s'):~P.2f}, {self.fj.T_oz.to('degC'):~P.1f}",
            x1_values=[self.fj.x_oz.to('m').m],
            y1_values=[self.fj.y_oz.to('m').m],
            style_props={'marker': 'o', 'color': 'tab:orange'}
        )
        # chart1.add_note(
        #     text=f"({self.fj.u_oz.to('m / s'):~P.2f}, {self.fj.T_oz.to('degC'):~P.1f})",
        #     x_pos=self.fj.x_oz.to('m').m + 0.3,
        #     y_pos=self.fj.y_oz.to('m').m,
        #     font_size=10,
        #     vert_align='bottom',
        #     hor_align='left',
        #     use_normalized_coordinates=False
        # )
        chart1.add_xy_data(
            label='room length / drop height',
            x1_values=[0.0, self.room.L.m, self.room.L.m],
            y1_values=[-self.sws.output.dy_max.m, -self.sws.output.dy_max.m, -self.room.H.m],
            style_props={'linestyle': '--', 'color': 'red'}
        )
        chart1.x1.add_title('horizontal distance, m')
        chart1.y1.add_title('vertical distance, m')
        chart1.x1.scale(0, 24, 1)
        chart1.y1.scale(-6.0, 1.5, 0.5)
        chart1.add_legend()
        return chart1

    def plot_jet_velocity(self) -> LineChart:
        # Draw centerline velocity of confined and free jet:
        chart2 = LineChart(size=(8, 6))
        chart2.add_xy_data(
            label='confined jet',
            x1_values=self.cj.x,
            y1_values=self.cj.u,  # self.room.H + self.cj.y,
            style_props={'color': 'tab:blue'}
        )
        chart2.add_xy_data(
            label='free jet',
            x1_values=self.fj.x,
            y1_values=self.fj.u,  # self.room.H + self.fj.y,
            style_props={'linestyle': '--', 'color': 'tab:orange'}
        )
        # chart2.add_xy_data(
        #     label='0.25 m/s',
        #     x1_values=[0.25, 0.25],
        #     y1_values=[0, 6.5],
        #     style_props={'color': 'red', 'linestyle': ':'}
        # )
        chart2.x1.add_title('horizontal distance, m')
        chart2.y1.add_title('centerline velocity, m/s')
        # chart2.x1.add_title('centerline velocity, m/s')
        # chart2.y1.add_title('vertical distance, m')
        # chart2.x1.scale(0, round(self.sws.output.U_o.to('m / s').m, 1), 0.05)
        # chart2.y1.scale(0, 7.0, 0.5)
        return chart2

    def plot_jet_temperature(self) -> LineChart:
        # Draw centerline temperature of confined and free jet as function of
        # height:
        chart3 = LineChart(size=(8, 6))
        chart3.add_xy_data(
            label='confined jet',
            x1_values=self.cj.x,
            y1_values=self.cj.T.to('degC'),  # self.room.H + self.cj.y,
            style_props={'color': 'tab:blue'}
        )
        chart3.add_xy_data(
            label='free jet',
            x1_values=self.fj.x,
            y1_values=self.fj.T.to('degC'),  # self.room.H + self.fj.y,
            style_props={'linestyle': '--', 'color': 'tab:orange'}
        )
        chart3.x1.add_title('horizontal distance, m')
        chart3.y1.add_title('centerline temperature, °C')
        # chart3.x1.add_title('centerline temperature, °C')
        # chart3.y1.add_title('vertical distance, m')
        # chart3.x1.scale(
        #     round(self.sws.output.T_o.to('degC').m, 0),
        #     round(self.room.T_r.to('degC').m) + 0.5,
        #     0.5
        # )
        # chart3.y1.scale(0, 7.0, 0.5)
        return chart3
