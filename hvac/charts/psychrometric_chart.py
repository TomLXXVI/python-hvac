from typing import Optional, List, Tuple
from dataclasses import dataclass
import os
import pickle
from CoolProp.CoolProp import HAPropsSI
import numpy as np
from matplotlib.patches import ConnectionPatch
from matplotlib.ticker import MultipleLocator
from .. import Quantity as Qty
from ..fluids import HumidAir, STANDARD_PRESSURE
from .matplotlibwrapper import LineChart


@dataclass
class StatePoint:
    Tdb: Qty
    W: Qty


class PsychrometricChart:
    CHARTS_PATH = os.path.abspath('./charts')
    P = STANDARD_PRESSURE

    def __init__(self, fig_size: Optional[Tuple[int, int]] = None):
        self.chart: Optional[LineChart] = None
        try:
            self._load(fig_size)
        except FileNotFoundError:
            self._create(fig_size if fig_size is not None else (12, 8), 96)
            self._dump()

    def _create(self, fig_size, dpi):
        self.chart = LineChart(size=fig_size, dpi=dpi)
        T_db_vec = np.linspace(-30, 55) + 273.15

        # lines of constant relative humidity
        for RH in np.arange(0.1, 1, 0.1):
            W = HAPropsSI("W", "R", RH, "P", self.P.to('Pa').magnitude, "T", T_db_vec)
            self.chart.axes.plot(T_db_vec - 273.15, W, color='k', lw=0.5)

        # saturation curve
        W_sat = HAPropsSI("W", "R", 1, "P", self.P.to('Pa').magnitude, "T", T_db_vec)
        self.chart.axes.plot(T_db_vec - 273.15, W_sat, color='k', lw=1.5)

        # lines of constant v_da
        for v_da in np.arange(0.69, 0.961, 0.01):
            R = np.linspace(0, 1)
            W = HAPropsSI("W", "R", R, "P", self.P.to('Pa').magnitude, "Vda", v_da)
            T_db = HAPropsSI("Tdb", "R", R, "P", self.P.to('Pa').magnitude, "Vda", v_da)
            self.chart.axes.plot(T_db - 273.15, W, color='b', lw=1.5 if abs(v_da % 0.05) < 0.001 else 0.5)

        # lines of constant T_wb
        for T_wb_C in np.arange(-16, 33, 2):
            if T_wb_C == 0:
                continue
            R = np.linspace(0.0, 1)
            T_db = HAPropsSI("Tdb", "R", R, "P", self.P.to('Pa').magnitude, "Twb", T_wb_C + 273.15)
            W = HAPropsSI("W", "R", R, "P", self.P.to('Pa').magnitude, "Tdb", T_db)
            self.chart.axes.plot(T_db - 273.15, W, color='r', lw=1.5 if abs(T_wb_C % 10) < 0.001 else 0.5)

        self.chart.x1.add_title('dry bulb temperature, °C')
        self.chart.y1.add_title('humidity ratio, kg_H2O/kg_da')
        self.chart.y1.scale(lower_limit=0, upper_limit=0.030, step=0.005)
        self.chart.x1.scale(lower_limit=-30, upper_limit=55, step=5)
        self.chart.x1.axes.xaxis.set_minor_locator(MultipleLocator())
        self.chart.x1.axes.yaxis.set_minor_locator(MultipleLocator(0.001))

    def plot_process(
        self,
        name: str,
        start_point: StatePoint | HumidAir,
        end_point: StatePoint | HumidAir,
        mix_point: Optional[StatePoint | HumidAir] = None,
    ):
        if mix_point is not None:
            x_data = [start_point.Tdb.to('degC').m, mix_point.Tdb.to('degC').m, end_point.Tdb.to('degC').m]
            y_data = [start_point.W.to('kg/kg').m, mix_point.W.to('kg/kg').m, end_point.W.to('kg/kg').m]
        else:
            x_data = [start_point.Tdb.to('degC').m, end_point.Tdb.to('degC').m]
            y_data = [start_point.W.to('kg/kg').m, end_point.W.to('kg/kg').m]

        self.chart.add_xy_data(
            label=name,
            x1_values=x_data,
            y1_values=y_data,
            style_props={'marker': 'o', 'color': 'orange'}
        )

        if mix_point is None:
            # add arrow to process line to show the direction of the process
            xyA = (start_point.Tdb.to('degC').m, start_point.W.to('kg/kg').m)
            xyB = (end_point.Tdb.to('degC').m, end_point.W.to('kg/kg').m)
            coordsA = "data"
            coordsB = "data"
            con = ConnectionPatch(
                xyA, xyB, coordsA, coordsB,
                arrowstyle="-|>", shrinkA=5, shrinkB=5, mutation_scale=20, color='orange', lw=2.0
            )
            self.chart.axes.add_artist(con)

    def plot_point(self, name: str, point: StatePoint | HumidAir):
        self.chart.add_xy_data(
            label=name,
            x1_values=[point.Tdb.to('degC').m],
            y1_values=[point.W.to('kg/kg').m],
            style_props={'marker': 'o', 'color': 'orange', 'lw': 2.0}
        )

    def plot_line(
        self,
        name: str,
        start_point: StatePoint | HumidAir,
        end_point: StatePoint | HumidAir
    ):
        x_data = [start_point.Tdb.to('degC').m, end_point.Tdb.to('degC').m]
        y_data = [start_point.W.to('kg/kg').m, end_point.W.to('kg/kg').m]

        self.chart.add_xy_data(
            label=name,
            x1_values=x_data,
            y1_values=y_data,
            style_props={'color': 'orange'}
        )

    def plot_space_condition_line(
        self,
        start_point: StatePoint | HumidAir,
        end_point: StatePoint | HumidAir,
        space_point: Optional[StatePoint | HumidAir] = None,
    ):
        self.plot_line('space condition line', start_point, end_point)
        self.chart.add_xy_data(
            label='space point',
            x1_values=[space_point.Tdb.to('degC').m],
            y1_values=[space_point.W.to('kg/kg').m],
            style_props={'marker': 'o', 'color': 'orange'}
        )

    def plot_curve(self, name: str, state_points: List[StatePoint | HumidAir]):
        x_data = [p.Tdb.to('degC').m for p in state_points]
        y_data = [p.W.to('kg/kg').m for p in state_points]

        self.chart.add_xy_data(
            label=name,
            x1_values=x_data,
            y1_values=y_data,
            style_props={'color': 'orange'}
        )

    def show(self):
        self.chart.show()

    def _dump(self):
        fp = os.path.join(self.CHARTS_PATH, 'psych_chart.pickle')
        pickle.dump(self.chart.figure, open(fp, 'wb'))

    def _load(self, fig_size: Tuple[int, int]):
        if not os.path.exists(self.CHARTS_PATH):
            os.makedirs(self.CHARTS_PATH)
        fp = os.path.join(self.CHARTS_PATH, 'psych_chart.pickle')
        fig = pickle.load(open(fp, 'rb'))
        self.chart = LineChart(constructs=(fig, fig.axes[0]))
        if fig_size is not None:
            self.chart.figure.set_figwidth(fig_size[0])
            self.chart.figure.set_figheight(fig_size[1])
