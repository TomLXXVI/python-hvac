from typing import Optional, Tuple, List
import numpy as np
from hvac import Quantity
from hvac.charts import LineChart
from .pump_curve import PumpCurve

Q_ = Quantity


class SystemCurve:

    def __init__(self):
        self.name: str = ''
        self._R_hyd: Quantity = Q_(float('nan'), 'Pa * s ** 2 / m ** 6')
        self._dP_tot: Quantity = Q_(float('nan'), 'Pa')
        self._dP_elev: Quantity = Q_(float('nan'), 'Pa')

    @classmethod
    def create(
            cls,
            R_hyd: Quantity,
            dP_tot: Optional[Quantity] = None,
            dP_elev: Optional[Quantity] = None,
            name: str = ''
    ) -> 'SystemCurve':
        obj = cls()
        obj._R_hyd = R_hyd.to('Pa * s ** 2 / m ** 6').magnitude
        obj._dP_tot = dP_tot.to('Pa').magnitude if dP_tot is not None else 0.0
        obj._dP_elev = dP_elev.to('Pa').magnitude if dP_elev is not None else 0.0
        obj.name = name
        return obj

    def axes(self, V_ini: Quantity, V_fin: Quantity, num: int = 50) -> Tuple[Quantity, Quantity]:
        """
        Get the coordinate-axes of the system curve.

        Parameters
        ----------
        V_ini: Quantity
            The start point of the system curve.
        V_fin : Quantity
            The end point of the system curve.
        num: int, optional
            The number of coordinates on the system curve to be calculated, start- and endpoint included.
            Default is 50.

        Returns
        -------
        Tuple[Quantity, Quantity]
            The volume flow rate axis and pressure difference axis of the system curve.
        """
        V_ini = V_ini.to('m ** 3 / s').magnitude
        V_fin = V_fin.to('m ** 3 / s').magnitude
        V_ax = np.linspace(V_ini, V_fin, num, endpoint=True)
        dP_ax = self._dP_tot + self._dP_elev + self._R_hyd * V_ax ** 2
        dP_ax = Quantity(dP_ax, 'Pa')
        V_ax = Quantity(V_ax, 'm ** 3 / s')
        return V_ax, dP_ax

    def pressure_loss(self, V: Quantity):
        V = V.to('m ** 3 / s').magnitude
        dP = self._dP_tot + self._dP_elev + self._R_hyd * V ** 2
        return Quantity(dP, 'Pa')

    def diagram(
            self,
            V_ini: Quantity,
            V_fin: Quantity,
            V_design: Optional[Quantity] = None,
            no_points: int = 50,
            **fig_kwargs
    ) -> LineChart:
        fig_size: Tuple[int, int] = fig_kwargs.get('fig_size', (6, 4))
        dpi: int = fig_kwargs.get('dpi', 96)
        V_step: Quantity = fig_kwargs.get('V_step')
        V_max: Quantity = fig_kwargs.get('V_max')
        dP_step: Quantity = fig_kwargs.get('dP_step')
        dP_max: Quantity = fig_kwargs.get('dP_max')
        V_unit: str = fig_kwargs.get('V_unit', 'm ** 3 / s')
        dP_unit: str = fig_kwargs.get('dP_unit', 'Pa')

        V_ax, dP_ax = self.axes(V_ini, V_fin, no_points)
        V_ax = V_ax.to(V_unit).magnitude
        dP_ax = dP_ax.to(dP_unit).magnitude

        g = LineChart(size=fig_size, dpi=dpi)
        g.add_xy_data(
            label=self.name if self.name else 'system curve',
            x1_values=V_ax,
            y1_values=dP_ax
        )
        if V_design is not None:
            dP_des = self.pressure_loss(V_design)
            g.add_xy_data(
                label='working point',
                x1_values=V_design.to(V_unit).magnitude,
                y1_values=dP_des.to(dP_unit).magnitude,
                style_props={'marker': 'o', 'linestyle': 'None', 'color': 'red'}
            )
        if V_max is not None and V_step is not None:
            g.x1.scale(lower_limit=0.0, upper_limit=V_max.to(V_unit).magnitude, step=V_step.to(V_unit).magnitude)
        if dP_max is not None and dP_step is not None:
            g.y1.scale(lower_limit=0.0, upper_limit=dP_max.to(dP_unit).magnitude, step=dP_step.to(dP_unit).magnitude)
        g.x1.add_title(f'volume flow rate, {V_unit}')
        g.y1.add_title(f'pressure loss, {dP_unit}')
        return g


def plot_curves(
        pump_curves: List[PumpCurve],
        system_curves: List[SystemCurve],
        working_point: Optional[Tuple[Quantity, Quantity]] = None,
        **fig_kwargs
) -> LineChart:
    fig_size: Tuple[int, int] = fig_kwargs.get('fig_size', (6, 4))
    dpi: int = fig_kwargs.get('dpi', 96)

    V_step: Quantity = fig_kwargs.get('V_step')
    V_max: Quantity = fig_kwargs.get('V_max', Quantity(10, 'm ** 3 / hr'))
    dP_step: Quantity = fig_kwargs.get('dP_step')
    dP_max: Quantity = fig_kwargs.get('dP_max')
    V_unit: str = fig_kwargs.get('V_unit', 'm ** 3 / s')
    dP_unit: str = fig_kwargs.get('dP_unit', 'Pa')

    g = LineChart(size=fig_size, dpi=dpi, constructs=fig_kwargs.get('figure_constructs'))
    for i, pump_curve in enumerate(pump_curves):
        V_ax, dP_ax = pump_curve.axes()
        g.add_xy_data(
            label=pump_curve.name if pump_curve.name else f"pump curve {i+1}",
            x1_values=V_ax.to(V_unit).magnitude,
            y1_values=dP_ax.to(dP_unit).magnitude
        )
    for i, system_curve in enumerate(system_curves):
        V_ax, dP_ax = system_curve.axes(V_ini=Quantity(0.0, 'm ** 3 / s'), V_fin=V_max)
        g.add_xy_data(
            label=system_curve.name if system_curve.name else f"system curve {i}",
            x1_values=V_ax.to(V_unit).magnitude,
            y1_values=dP_ax.to(dP_unit).magnitude
        )
    if working_point:
        g.add_xy_data(
            label='working point',
            x1_values=working_point[0].to(V_unit).magnitude,
            y1_values=working_point[1].to(dP_unit).magnitude,
            style_props={'marker': 'o', 'color': 'red', 'linestyle': 'None'}
        )
    if V_max is not None and V_step is not None:
        g.x1.scale(lower_limit=0.0, upper_limit=V_max.to(V_unit).magnitude, step=V_step.to(V_unit).magnitude)
    if dP_max is not None and dP_step is not None:
        g.y1.scale(lower_limit=0.0, upper_limit=dP_max.to(dP_unit).magnitude, step=dP_step.to(dP_unit).magnitude)
    g.x1.add_title(f'V, {V_unit}')
    g.y1.add_title(f'dP, {dP_unit}')
    g.add_legend(anchor='upper center', position=(0.5, -0.1), columns=fig_kwargs.get('column_count', 2))
    return g
