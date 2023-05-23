from typing import List, Tuple, Optional
from hvac.pint_setup import UNITS
import numpy as np
from pint import Unit
from scipy.optimize import curve_fit
from hvac import Quantity
from hvac.charts import LineChart


class PumpCurve:

    def __init__(self):
        self.name: str = ''
        self._a0 = Quantity(0.0, 'Pa')
        self._a1 = Quantity(0.0, 'Pa * s / m ** 3')
        self._a2 = Quantity(0.0, 'Pa * s ** 2 / m ** 6')

    @classmethod
    def create(cls, *coefficients: Quantity, name: str = '') -> 'PumpCurve':
        obj = cls()
        obj.name = name
        obj._a0, obj._a1, obj._a2 = coefficients
        return obj

    @classmethod
    def curve_fitting(
        cls, 
        coordinates: List[Tuple[Quantity, Quantity]], 
        name: str = '', analysis_flow_units: Unit = UNITS('m ** 3 / s')
    ) -> 'PumpCurve':
            V_ax = []
            dP_ax = []
            for V, dP in coordinates:
                V_ax.append(V.to(analysis_flow_units).magnitude)
                dP_ax.append(dP.to('Pa').magnitude)

            def objective(x, a, b, c):
                return a * x + b * x ** 2 + c

            popt, _ = curve_fit(objective, V_ax, dP_ax)
            a, b, c = popt

            pump_curve = cls()
            pump_curve.name = name
            pump_curve._a0 = Quantity(c, 'Pa')
            pump_curve._a1 = Quantity(a, 'Pa') / analysis_flow_units
            pump_curve._a2 = Quantity(a, 'Pa') / (analysis_flow_units * analysis_flow_units)
            return pump_curve

    def pump_pressure(self, V: Quantity) -> Quantity:
        dP = self._a0 + self._a1 * V + self._a2 * V ** 2
        return dP

    def volume_flow_rate(self, dP: Quantity) -> Quantity:
        dP = dP.to('Pa').magnitude
        coeff = [self._a2, self._a1, self._a0 - dP]
        roots = np.roots(coeff)
        V = roots[roots > 0]
        V = V[0]
        return Quantity(V, 'm ** 3 / s')

    @property
    def coefficients(self) -> Tuple[Quantity, ...]:
        return self._a0, self._a1, self._a2

    @coefficients.setter
    def coefficients(self, values: Tuple[Quantity, ...]):
        self._a0, self._a1, self._a2 = values

    def axes(
            self,
            V_ini: Optional[Quantity] = None,
            V_fin: Optional[Quantity] = None,
            num: int = 50
    ) -> Tuple[Quantity, Quantity]:
        if V_ini is None:
            V_ini = 0.0
        if V_fin is None:
            V_fin = self.volume_flow_rate(dP=Quantity(0.0, 'Pa'))  # find V for which pump pressure is zero
        V_ax = np.linspace(V_ini, V_fin, num, endpoint=True)
        dP_ax = self._a0 + self._a1 * V_ax + self._a2 * V_ax ** 2
        return V_ax, dP_ax

    def diagram(
            self,
            V_ini: Quantity,
            V_fin: Quantity,
            num: int = 50,
            working_point: Optional[Tuple[Quantity, Quantity]] = None,
            **fig_kwargs
    ) -> LineChart:
        """
        Plot pump curve.

        Parameters
        ----------
        V_ini: Quantity
            First volume flow rate for which the pump curve is calculated.
        V_fin: Quantity
            Final volume flow rate for which the pump curve is calculated.
        num: int
            Number of volume flow rates between `V_ini` and `V_fin` for which the pump curve is calculated.
        working_point: Tuple[Quantity, Quantity]
            Volume flow rate and corresponding pump pressure at given working point. The working point will be
            shown as a red dot on the diagram.
        fig_kwargs:
            Figure settings.

        Following figure settings are available:

        - fig_size: Tuple[int, int]
            The width and height of the figure in inches.
        - dpi: int
            Dots per inch.
        - V_step: Quantity
            The step size between successive volume flow rates on the diagram axis.
        - V_max: Quantity
            The last volume flow rate indicated on the diagram axis.
        - dP_step: Quantity
            The step size between successive pressures on the diagram axis.
        - dP_max: Quantity
            The last pressure on the diagram axis.
        - V_unit: str
            The unit in which volume flow rates are expressed on the diagram.
        - dP_unit: str
            The unit in which pressures are expressed on the diagram.

        Returns
        -------
        LineChart
        """
        fig_size: Tuple[int, int] = fig_kwargs.get('fig_size', (6, 4))
        dpi: int = fig_kwargs.get('dpi', 96)
        V_step: Quantity = fig_kwargs.get('V_step')
        V_max: Quantity = fig_kwargs.get('V_max')
        dP_step: Quantity = fig_kwargs.get('dP_step')
        dP_max: Quantity = fig_kwargs.get('dP_max')
        V_unit: str = fig_kwargs.get('V_unit', 'm ** 3 / s')
        dP_unit: str = fig_kwargs.get('dP_unit', 'Pa')

        V_ax, dP_ax = self.axes(V_ini, V_fin, num)
        V_ax = V_ax.to(V_unit).magnitude
        dP_ax = dP_ax.to(dP_unit).magnitude
        diagram = LineChart(size=fig_size, dpi=dpi)
        diagram.add_xy_data(
            label=self.name if self.name else 'pump curve',
            x1_values=V_ax,
            y1_values=dP_ax
        )
        if working_point is not None:
            diagram.add_xy_data(
                label='working point',
                x1_values=working_point[0].to(V_unit).magnitude,
                y1_values=working_point[1].to(dP_unit).magnitude,
                style_props={'marker': 'o', 'linestyle': 'None', 'color': 'red'}
            )
        if V_max is not None and V_step is not None:
            diagram.x1.scale(lower_limit=0.0, upper_limit=V_max.to(V_unit).magnitude, step=V_step.to(V_unit).magnitude)
        if dP_max is not None and dP_step is not None:
            diagram.y1.scale(lower_limit=0.0, upper_limit=dP_max.to(dP_unit).magnitude, step=dP_step.to(dP_unit).magnitude)
        diagram.x1.add_title(f'volume flow rate, {V_unit}')
        diagram.y1.add_title(f'pump pressure, {dP_unit}')
        return diagram


if __name__ == '__main__':

    from hvac.fluid_flow.utils.misc import head_to_pressure

    pc = PumpCurve.curve_fitting([
        (Quantity(0.0, 'm ** 3 / hr'), head_to_pressure(Quantity(60.0, 'm'))),
        (Quantity(2.4, 'm ** 3 / hr'), head_to_pressure(Quantity(52.0, 'm'))),
        (Quantity(4.2, 'm ** 3 / hr'), head_to_pressure(Quantity(48.0, 'm'))),
        (Quantity(6.0, 'm ** 3 / hr'), head_to_pressure(Quantity(36.0, 'm'))),
    ])

    pump_coeff = pc.coefficients
    for c in pump_coeff:
        print(c)

    pcd = pc.diagram(
        V_ini=Quantity(0.0, 'm ** 3 / hr'),
        V_fin=Quantity(6.0, 'm ** 3 / hr'),
        fig_size=(9, 6),
        dpi=96,
        V_step=Quantity(0.5, 'm ** 3 / hr'),
        V_max=Quantity(7.0, 'm ** 3 / hr'),
        dP_step=head_to_pressure(Quantity(5, 'm')),
        dP_max=head_to_pressure(Quantity(65, 'm')),
        V_unit='L / s',
        dP_unit='bar'
    )
    pcd.show()
