import numpy as np
from scipy.optimize import curve_fit
from hvac import Quantity
from hvac.charts import LineChart


Q_ = Quantity


class PumpCurve:

    def __init__(self):
        """Instantiates an empty `PumpCurve` object."""
        self.name: str = ''
        self._a0: Quantity = Q_(0.0, 'Pa')
        self._a1: Quantity = Q_(0.0, 'Pa * s / m ** 3')
        self._a2: Quantity = Q_(0.0, 'Pa * s ** 2 / m ** 6')

    @classmethod
    def create(
        cls,
        coeffs: list[Quantity, Quantity, Quantity],
        name: str = ''
    ) -> 'PumpCurve':
        """Creates a `PumpCurve` object.

        A pump curve is represented by 2nd order polynomial:
            dP = a0 + a1 * x + a2 * x ** 2

        Parameters
        ----------
        coeffs:
            List of the 3 polynomial coefficients of the pump curve model,
            ordered left to right from zeroth to second order (a0, a1, a2).
        name: optional
            A label to identify the pump curve.
        """
        obj = cls()
        obj.name = name
        obj._a0 = coeffs[0]
        obj._a1 = coeffs[1]
        obj._a2 = coeffs[2]
        return obj

    @classmethod
    def curve_fitting(
        cls,
        coordinates: list[tuple[Quantity, Quantity]],
        name: str = ''
    ) -> 'PumpCurve':
        """Creates a `PumpCurve` object whereby the polynomial coefficients
        a0, a1 and a2 of the pump curve model are determined by curve-fitting
        the given list of pump curve coordinates.

        Parameters
        ----------
        coordinates:
            List of pairs (V, dP) where V is the volume flow rate and dP is
            the corresponding pump pressure differential.
        name: optional
            A label to identify the pump curve.
        """
        V_ax, dP_ax = zip(*coordinates)
        unit_V, unit_dP = V_ax[0].units, dP_ax[0].units
        V_ax = [V.to(V.units).m for V in V_ax]
        dP_ax = [dP.to(dP.units).m for dP in dP_ax]

        def objective(x, a0, a1, a2):
            return a0 + a1 * x + a2 * x ** 2

        popt, *_ = curve_fit(objective, V_ax, dP_ax)
        a0, a1, a2 = popt

        pump_curve = cls()
        pump_curve.name = name
        pump_curve._a0 = Q_(a0, unit_dP)
        pump_curve._a1 = Q_(a1, unit_dP / unit_V)
        pump_curve._a2 = Q_(a2, unit_dP / unit_V ** 2)
        return pump_curve

    def pump_pressure(self, V: Quantity) -> Quantity:
        """Returns the pump pressure differential that corresponds with
        volume flow rate `V`.
        """
        dP = self._a0 + self._a1 * V + self._a2 * V ** 2
        return dP

    def volume_flow_rate(self, dP: Quantity) -> Quantity:
        """Returns the volume flow rate that corresponds with pump pressure
        differential `dP`.
        """
        dP = dP.to(self._a0.units)
        coeffs = [coeff.m for coeff in [self._a2, self._a1, self._a0 - dP]]
        roots = np.roots(coeffs)
        V = roots[roots > 0]
        unit_V = (1 / self._a1.units) * dP.units
        V = Q_(V[0], unit_V)
        return V

    @property
    def coefficients(self) -> tuple[Quantity, ...]:
        """Returns the polynomial coefficients of the 2-nd order pump curve
        model.
        """
        return self._a0, self._a1, self._a2

    @coefficients.setter
    def coefficients(self, values: tuple[Quantity, ...]):
        """Sets the polynomial coefficients of the 2-nd order pump curve
        model. Elements of `values` must be ordered left to right from 0-th
        to 2-nd order (a0, a1, a2).
        """
        self._a0, self._a1, self._a2 = values

    def axes(
        self,
        V_ini: Quantity | None = None,
        V_fin: Quantity | None = None,
        num: int = 50
    ) -> tuple[Quantity, Quantity]:
        """Returns a tuple (`V_ax`, `dP_ax`) of the pump curve coordinates.
        `V_ax` contains a range of volume flow rates (encapsulated within a single
        `Quantity` object, i.e. a `Quantity` array) and `dP_ax` contains the
        range of corresponding pump pressure differentials (also a `Quantity`
        array).

        Parameters
        ----------
        V_ini:
            Volume flow rate of the first coordinate.
        V_fin:
            Volume flow rate of the last coordinate.
        num:
            Number of coordinates.
        """
        if V_ini is None:
            V_ini = Q_(0.0, 'm ** 3 / s')
        if V_fin is None:
            V_fin = self.volume_flow_rate(dP=Q_(0.0, 'Pa'))  # find V for which pump pressure is zero

        V_ax = np.linspace(V_ini, V_fin, num, endpoint=True)
        dP_ax = self._a0 + self._a1 * V_ax + self._a2 * V_ax ** 2
        return V_ax, dP_ax

    def diagram(
        self,
        V_ini: Quantity,
        V_fin: Quantity,
        num: int = 50,
        working_point: tuple[Quantity, Quantity] | None = None,
        **fig_kwargs
    ) -> LineChart:
        """
        Plot pump curve.

        Parameters
        ----------
        V_ini:
            First volume flow rate for which the pump curve is calculated.
        V_fin:
            Final volume flow rate for which the pump curve is calculated.
        num:
            Number of volume flow rates between `V_ini` and `V_fin` for which
            the pump curve is calculated.
        working_point: optional
            Volume flow rate and corresponding pump pressure at given working
            point. The working point will be shown as a red dot on the diagram.
        fig_kwargs:
            Figure settings.

        Following figure settings are available:
        - fig_size: Tuple[int, int]
            The width and height of the figure in inches.
        - dpi: int
            Dots per inch.
        - V_step: Quantity
            The step size between successive volume flow rates on the diagram
            axis.
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
        fig_size: tuple[int, int] = fig_kwargs.get('fig_size', (6, 4))
        dpi: int = fig_kwargs.get('dpi', 96)
        V_step: Quantity = fig_kwargs.get('V_step')
        V_max: Quantity = fig_kwargs.get('V_max')
        dP_step: Quantity = fig_kwargs.get('dP_step')
        dP_max: Quantity = fig_kwargs.get('dP_max')
        V_unit: str = fig_kwargs.get('V_unit', 'm ** 3 / s')
        dP_unit: str = fig_kwargs.get('dP_unit', 'Pa')

        V_ax, dP_ax = self.axes(V_ini, V_fin, num)
        V_ax = V_ax.to(V_unit).m
        dP_ax = dP_ax.to(dP_unit).m
        diagram = LineChart(size=fig_size, dpi=dpi)
        diagram.add_xy_data(
            label=self.name if self.name else 'pump curve',
            x1_values=V_ax,
            y1_values=dP_ax
        )
        if working_point is not None:
            diagram.add_xy_data(
                label='working point',
                x1_values=working_point[0].to(V_unit).m,
                y1_values=working_point[1].to(dP_unit).m,
                style_props={
                    'marker': 'o',
                    'linestyle': 'None', 'color': 'red'
                }
            )
        if V_max is not None and V_step is not None:
            diagram.x1.scale(
                lower_limit=0.0,
                upper_limit=V_max.to(V_unit).m,
                step=V_step.to(V_unit).m
            )
        if dP_max is not None and dP_step is not None:
            diagram.y1.scale(
                lower_limit=0.0,
                upper_limit=dP_max.to(dP_unit).m,
                step=dP_step.to(dP_unit).m
            )
        diagram.x1.add_title(f'volume flow rate, {V_unit}')
        diagram.y1.add_title(f'pump pressure, {dP_unit}')
        return diagram


if __name__ == '__main__':

    from hvac.fluid_flow.utils.misc import head_to_pressure

    pc = PumpCurve.curve_fitting([
        (Q_(0.0, 'm ** 3 / hr'), head_to_pressure(Q_(60.0, 'm'))),
        (Q_(2.4, 'm ** 3 / hr'), head_to_pressure(Q_(52.0, 'm'))),
        (Q_(4.2, 'm ** 3 / hr'), head_to_pressure(Q_(48.0, 'm'))),
        (Q_(6.0, 'm ** 3 / hr'), head_to_pressure(Q_(36.0, 'm'))),
    ])

    pump_coeff = pc.coefficients
    for c in pump_coeff:
        print(c)

    pcd = pc.diagram(
        V_ini=Q_(0.0, 'm ** 3 / hr'),
        V_fin=Q_(6.0, 'm ** 3 / hr'),
        fig_size=(9, 6),
        dpi=96,
        V_step=Q_(0.5, 'm ** 3 / hr'),
        V_max=Q_(7.0, 'm ** 3 / hr'),
        dP_step=head_to_pressure(Q_(5, 'm')),
        dP_max=head_to_pressure(Q_(65, 'm')),
        V_unit='L / s',
        dP_unit='bar'
    )
    pcd.show()
