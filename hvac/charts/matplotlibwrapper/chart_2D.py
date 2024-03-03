from typing import Optional, Tuple, Dict, Any, Iterable, Union
from dataclasses import dataclass
from abc import ABC, abstractmethod
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from .axes import X1Axis, Y1Axis, X2Axis, Y2Axis


@dataclass
class Legend:
    y1: Y1Axis
    y2: Optional[Y2Axis] = None
    anchor: str = 'upper center'
    position: Tuple[float, float] = (0.5, -0.1)
    columns: int = 2

    def draw(self):
        y1_handles, y1_labels = self.y1.axes.get_legend_handles_labels()
        if isinstance(self.y2, Y2Axis):
            y2_handles, y2_labels = self.y2.axes.get_legend_handles_labels()
            y1_handles.extend(y2_handles)
            y1_labels.extend(y2_labels)
        self.y1.axes.legend(
            y1_handles, y1_labels,
            loc=self.anchor, ncol=self.columns, bbox_to_anchor=self.position
        )


class Chart(ABC):

    def __init__(
        self,
        size: Optional[Tuple[float, float]] = None,
        dpi: Optional[int] = None,
        constructs: Optional[Tuple[Figure, Axes]] = None
    ):
        if constructs is not None:
            self.figure = constructs[0]
            self.axes = constructs[1]
        else:
            self.figure, self.axes = plt.subplots(
                figsize=size,
                dpi=dpi,
                layout='constrained'
            )
        self.x1: X1Axis = X1Axis(self.axes)
        self.y1: Y1Axis = Y1Axis(self.axes)
        self.x2: Optional[X2Axis] = None
        self.y2: Optional[Y2Axis] = None
        self.datasets: Dict[str, Any] = {}
        self.legend: Optional[Legend] = None
        self._y1_img = None
        self._y2_img = None

    def add_x2_axis(self):
        """Add a secondary x-axis below the main x-axis."""
        self.x2 = X2Axis(self.axes)

    def add_y2_axis(self):
        """Add a secondary y-axis at the right of the chart."""
        self.y2 = Y2Axis(self.axes)

    def add_xy_data(
        self,
        label: str,
        x1_values: Optional[Union[Iterable, np.ndarray]] = None,
        y1_values: Optional[Union[Iterable, np.ndarray]] = None,
        x2_values: Optional[Union[Iterable, np.ndarray]] = None,
        y2_values: Optional[Union[Iterable, np.ndarray]] = None,
        style_props: Optional[Dict[str, Any]] = None
    ):
        """
        Add x- and y- data to chart for drawing.

        If the chart has a secondary x-axis (x2) and/or a secondary y-axis (y2),
        also `x2_values` and `y2_values` can be added to the chart.
        `style_props` can be a dictionary with values for properties that style
        the plot (e.g. {'marker': 'o', 'linestyle': 'none'})
        """
        self.datasets[label] = {
            'x1_values': x1_values,
            'y1_values': y1_values,
            'x2_values': x2_values,
            'y2_values': y2_values,
            'style_props': style_props or {}
        }

    @abstractmethod
    def _draw_xy_data(self):
        pass

    def add_legend(
        self,
        anchor: str = 'upper center',
        position: Tuple[float, float] = (0.5, -0.1),
        columns: int = 2
    ):
        """
        Add legend to chart.

        Parameters
        ----------
        anchor: {'upper left', 'upper center', 'upper right', 'center right',
                 'lower right', 'lower center', 'lower left', 'center left',
                 'center', 'best'}
            Reference point on the border of the legend for positioning the
            legend on the chart.
        position:
            Tuple with x and y coordinates of the anchor position with respect
            to the origin of the axes. By default, the legend is positioned at
            the bottom and at the center of the chart, under the x-axis.
        columns:
            Number of label columns the legend list is to be divided in.
            Default is 2 columns.
        """
        self.legend = Legend(self.y1, self.y2, anchor, position, columns)

    def add_title(self, title: str):
        """Add title at top of the chart."""
        self.axes.set_title(title)

    def add_note(
        self,
        text: str,
        x_pos: float = 0.05,
        y_pos: float = 0.95,
        font_size: int = 12,
        vert_align: str = 'top',
        hor_align: str = 'left',
        use_normalized_coordinates: bool = True,
        bbox_props: dict | None | bool = None
    ) -> None:
        """Add some text within the chart.

        Parameters
        ----------
        text:
            Text to display.
        x_pos:
            Horizontal position of the text box anchor point.
        y_pos:
            Vertical position of the text box anchor point.
        font_size:
            Size of the text font.
        vert_align: {'top', 'center', 'bottom', 'baseline', 'center_baseline'}
            Sets the vertical position of the text box anchor point.
        hor_align: {'left', 'center', 'right'}
            Sets the horizontal position of the text box anchor point.
        use_normalized_coordinates:
            Instead of data coordinates, use normalized axes coordinates, i.e.
            the bottom-left corner of the axes is at (0,0) and the top-right
            corner is at (1, 1).
        bbox_props:
            If `None` or `False`, no bounding box is drawn around the text. If
            `True` a bounding box with default properties is drawn around the
            text. Otherwise, pass the bounding box properties as a dict. The
            default properties of the bounding box are specified as:
            `{'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.5}`.
        """
        kwargs = {
            'fontsize': font_size,
            'horizontalalignment': hor_align,
            'verticalalignment': vert_align,
        }
        if use_normalized_coordinates:
            kwargs['transform'] = self.y1.axes.transAxes
        if isinstance(bbox_props, dict):
            kwargs['bbox'] = bbox_props
        elif bbox_props is True:
            kwargs['bbox'] = {
                'boxstyle': 'round',
                'facecolor': 'wheat',
                'alpha': 0.5
            }
        self.y1.axes.text(x_pos, y_pos, text, **kwargs)

    def draw(self, with_grid: bool = True):
        """Only draw the chart (but don't show it)."""
        self._draw_xy_data()
        if self.legend:
            self.legend.draw()
        self.axes.grid(with_grid)

    def show(self, with_grid: bool = True):
        """Draw and show the chart."""
        self.draw(with_grid)
        plt.show()

    def save(
        self, name: str,
        location: Optional[str] = None,
        fmt: str = 'png',
        with_grid: bool = True
    ):
        """Draw and save the chart on disk."""
        self.draw(with_grid)
        location = location or os.getcwd()
        path = os.path.join(location, f'{name}.{fmt}')
        self.figure.savefig(path, bbox_inches='tight')
        plt.close(self.figure)

    @property
    def y1_img(self):
        return self._y1_img

    @property
    def y2_img(self):
        return self._y2_img


class LineChart(Chart):

    def _draw_xy_data(self):
        for label, dataset in self.datasets.items():
            if dataset['x1_values'] is not None:
                if dataset['y1_values'] is not None:
                    self._y1_img = self.y1.axes.plot(
                        dataset['x1_values'],
                        dataset['y1_values'],
                        label=label,
                        **dataset['style_props']
                    )
                if dataset['y2_values'] is not None:
                    self._y2_img = self.y2.axes.plot(
                        dataset['x1_values'],
                        dataset['y2_values'],
                        label=label,
                        **dataset['style_props']
                    )


class FilledLineChart(Chart):

    def _draw_xy_data(self):
        for label, dataset in self.datasets.items():
            if dataset['x1_values'] is not None:
                if dataset['y1_values'] is not None:
                    self._y1_img = self.y1.axes.fill_between(
                        dataset['x1_values'],
                        dataset['y1_values'],
                        label=label,
                        **dataset['style_props']
                    )
                if dataset['y2_values'] is not None:
                    self._y2_img = self.y2.axes.fill_between(
                        dataset['x1_values'],
                        dataset['y2_values'],
                        label=label,
                        **dataset['style_props']
                    )


class BarChart(Chart):

    def _draw_xy_data(self):
        for label, dataset in self.datasets.items():
            if dataset['x1_values'] is not None:
                if dataset['y1_values'] is not None:
                    self._y1_img = self.y1.axes.bar(
                        dataset['x1_values'],
                        dataset['y1_values'],
                        label=label,
                        **dataset['style_props']
                    )
                if dataset['y2_values'] is not None:
                    self._y2_img = self.y2.axes.bar(
                        dataset['x1_values'],
                        dataset['y1_values'],
                        label=label,
                        **dataset['style_props']
                    )
