from abc import ABC, abstractmethod
import numpy as np
from matplotlib.axes import Axes
from matplotlib.ticker import FormatStrFormatter


class Axis(ABC):
    
    def __init__(self, axes: Axes):
        self._axes = axes
    
    @property
    def axes(self) -> Axes:
        return self._axes
    
    @abstractmethod
    def add_title(self, title: str):
        pass
    
    @abstractmethod
    def scale(self, lower_limit: float, upper_limit: float, step: float):
        pass
    
    @abstractmethod
    def format_ticks(self, fmt: str = '%.2f'):
        pass
        

class X1Axis(Axis):
    
    def add_title(self, title: str):
        self._axes.set_xlabel(title)

    def scale(self, lower_limit: float, upper_limit: float, step: float):
        ticks = np.arange(lower_limit, upper_limit, step)
        self._axes.set_xticks(ticks)
        self._axes.set_xlim(ticks[0], ticks[-1])

    def format_ticks(self, fmt: str = '%.2f'):
        self._axes.xaxis.set_major_formatter(FormatStrFormatter(fmt))


class Y1Axis(Axis):

    def add_title(self, title: str):
        self._axes.set_ylabel(title)

    def scale(self, lower_limit: float, upper_limit: float, step: float):
        ticks = np.arange(lower_limit, upper_limit, step)
        self._axes.set_yticks(ticks)
        self._axes.set_ylim(ticks[0], ticks[-1])

    def format_ticks(self, fmt: str = '%.2f'):
        self._axes.yaxis.set_major_formatter(FormatStrFormatter(fmt))


class X2Axis(Axis):
    
    def __init__(self, axes: Axes):
        super().__init__(axes)
        self._twin_axes = axes.twiny()
        self._twin_axes.xaxis.set_ticks_position('bottom')
        self._twin_axes.xaxis.set_label_position('bottom')
        self._twin_axes.spines['bottom'].set_position(('outward', 40))
    
    @property
    def axes(self) -> Axes:
        return self._twin_axes
    
    def add_title(self, title: str):
        self._twin_axes.set_xlabel(title)
    
    def scale(self, lower_limit: float, upper_limit: float, step: float):
        ticks = np.arange(lower_limit, upper_limit, step)
        self._twin_axes.set_xticks(ticks)
        self._twin_axes.set_xlim(ticks[0], ticks[-1])
    
    def format_ticks(self, fmt: str = '%.2f'):
        self._twin_axes.xaxis.set_major_formatter(FormatStrFormatter(fmt))


class Y2Axis(Axis):
    
    def __init__(self, axes: Axes):
        super().__init__(axes)
        self._twin_axes = axes.twinx()

    @property
    def axes(self) -> Axes:
        return self._twin_axes
    
    def add_title(self, title: str):
        self._twin_axes.set_ylabel(title)

    def scale(self, lower_limit: float, upper_limit: float, step: float):
        ticks = np.arange(lower_limit, upper_limit, step)
        self._twin_axes.set_yticks(ticks)
        self._twin_axes.set_ylim(ticks[0], ticks[-1])

    def format_ticks(self, fmt: str = '%.2f'):
        self._twin_axes.yaxis.set_major_formatter(FormatStrFormatter(fmt))
