from typing import Any
from copy import deepcopy
from hvac import Quantity


class Schedule:

    def __init__(self):
        self.ID: str = ''
        self.dt_hr: float = 0.0
        self.base_value: int | Quantity | bool = 0
        self.time_axis: list[float] = []
        self._time_intervals: list[tuple[float, float]] = []
        self.value_axis: list[int | Quantity] = []

    @classmethod
    def create(
        cls,
        ID: str,
        base_value: int | Quantity | bool,
        dt_hr: float = 1.0
    ) -> 'Schedule':
        """
        Create a `Schedule` instance.

        Parameters
        ----------
        ID:
            Identifier of schedule.
        base_value:
            The base value of the schedule.
        dt_hr: default 1.0 hr
            Time step of schedule in decimal hours.
        """
        schedule = cls()
        schedule.ID = ID
        schedule.base_value = base_value
        schedule.dt_hr = dt_hr
        schedule._set_up()
        return schedule

    def _set_up(self):
        k_max = int(24 / self.dt_hr)
        self.time_axis = [k * self.dt_hr for k in range(k_max)]
        self._time_intervals = [
            (t1, t2)
            for t1, t2 in zip(self.time_axis, self.time_axis[1:] + [24])
        ]
        self.value_axis = [deepcopy(self.base_value)] * len(self.time_axis)

    def set_value(self, t: float, val: int | Quantity | bool):
        """Set value in schedule at time t in decimal hours."""
        for t_ in self.time_axis:
            if t_ == t:
                i = self.time_axis.index(t_)
                self.value_axis[i] = val
                break

    def get_schedule(self) -> dict[str, list[float | int | Quantity | bool]]:
        return {
            't': self.time_axis,
            'value': self.value_axis
        }

    def __call__(
        self,
        t_sec: float | None = None,
        t_hr: float | None = None
    ) -> Any:
        if t_hr is None:
            t_hr = t_sec / 3600.0
        for t_interval in self._time_intervals:
            if t_interval[0] <= t_hr < t_interval[1]:
                i = self.time_axis.index(t_interval[0])
                value = self.value_axis[i]
                return value
            elif t_hr == t_interval[1]:
                i = self.time_axis.index(t_interval[1])
                value = self.value_axis[i]
                return value
        else:
            raise ValueError(f'no value at time {t_sec} seconds')


class OccupancySchedule(Schedule):
    """Schedule for number of persons in a space."""
    pass


class OnOffSchedule(Schedule):
    """Schedule to indicate when equipment or appliances are on or off."""
    pass


class TemperatureSchedule(Schedule):
    """Schedule for the setpoint temperature in a space."""

    def __call__(
        self,
        t_sec: float | None = None,
        t_hr: float | None = None
    ) -> float:
        T = super().__call__(t_sec, t_hr)
        return T.to('degC').m
