from typing import Optional, Tuple, List, Callable
from datetime import datetime as DateTime
from datetime import date as Date
from datetime import time as Time
import math
import csv
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from hvac import Quantity
from hvac.charts import LineChart
from .location import Location


Q_ = Quantity


class SunPath:

    def __init__(self, loc: Location, date: Date):
        """Define the sun path at the given location on the given date."""
        self.name = date.strftime('%b %d')
        self.date = date
        self.location = loc
        self._interpolate: Optional[Callable[[float], float]] = None
        self._table: Optional[pd.DataFrame] = None
        self._t_ax: List[DateTime] = []
        self._create_sun_path()

    def _create_sun_path(self):
        t_ax = [Time(h, 0, 0) for h in range(24)]
        t_ax.append(Time(23, 59, 59))
        azi_ax = []; elev_ax = []
        for t in t_ax:
            datetime = DateTime.combine(self.date, t)
            self._t_ax.append(datetime)
            sp = self.location.sun_position(datetime)
            azi_ax.append(sp.azimuth.to('deg').m)
            elev_ax.append(sp.elevation.to('deg').m)
        self._interpolate = interp1d(x=azi_ax, y=elev_ax, kind='cubic')
        self._table = pd.DataFrame(
            data=np.array([[t, azi, elev] for t, azi, elev in zip(t_ax, azi_ax, elev_ax)]),
            columns=['time', 'azimuth', 'elevation']
        )
        self._table.set_index('time', inplace=True)

    def elevation(self, azimuth: Quantity) -> Quantity:
        """Get elevation of the sun at the given azimuth."""
        azi = azimuth.to('deg').m
        elev = self._interpolate(azi)
        return Q_(elev, 'deg')

    @property
    def table(self) -> pd.DataFrame:
        """
        Returns a Pandas DataFrame with the coordinates of the sun path (azimuth and elevation) for each hour of the
        day.
        """
        return self._table

    @property
    def axes(self) -> Tuple[List[DateTime], List[Quantity], List[Quantity]]:
        """Get the azimuth and elevation coordinates of the sun path."""
        azi_ax = [Q_(v, 'deg') for v in self._table['azimuth'].values]
        elev_ax = [Q_(v, 'deg') for v in self._table['elevation'].values]
        return self._t_ax, azi_ax, elev_ax


class HorizonPoint:

    def __init__(self, azimuth: Quantity, elevation: Quantity, planar_distance: Quantity):
        """
        Define a point on a horizon profile. The position of the horizon point is determined by its azimuth angle and
        elevation angle as seen by an observer. The planar distance is the distance between the observer and the horizon
        point measured in a horizontal plane (top view).
        """
        self._azi = azimuth.to('rad').m
        self._elev = elevation.to('rad').m
        self._dist = planar_distance.to('m').m

    def move_origin(self, deltaS: Quantity, deltaE: Quantity, deltaH: Quantity):
        """
        Recalculates the position of the horizon point when the origin is moved.

        The position of a horizon point is relative to its observer whose standing in the origin of the reference frame.
        When the observer and a such the origin of the reference frame moves to a different position, the observed
        position of the horizon point will also change.

        Parameters
        ----------
        deltaS : movement along the North-South-axis. Movement is positive in the South direction.
        deltaE : movement along the West-East-axis. Movement is positive in the East direction.
        deltaH : change of height
        """
        S_orig, E_orig, H_orig = self._rectangular_coords(self._azi, self._elev, self._dist)
        S_new = S_orig - deltaS.to('m').m
        E_new = E_orig - deltaE.to('m').m
        H_new = H_orig - deltaH.to('m').m
        self._azi, self._elev, self._dist = self._polar_coords(S_new, E_new, H_new)

    @staticmethod
    def _rectangular_coords(azi, elev, dist) -> Tuple[float, ...]:
        azi = math.pi - azi
        S = dist * math.cos(azi)
        E = dist * math.sin(azi)
        H = dist * math.tan(elev)
        return S, E, H

    @staticmethod
    def _polar_coords(S: float, E: float, H: float) -> Tuple[float, ...]:
        dist = math.sqrt(S ** 2 + E ** 2)
        elev = math.atan(H / dist)
        cos_azi = abs(S) / dist
        if S >= 0.0 and E >= 0.0:  # obstacle point in quadrant SE
            azi_offset = math.pi
            sign = -1
        elif S >= 0.0 > E:  # quadrant SW
            azi_offset = math.pi
            sign = 1
        elif S <= 0.0 <= E:  # quadrant NE
            azi_offset = 0.0
            sign = 1
        else:  # quadrant NW
            azi_offset = 2 * math.pi
            sign = -1
        return azi_offset + sign * math.acos(cos_azi), elev, dist

    @property
    def azimuth(self) -> Quantity:
        """Get the azimuth angle of the horizon point."""
        return Q_(self._azi, 'rad')

    @property
    def elevation(self) -> Quantity:
        """Get the elevation angle of the horizon point."""
        return Q_(self._elev, 'rad')


class _Segment:

    def __init__(self, start_point: HorizonPoint, end_point: HorizonPoint):
        self.start_point = start_point
        self.end_point = end_point
        self._interpolate = interp1d(
            x=[start_point.azimuth.to('deg').m, end_point.azimuth.to('deg').m],
            y=[start_point.elevation.to('deg').m, end_point.elevation.to('deg').m]
        )

    def elevation(self, azimuth: Quantity) -> Quantity:
        """Get elevation [deg] of point on segment of which the azimuth [deg] is given."""
        elev = self._interpolate(azimuth.to('deg').m)
        return Q_(elev, 'deg')


class HorizonProfile:
    """
    Represents the contour of the horizon formed by obstacles that may block the sun rays at certain times during
    the day at the location under interest. A horizon profile is tied to a specific position from where the horizon
    profile is observed.
    """

    def __init__(self, name: str, points: List[HorizonPoint]):
        """
        A horizon profile is described by a list of `HorizonPoint` objects. These points are interconnected by
        line segments.

        Parameters
        ----------
        name: str
            Identifier for the horizon profile.
        points: List[HorizonPoint]
            A list of all horizon points that define the horizon profile. Horizon points will internally be sorted
            according to their azimuth angle in ascending order.
        """
        self.name = name
        self.points = sorted(points, key=lambda pt: pt.azimuth)
        num_of_segments = len(self.points) - 1
        self._segments = [_Segment(self.points[i], self.points[i + 1]) for i in range(num_of_segments)]

    def elevation(self, azimuth: Quantity) -> Quantity:
        """
        Get the elevation angle of the horizon profile at the specified azimuth angle.

        Using linear interpolation between two successive horizon points, it is possible to calculate the elevation of
        intermediary points on the horizon profile.
        """
        for segment in self._segments:
            if segment.start_point.azimuth <= azimuth <= segment.end_point.azimuth:
                return segment.elevation(azimuth)
        else:
            return Q_(0.0, 'deg')

    @property
    def axes(self) -> Tuple[List[float], List[float]]:
        """Get the azimuth and elevation coordinates of the horizon profile in degrees."""
        return (
            [pt.azimuth.to('deg').m for pt in self.points],
            [pt.elevation.to('deg').m for pt in self.points]
        )

    def move_origin(self, deltaS: Quantity, deltaE: Quantity, deltaH: Quantity):
        """
        Calculate the new coordinates of the horizon profile when the observer moves to a new position.

        The position of an obstacle point measured by an observer will depend on his position relative to the obstacle
        point. If the observer moves to a new position, the coordinates of the obstacle point will change also.

        Parameters
        ----------
        deltaS: Quantity
            The distance the observer moves along the S-axis relative to the original reference frame. A positive
            movement is in the direction of the S-axis; a movement in the opposite direction is negative.
        deltaE: Quantity
            The distance the observer moves along the E-axis relative to the original reference frame. A positive
            movement is in the direction of the E-axis; a movement in the opposite direction is negative.
        deltaH: Quantity
            The height difference between the original position and the new position. A positive value
            means that the observer is at a higher position than his original position; a negative value means that
            the observer is now at a lower position than before.

        """
        for pt in self.points:
            pt.move_origin(deltaS, deltaE, deltaH)
        # recreate line segments of horizon profile
        num_of_segments = len(self.points) - 1
        self._segments = [_Segment(self.points[i], self.points[i + 1]) for i in range(num_of_segments)]


def plot_sun_paths(sun_paths: List[SunPath], horizon_profile: Optional[HorizonProfile] = None, **kwargs):
    """
    Draw one or more sun paths on a cartesian diagram, possibly combined with a horizon profile as seen at the specific
    location.

    Parameters
    ----------
    sun_paths: List[SunPath]
        A list of `SunPath` objects.
    horizon_profile: HorizonProfile, optional
        A horizon profile.
    **kwargs: dict
        Optional keyword arguments for setting the size (keyword `fig_size`) and dpi (keyword `dpi`) of the diagram.
        Keyword `fig_size` expects a 2-tuple of *int* representing respectively the width and the height of the
        diagram. Keyword `dpi` expects an *int* representing the dots-per-inch of the diagram plot.
    """
    graph = LineChart(
        size=kwargs.get('fig_size', (10, 6)),
        dpi=kwargs.get('dpi', 96),
    )
    for path in sun_paths:
        t_axis, azi_axis, elev_axis = path.axes
        graph.add_xy_data(
            label=path.name,
            x1_values=[azi.to('deg').m for azi in azi_axis],
            y1_values=[elev.to('deg').m for elev in elev_axis],
            style_props={'marker': 'o'}
        )
    if horizon_profile is not None:
        azi_axis, elev_axis = horizon_profile.axes
        graph.add_xy_data(
            label=horizon_profile.name,
            x1_values=azi_axis,
            y1_values=elev_axis,
            style_props={'fill': {'color': 'black'}}
        )
    graph.add_legend(anchor='upper center', position=(0.5, -0.2), columns=7)
    graph.x1.add_title('azimuth')
    graph.x1.scale(0, 360, 20)
    graph.y1.add_title('elevation')
    graph.y1.scale(0, 90, 5)
    return graph


def create_horizon_profile(name: str, file_path: str) -> HorizonProfile:
    """
    Create a horizon profile from CSV-file.

    The CSV-file has 3 columns. The first row is considered as a header and will be skipped.

    - 1st column : azimuth angle of horizon points in decimal degrees.
    - 2nd column : elevation angle of horizon points in decimal degrees.
    - 3rd column : planar distance between horizon points and observer in metres.

    Parameters
    ----------
    name: str
        Identifier for the horizon profile.
    file_path: str
        Filepath of the CSV-file.

    Returns
    -------
    HorizonProfile
    """
    horizon_points = []
    with open(file_path, newline='') as fh:
        reader = csv.reader(fh)
        for i, row in enumerate(reader):
            if i == 0:
                continue  # skip first row (header)
            else:
                horizon_point = HorizonPoint(
                    azimuth=Q_(float(row[0]), 'deg'),
                    elevation=Q_(float(row[1]), 'deg'),
                    planar_distance=Q_(float(row[2]), 'm')
                )
                horizon_points.append(horizon_point)
    horizon_profile = HorizonProfile(name, horizon_points)
    return horizon_profile
