from datetime import date as Date
from datetime import time as Time
import numpy as np
from scipy.interpolate import interp1d
from hvac import Quantity
from hvac.sun.surface import Location
from hvac.charts.matplotlibwrapper import LineChart


Q_ = Quantity


class SunPath:
    """Represents a sun path on a given location and date.

    A sun path is a curve of solar altitude angle versus solar azimuth angle
    between sunrise and sunset at the given date, describing the position of
    the sun during daytime.

    Attributes
    ----------
    gamma: Quantity
        List of the solar azimuth angles between sunrise and sunset, set out
        along the horizontal axis of the sun chart.
    alpha: Quantity
        List of the solar altitude angles between sunrise and sunset, set out
        along the vertical axis of the sun chart.

    Methods
    -------
    __call__(gamma: Quantity):
        Returns the solar altitude angle that corresponds with the given
        solar azimuth angle using linear interpolation.
    """
    def __init__(self, location: Location, date: Date | None = None):
        """Creates a `SunPath` instance.

        Parameters
        ----------
        location:
            The location where the sun's trajectory is to be determined.
        date:
            The date the sun's trajectory is to be determined.
        """
        self.location = location
        if date is not None: self.location.date = date
        self.gamma: Quantity | None = None
        self.alpha: Quantity | None = None
        self._interp = None
        self.__create__()

    def __create__(self):
        sr = self.location.sunrise
        ss = self.location.sunset
        _gamma, _alpha = [], []
        tl = [sr]
        tl += [
            t
            for h in range(0, 24) for m in range(0, 60, 30)
            if sr < (t := Time(h, m, 0)) < ss
        ]
        tl += [ss]
        for t in tl:
            self.location.solar_time = t
            _gamma.append(self.location.sun.gamma.to('deg').m)
            _alpha.append(self.location.sun.alpha.to('deg').m)
        self._interp = interp1d(_gamma, _alpha)
        self.gamma = Q_(_gamma, 'deg')
        self.alpha = Q_(_alpha, 'deg')

    def __call__(self, gamma: Quantity) -> Quantity:
        """Returns the solar altitude angle that corresponds with the given
        solar azimuth angle using linear interpolation.
        """
        return Q_(self._interp(gamma.to('deg').m), 'deg')


class ShadingPoint:
    """Represents a point along the contour of an obstruction that may
    block the path from the sun to the point where an observer is standing.

    Attributes
    ----------
    gamma:
        shading point azimuth angle with respect to the observer
    alpha:
        shading point altitude angle with respect to the observer
    distance:
        horizontal distance between the shading point and the observer
    """
    def __init__(
        self,
        gamma: Quantity,
        alpha: Quantity,
        distance: Quantity
    ) -> None:
        """Creates a `ShadingPoint` instance.

        Parameters
        ----------
        gamma:
            the shading point azimuth angle as seen by the observer facing south;
            positive from south to west, negative from south to east
        alpha:
            shading point altitude angle as seen by the observer;
            always positive
        distance:
            horizontal distance between the observer and the shading point;
            always positive
        """
        self.gamma = gamma
        self.alpha = alpha
        self.distance = distance

    def move_origin(
        self,
        d_south: Quantity,
        d_west: Quantity,
        d_height: Quantity
    ) -> None:
        """Given the coordinates of the shading point as seen from an observer at
        a given position 1, recalculates the coordinates of the shading point as
        seen from an observer standing at a different position 2.

        Parameters
        ----------
        d_south:
            The distance between observer position 2 and observer position 1
            along the north-south axis; positive when position 2 is southwards
            from position 1, otherwise negative.
        d_west:
            The distance between observer position 2 and observer position 1
            along the east-west axis; positive when position 2 is westwards from
            position 1, otherwise negative.
        d_height:
            The height difference between observer position 2 and observer
            position 1; positive when position 2 is higher than position 1,
            otherwise negative.

        Returns
        -------
        None
        """
        _s1, _w1, _h1 = self._rect_coords(self.gamma, self.alpha, self.distance)
        _s2 = _s1 - d_south
        _w2 = _w1 - d_west
        _h2 = _h1 - d_height
        self.gamma, self.alpha, self.distance = self._polar_coords(_s2, _w2, _h2)

    @staticmethod
    def _rect_coords(
        _g: Quantity,
        _a: Quantity,
        _d: Quantity
    ) -> tuple[Quantity, ...]:
        """Returns the rectangular coordinates of the shading point with respect
        to the position of the observer facing south. The positive sense of the
        horizontal, rectangular coordinate system is directed to the south
        (north-south axis) and to the west (east-west axis).

        Parameters
        ----------
        _g:
            the shading point azimuth angle as seen by the observer facing south;
            positive from south to west, negative from south to east
        _a:
            the shading point altitude angle
        _d:
            distance between the observer and the shading point.

        Returns
        -------
        The south-north coordinate, the east-west coordinate, and the height
        of the shading point with respect to the observer facing south in the
        origin of the rectangular coordinate system.
        """
        _s = _d * np.cos(_g)
        _w = _d * np.sin(_g)
        _h = _d * np.tan(_a)
        return _s, _w, _h

    @staticmethod
    def _polar_coords(
        _s: Quantity,
        _w: Quantity,
        _h: Quantity
    ) -> tuple[Quantity, ...]:
        """Returns the polar coordinates of the shading point with respect
        to the position of the observer.

        Parameters
        ----------
        _s:
            coordinate on the north-south axis of the shading point with respect
            to the position of the observer facing south.
        _w:
            coordinate on the east-west axis of the shading point with respect
            to the position of the observer facing south.
        _h:
            height of the shading point as seen from the observer.

        Returns
        -------
        The shading point azimuth angle, altitude angle and horizontal distance
        between the observer and the shading point, in that order.
        """
        _d = np.sqrt(_w ** 2 + _s ** 2)
        _a = np.arctan(_h / _d)
        _g = np.arctan(_w / _s)
        return _g, _a, _d


class _ShadingSegment:
    """Represents a shading line segment along the contour of an obstruction
    that may block the path from the sun to the point where an observer is
    standing."""

    def __init__(self, start: ShadingPoint, end: ShadingPoint) -> None:
        """Creates a `ShadingSegment` instance.

        Parameters
        ----------
        start:
            start point of the segment
        end:
            end point of the segment
        """
        self.start = start
        self.end = end
        self._interp = interp1d(
            x=[start.gamma.to('deg').m, end.gamma.to('deg').m],
            y=[start.alpha.to('deg').m, end.alpha.to('deg').m]
        )

    def altitude(self, gamma: Quantity) -> Quantity:
        """Get the altitude of a shading point on the segment when its
        azimuth angle is given."""
        _a = self._interp(gamma.to('deg').m)
        return Q_(_a, 'deg')


class ShadingProfile:
    """Represents the contour of an obstruction that may block the path from the
    sun to the point where an observer is standing.

    A shading profile is formed by line segments that connect the given points
    on the contour of the obstruction.
    """

    def __init__(self, ID: str, pts: list[ShadingPoint]) -> None:
        """Creates a `ShadingProfile` instance.

        Parameters
        ----------
        ID:
            identifier for the shading profile
        pts:
            list of `ShadingPoint` objects situated along the contour of the
            obstruction
        """
        self.ID = ID
        self.pts = sorted(pts, key=lambda pt: pt.gamma)
        self._segments = [
            _ShadingSegment(self.pts[i], self.pts[i+1])
            for i in range(len(self.pts) - 1)
        ]

    def altitude(self, gamma: Quantity) -> Quantity:
        """Returns the altitude angle of the shading point with given azimuth
        angle on the contour of an obstruction. Returns 0° if the shading point
        with given azimuth angle"""
        for seg in self._segments:
            if seg.start.gamma <= gamma <= seg.end.gamma:
                return seg.altitude(gamma)
        else:
            return Q_(0.0, 'deg')

    @property
    def gamma(self) -> Quantity:
        """Returns azimuth angles of the shading profile."""
        _g = Q_([pt.gamma.to('deg').m for pt in self.pts], 'deg')
        return _g

    @property
    def alpha(self) -> Quantity:
        """Returns altitude angles of the shading profile."""
        _a = Q_([pt.alpha.to('deg').m for pt in self.pts], 'deg')
        return _a

    def move_origin(
        self,
        d_south: Quantity,
        d_west: Quantity,
        d_height: Quantity
    ) -> None:
        """Given the coordinates of the shading points as seen from an observer
        at a given position 1, recalculates the coordinates of the shading points
        as seen from an observer standing at a different position 2.

        Parameters
        ----------
        d_south:
            The distance between observer position 2 and observer position 1
            along the north-south axis; positive when position 2 is southwards
            from position 1, otherwise negative.
        d_west:
            The distance between observer position 2 and observer position 1
            along the east-west axis; positive when position 2 is westwards from
            position 1, otherwise negative.
        d_height:
            The height difference between observer position 2 and observer
            position 1; positive when position 2 is higher than position 1,
            otherwise negative.

        Returns
        -------
        None
        """
        for pt in self.pts:
            pt.move_origin(d_south, d_west, d_height)
        self._segments = [
            _ShadingSegment(self.pts[i], self.pts[i + 1])
            for i in range(len(self.pts) - 1)
        ]


class _SolarTimeIsoLines(list):

    def __init__(self, location: Location, dates: list[Date]):
        super().__init__()
        self.location = location
        self.dates = dates
        self.__create__()

    def __create__(self):
        hl = [h for h in range(6, 19)]
        for h in hl:
            _gamma, _alpha = [], []
            for d in self.dates:
                self.location.date = d
                self.location.solar_time = Time(h, 0, 0)
                _gamma.append(self.location.sun.gamma.to('deg').m)
                _alpha.append(self.location.sun.alpha.to('deg').m)
            self.append((h, Q_(_gamma, 'deg'), Q_(_alpha, 'deg')))


class SunChart:
    """Class for plotting a sun chart without or with a shading profile."""

    def __init__(self, location: Location):
        """Creates `SunChart` instance.

        Parameters
        ----------
        location:
            Location for which the sun chart will be drawn.
        """
        self.location = location
        self.dates = [
            Date(2023, 1, 17),
            Date(2023, 2, 16),
            Date(2023, 3, 16),
            Date(2023, 4, 15),
            Date(2023, 5, 15),
            Date(2023, 6, 11),
            Date(2023, 7, 17),
            Date(2023, 8, 16),
            Date(2023, 9, 15),
            Date(2023, 10, 15),
            Date(2023, 11, 14),
            Date(2023, 12, 10)
        ]
        self.sun_paths = [
            SunPath(location, date) for date in self.dates
        ]
        self.solar_time_isolines = _SolarTimeIsoLines(location, self.dates)
        self.shading_profile: ShadingProfile | None = None

    def add_shading_profile(self, spf: ShadingProfile):
        """Adds a shading profile to be drawn on the sun chart."""
        self.shading_profile = spf

    def plot(
        self,
        fig_size: tuple[int, int] = (10, 6),
        dpi: int = 150,
        legend_on: bool = True
    ) -> LineChart:
        """Returns a `LineChart` object of the sun chart."""
        sun_chart = LineChart(fig_size, dpi)
        for d, sun_path in zip(self.dates, self.sun_paths):
            sun_chart.add_xy_data(
                label=f'sun path {d.month}/{d.day}',
                x1_values=sun_path.gamma,
                y1_values=sun_path.alpha,
            )
        for isoline in self.solar_time_isolines:
            sun_chart.add_xy_data(
                label=f'solar time {isoline[0]}',
                x1_values=isoline[1],
                y1_values=isoline[2],
                style_props={'linestyle': '--'}
            )
        if self.shading_profile is not None:
            sun_chart.y1.axes.fill_between(
                x=self.shading_profile.gamma,
                y1=self.shading_profile.alpha,
                color='black',
                edgecolor='black'
            )
        sun_chart.x1.add_title('solar azimuth, deg')
        sun_chart.y1.add_title('solar altitude, deg')
        sun_chart.x1.scale(-120, 150, 30)
        sun_chart.y1.scale(0, 95, 5)
        if legend_on:
            sun_chart.add_legend(columns=6)
        return sun_chart


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.get_cmap('tab20').colors)

    loc = Location(fi=Q_(40, 'deg'))

    # shading profile of an east-west oriented wall with the observer to the
    # north of the wall at a distance of 10 m.
    shading_profile = ShadingProfile(
        ID='south wall',
        pts=[
            ShadingPoint(Q_(-71.6, 'deg'), Q_(4.52, 'deg'), Q_(31.6, 'm')),
            ShadingPoint(Q_(-45.0, 'deg'), Q_(10.0, 'deg'), Q_(14.1, 'm')),
            ShadingPoint(Q_(0.0, 'deg'), Q_(14.0, 'deg'), Q_(10.0, 'm')),
            ShadingPoint(Q_(45.0, 'deg'), Q_(10.0, 'deg'), Q_(14.1, 'm')),
            ShadingPoint(Q_(71.6, 'deg'), Q_(4.52, 'deg'), Q_(31.6, 'm'))
        ]
    )

    sun_chart = SunChart(loc)
    sun_chart.add_shading_profile(shading_profile)

    # move the observer 5 m closer to the wall (i.e. move the origin 5 m towards
    # the south)
    sun_chart.shading_profile.move_origin(
        d_south=Q_(5, 'm'),
        d_west=Q_(0, 'm'),
        d_height=Q_(0, 'm')
    )

    drawing = sun_chart.plot()
    drawing.show()
