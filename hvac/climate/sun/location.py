from typing import Optional
from dataclasses import dataclass
from datetime import datetime as DateTime
from datetime import date as Date
from datetime import time as Time
from datetime import timedelta as TimeDelta
import pytz
import astral
import astral.sun
from hvac import Quantity
from .solar_time import solar_time, hour_angle


DEFAULT_TZ = "Etc/UTC"
Q_ = Quantity


@dataclass
class SunPosition:
    azimuth: Quantity | None = None
    elevation: Quantity | None = None
    zenith: Quantity | None = None

    def __str__(self):
        return (
            f"azimuth {self.azimuth.to('deg'):~P.4f}° (CW N) | "
            f"elevation {self.elevation.to('deg'):~P.4f}° (zenith {self.zenith.to('deg'):~P.4f}°)"
        )


class Location:
    """
    Represents a geographic location on Earth.
    """

    def __init__(
        self,
        name: str,
        lat: Quantity,
        lon: Quantity,
        alt: Quantity = Q_(0.0, 'm'),
        tz: str = DEFAULT_TZ
    ):
        """
        Parameters
        ----------
        name: str
            A name for the location
        lat: Quantity
            Latitude of the location (N: positive, S: negative).
        lon: Quantity
            Longitude of the location (E: positive, W: negative).
        alt: Quantity
            Altitude of the location.
        tz: str
            Timezone of the location (see tz database).
        """
        self.name = name
        self.lat = lat
        self.lon = lon
        self.alt = alt
        self.tz = tz
        self._tzinfo = pytz.timezone(tz)
        self._observer = astral.Observer(
            latitude=lat.to('deg').m,
            longitude=lon.to('deg').m,
            elevation=alt.to('m').m
        )
        self._datetime: Optional[DateTime] = None

    def __str__(self):
        return (
            f"{self.name} | "
            f"lat = {self.lat('deg'):.4f}° | "
            f"lon = {self.lon('deg'):.4f}° | "
            f"alt = {self.alt('m'):.0f} m | "
            f"tz = {self.tz}"
        )

    def _localize(self, datetime: DateTime) -> DateTime:
        if datetime.tzinfo is None:
            # make naive datetime timezone aware using the timezone of the location
            return self._tzinfo.localize(datetime, is_dst=True)
        elif datetime.tzinfo.zone != self._tzinfo.zone:
            # convert datetime to timezone of the location
            return datetime.astimezone(self._tzinfo)
        else:
            return datetime

    def sun_position(self, datetime: DateTime, tz_aware: bool = True) -> SunPosition:
        """
        Get position of the sun at the given datetime.

        Parameters
        ----------
        datetime: DateTime
            Naive Python datetime at which the position of the sun is to be calculated.
        tz_aware: bool, optional
            Flag to indicate if the datetime should be considered as a local standard time using the time zone of the
            location (`True`) or in UTC (`False`). Default value is `True`.

        Returns
        -------
        SunPosition

        The position of the sun is specified by its azimuth angle, measured clockwise from North, and its elevation or
        zenith angle.
        """
        if tz_aware: datetime = self._localize(datetime)
        azi = astral.sun.azimuth(self._observer, datetime)
        elev = astral.sun.elevation(self._observer, datetime)
        zen = 90.0 - elev
        return SunPosition(
            azimuth=Q_(azi, 'deg'),  # measured clockwise from North
            elevation=Q_(elev, 'deg'),
            zenith=Q_(zen, 'deg')
        )

    def sunrise(self, date: Date, tz_aware: bool = True) -> DateTime:
        """
        Get sunrise time at given date.

        If `tz_aware` is `True`, time will be returned in local standard time of the location, otherwise in UTC.
        """
        tzinfo = self._tzinfo if tz_aware else None
        return astral.sun.sunrise(self._observer, date, tzinfo=tzinfo)

    def noon(self, date: Date, tz_aware: bool = True) -> DateTime:
        """
        Get time of solar noon at given date.

        If `tz_aware` is `True`, time will be returned in local standard time of the location, otherwise in UTC.
        """
        tzinfo = self._tzinfo if tz_aware else None
        return astral.sun.noon(self._observer, date, tzinfo=tzinfo)

    def sunset(self, date: Date, tz_aware: bool = True) -> DateTime:
        """
        Get time of sunset at current day.

        If `tz_aware` is `True`, time will be returned in local standard time of the location, otherwise in UTC.
        """
        tzinfo = self._tzinfo if tz_aware else None
        return astral.sun.sunset(self._observer, date, tzinfo=tzinfo)

    def daylight_duration(self, date: Date) -> TimeDelta:
        """
        Get daylight duration at given date.
        """
        sunrise, sunset = astral.sun.daylight(self._observer, date)
        return sunset - sunrise

    def solar_time(self, datetime: DateTime) -> Time:
        """
        Get solar time at given date and local standard time at the location.
        """
        datetime = self._localize(datetime)
        return solar_time(datetime, self.lon.to('deg').m)

    def hour_angle(self, datetime: DateTime) -> Quantity:
        """
        Get solar hour angle at given date and local standard time at the location.
        """
        datetime = self._localize(datetime)
        hra = hour_angle(datetime, self.lon.to('deg').m)
        return Q_(hra, 'deg')

    def sunset_hour_angle(self, date: Date) -> Quantity:
        """
        Get solar hour angle at sunset on the given date.
        """
        sunset = self.sunset(date)
        hra_sunset = hour_angle(sunset, self.lon.to('deg').m)
        return Q_(hra_sunset, 'deg')
