"""Module for reading the TMY-data from a csv-file into a Pandas-DataFrame
object.
"""
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.climate.sun.solar_time import solar_time

Q_ = Quantity


class TMY:
    """Class for reading TMY datafiles."""

    def __init__(
        self,
        file_path: str,
        dt_fmt: str = '%Y%m%d:%H%M',
        tz: str = 'UTC',
    ) -> None:
        """
        Creates a `TMY` instance by reading the csv-file pointed to by
        `file_path` into a Pandas DataFrame object.

        The implementation of this class is based on the csv-file with TMY-data
        that is downloaded from PVGIS.
        (https://re.jrc.ec.europa.eu/pvg_tools/en/#api_5.2)
        This means that date and time are grouped in a single column, which is
        the first column of the csv-file. Also, the cvs-file must only contain
        the TMY-data (other rows or columns need to be removed from the
        csv-file).

        Parameters
        ----------
        file_path:
            The path to the csv-file containing the TMY-data.
        dt_fmt: default '%Y%m%d:%H%M'
            The format in which date-times are written in the csv-file.
        tz: default 'UTC'
            The timezone of the date-times in the csv-file, expressed in
            its tz-database notation.
            (https://en.wikipedia.org/wiki/List_of_tz_database_time_zones).
        """
        self.file_path = file_path
        self.dataframe = self._read_csv_file(dt_fmt, tz)

    def _read_csv_file(self, dt_fmt: str, tz: str) -> pd.DataFrame:
        df = pd.read_csv(
            self.file_path,
            parse_dates=[0],
            dayfirst=True,
            date_format=dt_fmt
        )
        df = df.set_index(df.columns[0])
        df.index = df.index.tz_localize(tz, nonexistent='shift_forward')
        return df

    @property
    def header(self) -> pd.Index:
        """Returns the column names of the csv-file."""
        return self.dataframe.columns

    def get_month(self, month: int) -> pd.DataFrame:
        """Returns only the TMY-data for the given month, indicated by its
        integer index.
        """
        months = self.dataframe.groupby(self.dataframe.index.month)
        month = months.get_group(month)
        return month

    def get_day(self, month: int, day: int) -> pd.DataFrame:
        """Returns only the TMY-data for the given month and day, both indicated
        by their integer index.
        """
        month = self.get_month(month)
        days = month.groupby(month.index.day)
        day = days.get_group(day)
        return day

    def convert_to_timezone(self, tz: str) -> None:
        """Convert the date-time index to time zone `tz`."""
        self.dataframe.index = self.dataframe.index.tz_convert(tz)
        self.dataframe.index.name = f'Time ({tz})'
        self.dataframe.sort_index(inplace=True)

    def convert_to_solar_time(self, tz: str, L_loc: float) -> None:
        """Convert the date-time index to local solar time. Parameter `tz` is
        the local timezone with respect to GMT or UTC (e.g. 'Etc/GMT-1') and is
        used to determine the local standard time meridian of the time zone.
        Do not specify a timezone with DST, should it be in use at the location
        under consideration, as this may give errors at times when DST becomes
        active or inactive.
        Longitude `L_loc` in decimal degrees must be positive east of UTC and
        negative west of UTC.
        """
        if self.dataframe.index.tz != tz:
            self.convert_to_timezone(tz)
        self.dataframe.index = pd.DatetimeIndex([
            solar_time(dt_local, L_loc)
            for dt_local in self.dataframe.index
        ])
        self.dataframe.index.name = 'Time (solar time)'


def frequency_table(
    dataset: np.ndarray,
    bins: np.ndarray | str | int = 'auto'
) -> pd.DataFrame:
    """Returns a Pandas-DataFrame with the frequency table of the given dataset.

    Parameters
    ----------
    dataset:
        Numpy array of values.
    bins:
        See Numpy documentation about `numpy.histogram`:
        If bins is an int, it defines the number of equal-width bins in the
        given range. If bins is a sequence, it defines a monotonically
        increasing array of bin edges, including the rightmost
        edge, allowing for non-uniform bin widths. If bins is a string, it
        defines the method used to calculate the optimal bin width, as defined
        by `histogram_bin_edges`. By default, `bins` is set to 'auto', which
        means that Numpy will determine the optimal bin width and consequently
        the number of bins.
    """
    freq, bin_edges = np.histogram(dataset, bins=bins)
    bin_edges = bin_edges.round(2)
    bins = [
        f"[{bin_edges[i]}, {bin_edges[i+1]})"
        for i in range(len(bin_edges) - 2)
    ]
    bins += [f"[{bin_edges[-2]}, {bin_edges[-1]}]"]
    cum_freq = np.cumsum(freq)
    df = pd.DataFrame({
        'bin': bins,
        'mid': [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(bin_edges) - 1)],
        'freq': freq,
        'cum_freq': cum_freq
    })
    return df
