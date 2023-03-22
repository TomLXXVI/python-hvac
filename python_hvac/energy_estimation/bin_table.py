from typing import Optional, Dict, Union, Tuple, List
from enum import IntEnum
from datetime import datetime
from dataclasses import dataclass
import pandas as pd
import numpy as np
from hvac import Quantity


Q_ = Quantity


class Month(IntEnum):
    Jan = 1
    Feb = 2
    Mar = 3
    Apr = 4
    May = 5
    Jun = 6
    Jul = 7
    Aug = 8
    Sep = 9
    Oct = 10
    Nov = 11
    Dec = 12


@dataclass
class TimeSegment:
    """
    Helper class to add time segments to a bin table.

    Parameters
    ----------
    label: str
        A label to identify the time segment.
    limits: Tuple[int, int]
        The lower and upper hour of the time segment.
    """
    label: str
    limits: Tuple[int, int]


class BinTableCreator:
    """Create monthly temperature bin tables from TMY-file (csv) containing
    the hourly average outside air temperature of each hour of every day of the
    year."""
    _default_time_segments: Dict[str, TimeSegment] = {
        '0 - 3 h': TimeSegment('0 - 3 h', (0, 4)),
        '4 - 7 h': TimeSegment('4 - 7 h', (4, 8)),
        '8 - 11 h': TimeSegment('8 - 11 h', (8, 12)),
        '12 - 15 h': TimeSegment('12 - 15 h', (12, 16)),
        '16 - 19 h': TimeSegment('16 - 19 h', (16, 20)),
        '20 - 23 h': TimeSegment('20 - 23 h', (20, 24))
    }

    def __init__(
        self,
        file_path: str,
        date_time_column_name: str,
        temperature_column_name: str,
        bin_limits: Tuple[Quantity, Quantity] = (Q_(-10, 'degC'), Q_(40, 'degC')),
        bin_width: Quantity = Q_(2, 'K'),
        date_time_fmt: str = '%Y%m%d:%H%M',
        time_segments: Optional[List[TimeSegment]] = None,
        T_unit: str = 'degC',
    ) -> None:
        """
        Creates a dictionary of monthly temperature bin tables from a csv-file
        with TMY hourly averaged outdoor air temperatures. After instantiation,
        use method `get_bin_table` to get the bin table of a given month. A bin
        table is represented as a Pandas DataFrame object. Each temperature bin
        contains the number of hours the hourly averaged outdoor air temperature
        resides in that temperature bin.

        Parameters
        ----------
        file_path: str
            Path of the csv-file with TMY data.
        date_time_column_name: str
            The column title name in which date-time values are present.
        temperature_column_name: str
            The column title name in which temperature values are present.
        bin_limits: Tuple[Quantity, Quantity], default (-10 °C, 40 °C)
            The lower and upper limit of the temperature bins.
        bin_width: Quantity, default 2 K
            The width of one temperature bin in units of temperature difference.
        date_time_fmt: str, default '%Y%m%d:%H%M'
            The date-time format in which the date-time values are expressed
        time_segments: List[TimeSegment], optional
            Time segments of the day for which the hours the outdoor air
            temperature is in a certain temperature bin are counted. If not
            specified or None, the default time segments of the class will
            be used: 0 to 3 h, 4 to 7 h, 8 to 11 h, 12 to 15 h, 16 to 19 h and
            20 to 23 h.
        T_unit: str, default 'degC'
            The measurement unit of the temperature values in the TMY datafile.

        Notes
        -----
        The temperature measurement unit used to specify the bin limits is
        considered to be the unit desired by the user. If `T_unit` is different
        from the desired unit, the temperature values in the TMY datafile are
        converted from `T_unit` to the desired unit.
        """
        self._file_path = file_path
        self._date_time_column_name = date_time_column_name
        self._temperature_column_name = temperature_column_name
        self._bin_limits = bin_limits
        self._bin_width = bin_width
        self._date_time_fmt = date_time_fmt
        self._T_unit = T_unit
        if time_segments is not None:
            self.time_segments = {
                time_segment.label: time_segment
                for time_segment in time_segments
            }
        else:
            self.time_segments = self._default_time_segments
        self._yearly_tmy_data: Optional[pd.DataFrame] = None
        self._monthly_tmy_data: Optional[pd.DataFrameGroupBy] = None
        self._monthly_bin_tables: Dict[int, pd.DataFrame] = {}

        self._create_bin_tables()

    def get_bin_table(self, month: Union[Month, int]) -> pd.DataFrame:
        """Get the temperature bin table of the specified month."""
        return self._monthly_bin_tables[month]

    def _create_bin_tables(self) -> None:
        """Reads the TMY-data file and creates the temperature bin table for
         each month of the year."""
        self._read_file()
        self._group_by_month()
        for m in range(1, 13):
            self._monthly_bin_tables[m] = self._create_monthly_bin_table(m)

    def _read_file(self) -> None:
        """Read TMY data from csv-file."""
        self._yearly_tmy_data = pd.read_csv(
            self._file_path,
            parse_dates=[0],
            infer_datetime_format=True,
            dayfirst=True,
            date_parser=lambda dt_str: datetime.strptime(dt_str, self._date_time_fmt)
        )
        # convert temperature values if wanted
        T_unit_desired = self._bin_limits[0].units
        if T_unit_desired != Q_(0, self._T_unit).units:
            T_column = self._yearly_tmy_data.loc[:, self._temperature_column_name]
            T_column_converted = [Q_(T, self._T_unit).to(T_unit_desired).m for T in T_column]
            self._yearly_tmy_data[self._temperature_column_name] = T_column_converted

    def _group_by_month(self) -> None:
        """Split yearly TMY data into monthly groups."""
        self._monthly_tmy_data = self._yearly_tmy_data.groupby(
            self._yearly_tmy_data[self._date_time_column_name].dt.month
        )

    def _create_monthly_bin_table(self, month_index: int) -> pd.DataFrame:
        """Create a temperature bin table for the month specified by its
        month of the year index (1...12)."""

        # get dataframe of month specified by `month_index`
        month = self._monthly_tmy_data.get_group(month_index)

        # group the monthly tmy data by hour of the day (ignoring the dates)
        hour_groups = month.groupby(month[self._date_time_column_name].dt.hour)

        # arrange hour groups into time segments
        time_segments = self._get_time_segment_dict(hour_groups)

        # create a list with the boundaries of the temperature bins
        T_start = self._bin_limits[0].m
        T_end = self._bin_limits[1].m
        dT_bin = self._bin_width.m
        T_bins = np.arange(T_start, T_end, dT_bin)

        # for each time segment put the temperature values into bins and count
        # how many temperature values are in each bin; put the bin series of
        # each time segment in a list
        bin_table = []
        for key, time_segment in time_segments.items():
            bin_series = pd.value_counts(
                pd.cut(x=time_segment[self._temperature_column_name], bins=T_bins)
            )
            bin_series.name = key
            bin_table.append(bin_series)

        # merge the list of bin series of each hourly segment into a single
        # dataframe and sort the bins from low to high
        bin_table = pd.concat(bin_table, axis=1)
        bin_table = bin_table.sort_index()

        return bin_table

    def _get_time_segment_dict(self, hour_groups) -> Dict[str, pd.DataFrame]:
        """Arranges the hour groups into time segments."""

        # Initialize a dictionary of empty lists of which the keys refer to a
        # time segment
        time_segment_dict = {}
        for time_segment in self.time_segments.values():
            time_segment_dict[time_segment.label] = []

        # Append each hour group to the correct time segment list
        for hour, hour_group in hour_groups:
            for time_segment in self.time_segments.values():
                if hour in range(time_segment.limits[0], time_segment.limits[1]):
                    time_segment_dict[time_segment.label].append(hour_group)
                    break

        # Concatenate the hour groups in each list into a single dataframe
        # for each hourly segment
        for key in time_segment_dict.keys():
            time_segment_dict[key] = pd.concat(time_segment_dict[key])

        return time_segment_dict
