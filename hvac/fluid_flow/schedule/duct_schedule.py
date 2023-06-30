from typing import List, Dict, Optional, Union
import pandas as pd
from hvac import Quantity


class DuctSchedule:
    """Base class for schedules of commercially available duct sizes."""

    def __init__(self, unit: str = 'mm'):
        """Creates a `DuctSchedule` instance.

        Parameters
        ----------
        unit:
            The length unit in which the sizes of the schedule are expressed.
        """
        self.lookup_list: List[float] = []
        self.unit = unit

    def get_closest_internal_size(self, calculated_size: Quantity) -> Quantity:
        """Get the internal size that is closest to the calculated size."""
        calculated_size = calculated_size.to(self.unit).magnitude
        delta = [
            abs(calculated_size - lookup_size)
            for lookup_size in self.lookup_list
        ]
        index = delta.index(min(delta))
        available_size = self.lookup_list[index]
        return Quantity(available_size, self.unit)


# default schedule for circular ducts
circular_duct_schedule = DuctSchedule()
circular_duct_schedule.lookup_list = [
    50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0,
    250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,
    1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0
]  # mm

# default schedule for rectangular ducts
rectangular_duct_schedule = DuctSchedule()
rectangular_duct_schedule.lookup_list = [i for i in range(100, 300, 25)]
rectangular_duct_schedule.lookup_list += [i for i in range(300, 800, 50)]
rectangular_duct_schedule.lookup_list += [i for i in range(800, 3000, 100)]  # mm

# default schedule for flat oval ducts
flat_oval_duct_schedule = DuctSchedule()
flat_oval_duct_schedule.lookup_list = [i for i in range(175, 1500, 25)]
flat_oval_duct_schedule.lookup_list += [i for i in range(1500, 1850, 50)]  # mm


class DuctScheduleFactory:
    """Class for storing duct schedules at run-time and read user defined duct
    schedules from file.
    """
    schedules: Dict[str, DuctSchedule] = {
        'default_circular': circular_duct_schedule,
        'default_rectangular': rectangular_duct_schedule,
        'default_flat_oval': flat_oval_duct_schedule
    }

    @classmethod
    def get(
        cls,
        name: str,
        file_path: Optional[str] = None,
        sheet_name: Union[str, int] = 0,
        unit: str = 'mm'
    ) -> DuctSchedule:
        """Get a schedule from a file or one of the default schedules.

        Parameters
        ----------
        name:
            Identifier of the schedule. The names of the default schedules
            already present are:
            - 'default_circular' for circular ducts
            - 'default_rectangular' for rectangular ducts
            - 'default_flat_oval' for flat oval ducts
        file_path: optional, default None
            Path to a spreadsheet file with the commercially available sizes of
            a duct schedule.
        sheet_name:
            The name or index of the sheet.
        unit:
            The measuring unit in which the values in the spreadsheet are
            expressed.

        Returns
        -------
        DuctSchedule

        Notes
        -----
        The spreadsheet must only contain a single list of values (i.e.
        diameters or lengths of a side/axis of the duct cross-sections) in the
        first column of the sheet, starting at the first row (so don't use a
        header).

        In case of rectangular or flat oval ducts only the length of one side
        or axis can be specified. When sizing the cross-section of a duct
        the program assumes that the height of the cross-section will always be
        known by the user and that the width has to be determined, which means
        the program will select the width that is closest to the calculated
        value from the given schedule.
        """
        if name not in cls.schedules and file_path is not None:
            df = pd.read_excel(
                io=file_path,
                sheet_name=sheet_name,
                header=None
            )
            duct_schedule = DuctSchedule(unit)
            duct_schedule.lookup_list = df.iloc[:, 0].tolist()
            cls.schedules[name] = duct_schedule
        return cls.schedules.get(name)
