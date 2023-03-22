from typing import List, Dict, Optional
import pandas as pd
from hvac import Quantity


class DuctSchedule:
    """Base class for standardized duct schedules."""

    def __init__(self):
        self.lookup_list: List[float] = []

    def get_closest_internal_size(self, calculated_size: Quantity) -> Quantity:
        """Get the internal size that is closest to the calculated size."""
        calculated_size = calculated_size.to('mm').magnitude
        delta = [abs(calculated_size - lookup_size) for lookup_size in self.lookup_list]
        index = delta.index(min(delta))
        available_size = self.lookup_list[index]
        return Quantity(available_size, 'mm')


circular_duct_schedule = DuctSchedule()
circular_duct_schedule.lookup_list = [
    50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0,
    250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,
    1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0
]  # mm

rectangular_duct_schedule = DuctSchedule()
rectangular_duct_schedule.lookup_list = [i for i in range(100, 300, 25)]
rectangular_duct_schedule.lookup_list += [i for i in range(300, 800, 50)]
rectangular_duct_schedule.lookup_list += [i for i in range(800, 3000, 100)]  # mm

flat_oval_duct_schedule = DuctSchedule()
flat_oval_duct_schedule.lookup_list = [i for i in range(175, 1500, 25)]
flat_oval_duct_schedule.lookup_list += [i for i in range(1500, 1850, 50)]  # mm


class DuctScheduleFactory:
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
    ) -> DuctSchedule:
        if name not in cls.schedules and file_path is not None:
            df = pd.read_excel(
                io=file_path,
                sheet_name=0,
                header=0
            )
            duct_schedule = DuctSchedule()
            duct_schedule.lookup_list = df
            cls.schedules[name] = duct_schedule
        return cls.schedules.get(name)
