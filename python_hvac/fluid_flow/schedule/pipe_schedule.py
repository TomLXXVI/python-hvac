from typing import Optional, Dict
import pandas as pd
from hvac import Quantity


class PipeSchedule:
    """Base class for standardized pipe schedules."""

    def __init__(self, unit: str = 'mm'):
        self.lookup_table: Optional[pd.DataFrame] = None
        self.unit = unit

    def _lookup(self, nominal_diameter: Quantity, item: str) -> Quantity:
        D_nom = int(nominal_diameter.to(self.unit).magnitude)
        # D_nom = nominal_diameter.to(self.unit).magnitude
        try:
            result = self.lookup_table.loc[D_nom, item]
        except KeyError:
            raise LookupError(f'nominal diameter {D_nom} not found')
        return Quantity(result, self.unit)

    def get_external_diameter(self, nominal_diameter: Quantity) -> Quantity:
        """Get the external diameter that corresponds with the given nominal diameter."""
        return self._lookup(nominal_diameter, 'D_ext')

    def get_wall_thickness(self, nominal_diameter: Quantity) -> Quantity:
        """Get the pipe wall thickness that corresponds with the given nominal diameter."""
        return self._lookup(nominal_diameter, 't')

    def get_internal_diameter(self, nominal_diameter: Quantity) -> Quantity:
        """Get the internal diameter that corresponds with the given nominal diameter."""
        return self._lookup(nominal_diameter, 'D_int')

    def get_closest_nominal_diameter(self, internal_diameter: Quantity) -> Quantity:
        """Get the nominal diameter of which the corresponding internal diameter is closest
        to the given internal diameter."""
        D_int = internal_diameter.to(self.unit).magnitude
        delta = [abs(D_int - D_int_lookup) for D_int_lookup in self.lookup_table['D_int']]
        i = delta.index(min(delta))
        D_nom = self.lookup_table.index[i]
        return Quantity(D_nom, self.unit)


pipe_schedule_40 = PipeSchedule()
pipe_schedule_40.lookup_table = pd.DataFrame({
    'D_nom': [6, 8, 10, 15, 20, 25, 32, 40, 50, 65, 80, 90, 100],                               # mm
    'D_ext': [10.3, 13.7, 17.1, 21.3, 26.7, 33.4, 42.2, 48.3, 60.3, 73.0, 88.9, 101.6, 114.3],  # mm
    't': [1.73, 2.24, 2.31, 2.77, 2.87, 3.38, 3.56, 3.68, 3.91, 5.16, 5.49, 5.74, 6.02],        # mm
    'D_int': [6.84, 9.22, 12.5, 15.8, 21.0, 26.6, 35.1, 40.9, 52.5, 62.7, 77.9, 90.1, 102.3]    # mm
})
pipe_schedule_40.lookup_table.set_index('D_nom', inplace=True)


class PipeScheduleFactory:
    schedules: Dict[str, PipeSchedule] = {
        'default': pipe_schedule_40
    }

    @classmethod
    def get(
            cls,
            name: str,
            file_path: Optional[str] = None,
            unit: str = 'mm'
    ) -> PipeSchedule:
        if name not in cls.schedules and file_path is not None:
            df = pd.read_excel(
                io=file_path,
                sheet_name=0,
                header=0,
                index_col=0
            )
            pipe_schedule = PipeSchedule(unit)
            pipe_schedule.lookup_table = df
            cls.schedules[name] = pipe_schedule
        return cls.schedules.get(name)
