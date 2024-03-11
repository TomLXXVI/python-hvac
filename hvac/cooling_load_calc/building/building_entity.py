from __future__ import annotations
from typing import TYPE_CHECKING
import pandas as pd
from hvac import Quantity

if TYPE_CHECKING:
    from .building import Building
    from .ventilation_zone import VentilationZone


class BuildingEntity:
    """Represents a single building entity, e.g. an apartment, in a larger
    building.
    """
    def __init__(self):
        self.ID: str = ''
        self.building: Building | None = None
        self.ventilation_zones: dict[str, VentilationZone] = {}
        self._is_solved: bool = False
        self._cooling_load_table: pd.DataFrame | None = None

    @classmethod
    def create(cls, ID: str) -> BuildingEntity:
        """Creates a `BuildingEntity` object.

        Parameters
        ----------
        ID:
            Identifies the building entity in the building.
        """
        obj = cls()
        obj.ID = ID
        return obj

    def add_ventilation_zone(self, *vez: VentilationZone) -> None:
        """Adds a ventilation zone (or multiple zones) to the building entity."""
        for vez_ in vez: vez_.building_entity = self
        self.ventilation_zones.update({vez_.ID: vez_ for vez_ in vez})

    def get_cooling_load_table(
        self,
        dt_hr: float = 1.0,
        num_cycles: int = 5,
        unit: str = 'W',
        do_recalculation: bool = False
    ) -> pd.DataFrame:
        """Returns the total cooling load of the building entity.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step for the calculations is one quarter of an hour.
            The default value is 1 hour.
        num_cycles:
            Number of diurnal calculation cycles before the results of the last
            diurnal cycle are returned.
        unit:
            The measuring unit in which heat flows need to be expressed. The
            default unit is Watts (W).
        do_recalculation:
            If True, indicates that a full recalculation must be performed of
            the heat gains and cooling loads of all thermal zones in the
            building entity, e.g. after modifications to thermal zones in the
            building entity or after modifications to ventilation zones in
            the building entity.
            By default, if the cooling loads of the thermal zones were already
            calculated once, the original table is returned.

        Returns
        -------
        A Pandas DataFrame with the heat gains and the total cooling load
        of the building entity for each time index k of the design day.
        """
        if not self._is_solved or do_recalculation:
            self._cooling_load_table = sum(
                vez.get_cooling_load_table(dt_hr, num_cycles, unit, do_recalculation)
                for vez in self.ventilation_zones.values()
            )
            self._is_solved = True
        return self._cooling_load_table

    @property
    def floor_area(self) -> Quantity:
        """Returns the total floor area of the building entity."""
        A_floor = sum(
            vz.floor_area
            for vz in self.ventilation_zones.values()
        )
        return A_floor

    @property
    def envelope_area(self) -> Quantity:
        """Returns the total exterior building envelope surface area of the
        building entity.
        """
        A_env = sum(
            vz.envelope_area
            for vz in self.ventilation_zones.values()
        )
        return A_env

    @property
    def volume(self) -> Quantity:
        """Returns the total volume of the building entity."""
        V = sum(
            vz.volume
            for vz in self.ventilation_zones.values()
        )
        return V
