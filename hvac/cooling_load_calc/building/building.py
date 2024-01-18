from __future__ import annotations
from typing import TYPE_CHECKING
import pandas as pd
from hvac import Quantity

if TYPE_CHECKING:
    from .building_entity import BuildingEntity


class Building:

    def __init__(self):
        self.ID: str = ''
        self.building_entities: dict[str, BuildingEntity] = {}
        self._is_solved: bool = False
        self._cooling_load_table: pd.DataFrame | None = None

    @classmethod
    def create(cls, ID: str) -> Building:
        """Creates a `BuildingEntity` object.

        Parameters
        ----------
        ID:
            Identifies the building entity in the building.
        """
        obj = cls()
        obj.ID = ID
        return obj

    def add_building_entity(self, be: BuildingEntity):
        """Adds a building entity to the building."""
        self.building_entities[be.ID] = be
        be.building = self

    def get_cooling_load_table(
        self,
        dt_hr: float = 1.0,
        num_cycles: int = 5,
        unit: str = 'W',
        do_recalculation: bool = False
    ) -> pd.DataFrame:
        """Returns the total cooling load of the building.

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
            building, e.g. after modifications to thermal zones in the
            building or after modifications to the ventilation zones.
            By default, if the cooling loads of the thermal zones were already
            calculated once, the original table is returned.

        Returns
        -------
        A Pandas DataFrame with the heat gains and the total cooling load
        of the building for each time index k of the design day.
        """
        if not self._is_solved or do_recalculation:
            self._cooling_load_table = sum(
                be.get_cooling_load_table(dt_hr, num_cycles, unit, do_recalculation)
                for be in self.building_entities.values()
            )
            self._is_solved = True
        return self._cooling_load_table

    @property
    def floor_area(self) -> Quantity:
        """Returns the total floor area of the building."""
        A_floor = sum(
            be.floor_area
            for be in self.building_entities.values()
        )
        return A_floor

    @property
    def envelope_area(self) -> Quantity:
        """Returns the total exterior building envelope surface area of the
        building.
        """
        A_env = sum(
            be.envelope_area
            for be in self.building_entities.values()
        )
        return A_env

    @property
    def volume(self) -> Quantity:
        """Returns the total volume of the building."""
        V = sum(
            be.volume
            for be in self.building_entities.values()
        )
        return V
