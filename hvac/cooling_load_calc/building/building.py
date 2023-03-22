from typing import TYPE_CHECKING
import pandas as pd
from hvac import Quantity

if TYPE_CHECKING:
    from .building_entity import BuildingEntity


class Building:

    def __init__(self):
        self.ID: str = ''
        self.building_entities: dict[str, BuildingEntity] = {}
        self.heat_gains: pd.DataFrame | None = None

    @classmethod
    def create(cls, ID: str) -> 'Building':
        obj = cls()
        obj.ID = ID
        return obj

    def add_building_entity(self, be: 'BuildingEntity'):
        self.building_entities[be.ID] = be
        be.building = self

    def get_heat_gains(self, unit: str = 'W') -> pd.DataFrame:
        self.heat_gains = sum(
            be.get_heat_gains(unit)
            for be in self.building_entities.values()
        )
        return self.heat_gains

    @property
    def floor_area(self) -> Quantity:
        A_floor = sum(
            be.floor_area
            for be in self.building_entities.values()
        )
        return A_floor

    @property
    def envelope_area(self) -> Quantity:
        A_env = sum(
            be.envelope_area
            for be in self.building_entities.values()
        )
        return A_env

    @property
    def volume(self) -> Quantity:
        V = sum(
            be.volume
            for be in self.building_entities.values()
        )
        return V
