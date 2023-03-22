from typing import TYPE_CHECKING
import pandas as pd
from hvac import Quantity

if TYPE_CHECKING:
    from .building import Building
    from .ventilation_zone import VentilationZone


class BuildingEntity:

    def __init__(self):
        self.ID: str = ''
        self.building: Building | None = None
        self.ventilation_zones: dict[str, VentilationZone] = {}
        self.heat_gains: pd.DataFrame | None = None

    @classmethod
    def create(cls, ID: str) -> 'BuildingEntity':
        obj = cls()
        obj.ID = ID
        return obj

    def add_ventilation_zone(self, ventilation_zone: 'VentilationZone'):
        self.ventilation_zones[ventilation_zone.ID] = ventilation_zone
        ventilation_zone.building_entity = self

    def get_heat_gains(self, unit: str = 'W') -> pd.DataFrame:
        self.heat_gains = sum(
            vz.get_heat_gains(unit)
            for vz in self.ventilation_zones.values()
        )
        return self.heat_gains

    @property
    def floor_area(self) -> Quantity:
        A_floor = sum(
            vz.floor_area
            for vz in self.ventilation_zones.values()
        )
        return A_floor

    @property
    def envelope_area(self) -> Quantity:
        A_env = sum(
            vz.envelope_area
            for vz in self.ventilation_zones.values()
        )
        return A_env

    @property
    def volume(self) -> Quantity:
        V = sum(
            vz.volume
            for vz in self.ventilation_zones.values()
        )
        return V
