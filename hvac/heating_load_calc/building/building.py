from __future__ import annotations
from dataclasses import dataclass
import pandas as pd
from hvac import Quantity

from . import BuildingEntity

Q_ = Quantity


@dataclass
class ClimateDesignData:
    """Holds the climatic design information needed to perform a heating loss
    calculation (for Belgium: see NBN EN 12831-1 ANB, NA.2 - Table NA.1).

    Parameters
    ----------
    T_ext_d: Quantity
        External design temperature for the considered building.
    T_ext_an: Quantity
        Annual mean external temperature.
    T_ext_min: Quantity
        Average minimal external temperature during the coldest month.
    """
    T_ext_d: Quantity
    T_ext_an: Quantity
    T_ext_min: Quantity


class Building:

    def __init__(self):
        self.ID: str = ''
        self.climate_data: ClimateDesignData | None = None
        self.building_entities: dict[str, BuildingEntity] = {}

    @classmethod
    def create(cls, ID: str, climate_data: ClimateDesignData) -> Building:
        """Creates a new building.

        Parameters
        ----------
        ID: str
            Name of the building.
        climate_data: ClimateDesignData
            Instance of `ClimateDesignData` with the climate data at the
            geographical location of the building needed to a perform a
            heat load calculation.

        Returns
        -------
        Instance of class `Building`.
        """
        self = cls()
        self.ID = ID
        self.climate_data = climate_data
        return self

    def add_building_entity(self, ID: str) -> BuildingEntity:
        """Adds a new building entity to the building.

        Parameters
        ----------
        ID: str
            Name of the building entity.

        Returns
        -------
        Instance of class `BuildingEntity`.
        """
        be = BuildingEntity.create(ID, self.climate_data.T_ext_d)
        self.building_entities[be.ID] = be
        return be

    def get_transmission_heat_loss(self) -> Quantity:
        """Returns the heat loss of the building due to heat conduction through
        the building elements.
        """
        Q_trm = sum(
            be.get_transmission_heat_loss()
            for be in self.building_entities.values()
        )
        return Q_trm

    def get_ventilation_heat_loss(self) -> Quantity:
        """Returns the heat loss of the building due to space ventilation and
        outdoor air infiltration.
        """
        Q_ven = sum(
            vz.get_ventilation_heat_loss()
            for be in self.building_entities.values()
            for vz in be.ventilation_zones.values()
        )
        return Q_ven

    def get_additional_heating_up_power(self) -> Quantity:
        """Returns the total heat rate needed to heat-up the spaces in the
        building after a temperature set-back period.
        """
        Q_hu = sum(
            be.get_additional_heating_up_power()
            for be in self.building_entities.values()
        )
        return Q_hu

    def get_heat_load(self) -> Quantity:
        """Returns the total heating load of the building."""
        Q_trm = self.get_transmission_heat_loss()
        Q_ven = self.get_ventilation_heat_loss()
        Q_hu = self.get_additional_heating_up_power()
        return Q_trm + Q_ven + Q_hu

    def get_summary(self, unit: str = 'kW', n_digits: int = 3) -> pd.DataFrame:
        """Returns a Pandas Dataframe with an overview of the building entities
        in the building.

        Parameters
        ----------
        unit: str
            Desired unit of thermal power.
        n_digits: int
            The number of decimals displayed in the numeric results.
        """
        col_1 = 'building entity'
        col_2 = f'Q transmission [{unit}]'
        col_3 = f'Q ventilation [{unit}]'
        col_4 = f'Q heating-up [{unit}]'
        col_5 = f'Q total [{unit}]'
        d = {
            col_1: [],
            col_2: [],
            col_3: [],
            col_4: [],
            col_5: []
        }
        for be in self.building_entities.values():
            d[col_1].append(be.ID)
            d[col_2].append(round(be.get_transmission_heat_loss().to(unit).m, n_digits))
            d[col_3].append(round(be.get_ventilation_heat_loss().to(unit).m, n_digits))
            d[col_4].append(round(be.get_additional_heating_up_power().to(unit).m, n_digits))
            d[col_5].append(round(be.get_heat_load().to(unit).m, n_digits))
        df = pd.DataFrame(d)
        return df
