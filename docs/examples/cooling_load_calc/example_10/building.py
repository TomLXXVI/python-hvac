"""
COOLING LOAD CALCULATION OF A TWO-ZONE BUILDING

This example demonstrates how the creation of the model of a two-zone building
can be organized to perform a cooling load calculation of this building.
Each zone of the building is programmed in a separate module. This module then
imports these "zone-modules" to create the model of the whole building.

A building model has the hierarchy of a tree. A `Building` instance contains
at least one `BuildingEntity` instance. A `BuildingEntity` instance contains
at least one `VentilationZone` instance. And a `VentilationZone` instance
contains at least one `FixedTemperatureZone` instance.

In module `zone_1.py`, the space/thermal zone "Zone 1" of the two-zone building
is created and configured.
In module `zone_2.py`, the space/thermal zone "Zone 2" of the two-zone building
is created and configured.
In this module, the two-zone building is created and configured.

The Jupyter Notebook `cooling_load_calc.ipynb` instantiates the
`TwoZoneBuilding` class in this module and runs the cooling load calculation
by calling the method `cooling_load()` on this instance.
"""
import pandas as pd
from hvac import Quantity
from hvac.charts import LineChart
from hvac.cooling_load_calc import (
    WeatherData,
    Building,
    BuildingEntity,
    VentilationZone
)
import zone_1
import zone_2


Q_ = Quantity


class TwoZoneBuilding:

    def __init__(
        self,
        weather_data: WeatherData,
        T_int_des: Quantity
    ) -> None:
        self.weather_data = weather_data
        self.T_int_des = T_int_des

        # Creating the building hierarchy tree.
        # Create the `Building` object, the `BuildingEntity` object (the
        # building has only one building entity that coincides with the
        # building), and the `VentilationZone` object (the building has only one
        # ventilation zone that coincides with the building entity).
        self._building = Building.create('two-zone-building')
        self._building_entity = BuildingEntity.create('building-entity')
        self._ventilation_zone = VentilationZone.create(
            ID='ventilation-zone',
            q_env_50=Q_(6, 'm**3 / (hr * m**2)'),
            f_iz=1.0
        )
        # Link the ventilation zone to the building entity:
        self._building_entity.add_ventilation_zone(self._ventilation_zone)
        # Link the building entity with the building.
        self._building.add_building_entity(self._building_entity)

        # Instantiate zone 1 of the building (also passing the ventilation zone
        # to which zone 1 belongs):
        self.zone_1 = zone_1.Zone1(
            self.weather_data,
            self.T_int_des,
            self._ventilation_zone
        )
        # Instantiate zone 2 of the building (also passing the ventilation zone
        # to which zone 2 belongs):
        self.zone_2 = zone_2.Zone2(
            self.weather_data,
            self.T_int_des,
            self._ventilation_zone
        )

    def cooling_load(self, unit: str = 'W') -> pd.DataFrame:
        """Returns the hourly heat gains and cooling load of the whole building
        on the design day."""
        df = self._building.get_cooling_load_table(unit=unit)
        return df

    @staticmethod
    def plot_cooling_load(df: pd.DataFrame) -> LineChart:
        """Returns a `LineChart` object with the total cooling load, the
        sensible cooling load, and the latent cooling load of the whole
        building. The loads are plotted along the x-axis in solar time.
        """
        lch = LineChart()
        columns = {
            'Q_dot_sen_zone': 'sensible cooling load',
            'Q_dot_lat_zone': 'latent cooling load',
            'Q_dot_zone': 'total cooling load'
        }
        for k, v in columns.items():
            lch.add_xy_data(
                label=v,
                x1_values=df.index[:-1],
                y1_values=df[k][:-1]
            )
        lch.x1.add_title('solar time index')
        lch.y1.add_title('thermal power, W')
        lch.add_legend(columns=3)
        return lch
