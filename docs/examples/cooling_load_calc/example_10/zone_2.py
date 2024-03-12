"""
COOLING LOAD CALCULATION OF ZONE 2 IN THE TWO-ZONE BUILDING.

In this module the space/thermal zone "Zone 2" of the two-zone building is
created and configured. The same procedure as in module `zone_1.py` is followed.
"""
from typing import Callable
from datetime import time
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.charts import LineChart
from hvac.cooling_load_calc.core.utils import convert_to_solar_seconds
from hvac.cooling_load_calc import (
    WeatherData,
    wtcb,
    ExteriorBuildingElement,
    VentilationZone,
    FixedTemperatureZone,
    shelves,
    OfficeEquipment,
    EquipmentHeatGain,
    PeopleHeatGain,
    ExteriorShadingDevice,
    InteriorShadingDevice
)

Q_ = Quantity


class _SpaceConstruction:

    def __init__(
        self,
        weather_data: WeatherData,
        T_int_des: Quantity
    ) -> Quantity:
        self.weather_data = weather_data
        self.T_int_des = T_int_des

        self.exterior_walls = wtcb.ExteriorWallCatalog(
            t_ins=Q_(12, 'cm'),
            T_ext=self.weather_data.T_db_avg,
            T_int=self.T_int_des
        )
        self.roofs = wtcb.RoofCatalog(
            t_ins=Q_(18, 'cm'),
            T_ext=self.weather_data.T_db_avg,
            T_int=self.T_int_des
        )

        # Exterior space dimensions:
        self.space_width = Q_(4, 'm')
        self.space_length = Q_(10, 'm')
        self.space_height = Q_(4, 'm')

        # Declare exterior building elements of the space:
        self.south_wall: ExteriorBuildingElement | None = None
        self.east_wall: ExteriorBuildingElement | None = None
        self.north_wall: ExteriorBuildingElement | None = None
        self.east_roof: ExteriorBuildingElement | None = None

        self.zone: FixedTemperatureZone | None = None

    @staticmethod
    def _gable(w: Quantity, alpha_roof: Quantity) -> tuple[Quantity, Quantity]:
        """Calculates the surface area of the gable and the sloped length of
        one roof half.
        """
        alpha_roof = alpha_roof.to('rad')
        l_slope = 0.5 * w / np.cos(alpha_roof)
        h = l_slope * np.sin(alpha_roof)
        A_gable = (w * h) / 2
        return A_gable.to('m**2'), l_slope.to('m')

    def _interior_floor_area(self) -> Quantity:
        """Calculates the interior floor area of the space/thermal zone."""
        w = self.space_width
        w -= self.east_wall.constr_assem.thickness
        l = self.space_length
        l -= self.south_wall.constr_assem.thickness
        l -= self.north_wall.constr_assem.thickness
        A = w * l
        return A.to('m**2')

    def _create_south_wall(self) -> ExteriorBuildingElement:
        # Get the construction assembly from the WTCB catalog with the exterior
        # walls. To select the construction assembly use its ID in the WTCB
        # catalog.
        constr_assem = self.exterior_walls('F1')
        # Calculate the gross area of the wall:
        w = self.space_width
        h = self.space_height
        A_wall = w * h
        A_gable, _ = self._gable(w, alpha_roof=Q_(30, 'deg'))
        A_wall += A_gable
        # Create the exterior wall:
        wall = ExteriorBuildingElement.create(
            ID='south-wall',
            T_zone=lambda t: self.T_int_des,
            constr_assem=constr_assem,
            gross_area=A_wall,
            weather_data=self.weather_data,
            gamma=Q_(0, 'deg'),
            beta=Q_(90, 'deg'),
        )
        overhang = ExteriorShadingDevice(
            hor_proj_dist=Q_(0.75, 'm'),
            vert_offset=Q_(0.3, 'm')
        )
        wall.add_window(
            ID='window-1',
            width=Q_(1.5, 'm'),
            height=Q_(2, 'm'),
            props=shelves.WindowPropertiesShelf.load('window-5a-operable-wood/vinyl'),
            ext_shading=overhang
        )
        return wall

    def _create_east_wall(self) -> ExteriorBuildingElement:
        # Get the construction assembly:
        constr_assem = self.exterior_walls('F1')
        # Calculate the gross area of the wall:
        w = self.space_length
        h = self.space_height
        A_wall = w * h
        # Create the exterior wall:
        wall = ExteriorBuildingElement.create(
            ID='east-wall',
            T_zone=lambda t: self.T_int_des,
            constr_assem=constr_assem,
            gross_area=A_wall,
            weather_data=self.weather_data,
            gamma=Q_(270, 'deg'),
            beta=Q_(90, 'deg')
        )
        # Add window with a louvered shade on the inside (louver reflection
        # 0.5, slat angle 0°):
        louvered_shade = InteriorShadingDevice(
            IAC_dif=0.89,
            F_rad=0.71,
            IAC_0=0.98,
            IAC_60=0.80
        )
        wall.add_window(
            ID='window-1',
            width=Q_(2, 'm'),
            height=Q_(2, 'm'),
            props=shelves.WindowPropertiesShelf.load('window-5a-operable-wood/vinyl'),
            int_shading=louvered_shade
        )
        return wall

    def _create_north_wall(self) -> ExteriorBuildingElement:
        # Get the construction assembly:
        constr_assem = self.exterior_walls('F1')
        # Calculate the gross area of the wall:
        w = self.space_width
        h = self.space_height
        A_wall = w * h
        A_gable, _ = self._gable(w, alpha_roof=Q_(30, 'deg'))
        A_wall += A_gable
        # Create the exterior wall:
        wall = ExteriorBuildingElement.create(
            ID='north-wall',
            T_zone=lambda t: self.T_int_des,
            constr_assem=constr_assem,
            gross_area=A_wall,
            weather_data=self.weather_data,
            gamma=Q_(180, 'deg'),
            beta=Q_(90, 'deg'),
        )
        return wall

    def _create_east_roof(self) -> ExteriorBuildingElement:
        # Get the construction assembly:
        constr_assem = self.roofs('F1')
        # Calculate the gross area of the roof part:
        l = self.space_length
        _, w = self._gable(w=self.space_width, alpha_roof=Q_(30, 'deg'))
        A_roof = l * w
        # Create the roof part:
        roof = ExteriorBuildingElement.create(
            ID='east-roof',
            T_zone=lambda t: self.T_int_des,
            constr_assem=constr_assem,
            gross_area=A_roof,
            weather_data=self.weather_data,
            gamma=Q_(270, 'deg'),
            beta=Q_(30, 'deg')
        )
        # Add windows with dark opaque roller shades to the east-side roof half:
        roller_shade = InteriorShadingDevice(IAC_dif=0.76, F_rad=0.44)
        for i in range(2):
            roof.add_window(
                ID=f'window-{i+1}',
                width=Q_(1, 'm'),
                height=Q_(1, 'm'),
                props=shelves.WindowPropertiesShelf.load('window-5a-operable-wood/vinyl'),
                F_rad=0.33,
                int_shading=roller_shade
            )
        return roof

    def build(self) -> FixedTemperatureZone:
        """(1) Creates the building elements that surround the space/thermal
        zone, (2) creates the space/thermal zone (`FixedTemperatureZone`
        object), (3) adds the building elements to the space/thermal zone, (4)
        adds interior thermal mass to the space/thermal zone, (5) returns the
        space/thermal zone.
        """
        # Create the exterior building elements that surround the space:
        self.south_wall = self._create_south_wall()
        self.east_wall = self._create_east_wall()
        self.north_wall = self._create_north_wall()
        self.east_roof = self._create_east_roof()
        # Create the `FixedTemperatureZone` object that represents the space:
        A_flr = self._interior_floor_area()
        self.zone = FixedTemperatureZone.create(
            ID='zone-2',
            weather_data=self.weather_data,
            T_zone=lambda t: self.T_int_des,
            floor_area=A_flr,
            height=self.space_height
        )
        # Add the exterior building elements to the space:
        self.zone.add_ext_build_elem(
            self.south_wall, self.east_wall,
            self.north_wall, self.east_roof
        )
        # Add interior thermal mass to the space:
        # It is estimated that the floor area + half the wall area contributes
        # to the interior thermal mass.
        walls = [self.south_wall, self.north_wall, self.east_wall]
        A_walls = sum(w.net_area for w in walls)
        self.zone.add_thermal_storage_node(
            C=Q_(200, 'kJ / (K * m**2)'),
            A=A_flr + A_walls / 2,
            R_tz=Q_(0.1, 'K * m**2 / W')
        )
        return self.zone


class _SpaceVentilation:
    """
    Addition of local ventilation/air infiltration to the space/thermal zone.
    Calling method `configure` configures the local ventilation/air infiltration
    of the space/thermal zone.
    """
    def __init__(self, zone: FixedTemperatureZone, vez: VentilationZone) -> None:
        """Creates `_SpaceVentilation` instance.

        Parameters
        ----------
        zone:
            The space/thermal zone to which local ventilation/air infiltration
            is to be added.
        vez:
            The ventilation zone of the building to which the space/thermal zone
            belongs.
        """
        self.zone = zone
        # Link the building's ventilation zone to the space (zone):
        self.zone.ventilation_zone = vez
        # Link the space (zone) to the building's ventilation zone:
        self.zone.ventilation_zone.add_thermal_zone(self.zone)

    def configure(self) -> None:
        """Configures the local air infiltration/ventilation of the
        space/thermal zone.
        """
        self.zone.add_ventilation()  # we keep all the default values


class _InternalHeatGains:
    """Addition of internal heat gains to the space/thermal zone.
    Equipment heat gain and people heat gain are being considered.
    A distinction is made between an occupied time period and an unoccupied time
    period. The building is occupied from 9 o'clock till 17 o'clock. During the
    occupied period all equipment is turned on and 5 persons are present in the
    space/thermal zone. Outside the occupied period all equipment is turned off
    and nobody is present in the space/thermal zone.
    Calling method `add_internal_heat_gains` adds the internal heat gains to
    the space/thermal zone.
    """
    def __init__(self, weather_data: WeatherData, zone: FixedTemperatureZone) -> None:
        """Creates `_InternalHeatGains` instance.

        Parameters
        ----------
        weather_data:
            `WeatherData` object with the geographic location, the date of the
            design day, and the climatic design information. We need this here
            to create time schedules based on clock time.
        zone:
            The space/thermal zone to which the internal heat gains will be
            added.
        """
        self.zone = zone
        # The building is occupied from 9 o'clock till 17 o'clock. These times
        # must be converted to solar time for the cooling load calculations.
        self.t_start = convert_to_solar_seconds(
            clock_time=time(9),
            date=weather_data.date,
            L_loc=weather_data.location.L_loc,
            tz_loc='Europe/Brussels'
        )
        self.t_end = convert_to_solar_seconds(
            clock_time=time(17),
            date=weather_data.date,
            L_loc=weather_data.location.L_loc,
            tz_loc='Europe/Brussels'
        )

    def _create_office_equipment(self) -> OfficeEquipment:
        """Add office equipment to the space/thermal zone. All office equipment
        is turned on during the occupied period and all turned off outside the
        occupied period.
        """
        def _create_equipment_schedule() -> Callable[[float], float]:
            def schedule(t_sol_sec: float) -> float:
                if self.t_start <= t_sol_sec <= self.t_end:
                    return 1.0
                else:
                    return 0.0
            return schedule

        equipment = OfficeEquipment.create(
            ID='equipment',
            heat_gain_density=Q_(5, 'W / m**2'),
            floor_area=self.zone.floor_area,
            schedule=_create_equipment_schedule(),
            F_rad=Q_(10, 'pct')
        )
        return equipment

    def _create_equipment_heat_gain(self) -> EquipmentHeatGain:
        """Creates an `EquipmentHeatGain` object to which the office equipment
        is added.
        """
        equipment_ihg = EquipmentHeatGain(ID='equipment-heat-gain')
        office_equipment = self._create_office_equipment()
        equipment_ihg.add_equipment(office_equipment)
        return equipment_ihg

    def _create_people_heat_gain(self) -> PeopleHeatGain:
        """Creates a `PeopleHeatGain` object. During the occupied period 5
        persons are present. Outside the occupied period nobody is present in
        the space/thermal zone.
        """
        def _create_occupancy_schedule() -> Callable[[float], float]:
            def schedule(t_sol_sec: float) -> float:
                if self.t_start <= t_sol_sec <= self.t_end:
                    return 5
                else:
                    return 0
            return schedule

        people_ihg = PeopleHeatGain.create(
            ID='people-heat-gain',
            Q_dot_sen_person=Q_(70, 'W'),
            Q_dot_lat_person=Q_(45, 'W'),
            F_rad=Q_(60, 'pct'),
            schedule=_create_occupancy_schedule()
        )
        return people_ihg

    def add_internal_heat_gains(self) -> None:
        """Creates the internal heat gain objects and adds them to the
        space/thermal zone.
        """
        equipment_ihg = self._create_equipment_heat_gain()
        people_ihg = self._create_people_heat_gain()
        self.zone.add_internal_heat_gain(equipment_ihg, people_ihg)


class Zone2:

    def __init__(
        self,
        weather_data: WeatherData,
        T_int_des: Quantity,
        vez: VentilationZone
    ) -> None:
        self.weather_data = weather_data
        self.T_int_des = T_int_des
        self.vez = vez

        # STEP 1
        # Create the space envelope construction and the space (zone) itself:
        self.construction = _SpaceConstruction(weather_data, T_int_des)
        self.zone = self.construction.build()

        # STEP 2
        # Add local ventilation to the space.
        self.local_ventilation = _SpaceVentilation(self.zone, self.vez)
        self.local_ventilation.configure()

        # STEP 3
        # Add internal heat gains to the space.
        self.internal_heat_gains = _InternalHeatGains(weather_data, self.zone)
        self.internal_heat_gains.add_internal_heat_gains()

    def cooling_load(self, unit: str = 'W') -> pd.DataFrame:
        """Returns the hourly heat gains and cooling load in the single zone
        building on the design day.
        """
        df = self.zone.solve(unit=unit)
        return df

    @staticmethod
    def plot_cooling_load(df: pd.DataFrame) -> LineChart:
        """Returns a `LineChart` object with the total cooling load, the
        sensible cooling load, and the latent cooling load of the space/thermal
        zone. The loads are plotted along the x-axis in solar time.
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
