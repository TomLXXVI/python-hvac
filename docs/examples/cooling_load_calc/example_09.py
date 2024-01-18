"""
EXAMPLE 09
----------
SIMULATION OF AN "UNCONDITIONED" ZONE WITH AN AIR COOLING COIL.
The class `UnconditionedZone` represents a thermal model of a zone (a space or
a single-zone building) in which the zone air temperature is not fixed, but
depends on the heat gains in the zone. The actual zone air temperature at a
given time during the day being considered will follow from an energy balance
containing the momentary heat gains to the zone air and the momentary rate at
which heat is extracted from the zone air by the cooling system in the zone.
In this example the cooling system is represented by an air-cooling coil in
the zone. To model this air-cooling coil, being an air-to-water coil, we will
use in this example the `PlainFinTubeAirToWaterCounterFlowHeatExchanger` class.
The cooling capacity of the air-cooling coil depends on the entering water
temperature, the volume flow rate of water, the entering air temperature, which
is also the zone air temperature, and the volume flow rate of air through the
air cooling coil. Except the entering air temperature or zone air temperature,
we consider the other three operating conditions as being fixed.
The air-cooling coil is modeled inside the class `SomeAirCoil`. This class has
a method that takes the zone air temperature and returns the cooling capacity
of the air-cooling coil.
The zone is modeled inside the class `SomeZone`. Inside this class we will
instantiate an object of class `SomeAirCoil`, representing the air-cooling
coil being installed in the zone. For solving the energy balance of the zone for
the zone air temperature at each hour of the considered day, we will call the
method on the `SomeAirCoil` instance to get the rate at which the air-cooling
coil extracts heat from the zone air.

Notes
-----
To be able to solve the energy balance of the zone air at a given time moment `t`,
the zone air temperature at the previous time moment `t - dt` is used to calculate
the cooling capacity of the air-cooling coil at time moment `t`, and this value
is then used to calculate the zone air temperature at time moment `t`.
The implementation to solve the set of node equations in the nodal thermal model
of the zone uses matrix algebra, i.e., `[A] * [X] = [B]` with matrix `[A]` being
the coefficient matrix and `[X]` the column vector with the unknown node
temperatures at time moment `t`. As it is not possible to write the function of
the air-cooling coil's cooling capacity `Q_dot_cc = f(T_zone)` as `a * T_zone`
with `a` being a constant in the coefficient matrix [A], this function is put
in the input matrix [B], where it can only be solved with the zone air
temperature at the previous time moment `t - dt`.
"""
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.fluids import Fluid, FluidState, HumidAir
from hvac.heat_exchanger.fintube.continuous_fin import (
    PlainFinTubeAirToWaterCounterFlowHeatExchanger as AirCoil
)
from hvac.sun import Location, ClimateType, ReferenceDates
from hvac.cooling_load_calc import (
    WeatherData,
    ExteriorBuildingElement,
    wtcb,
    ConstructionAssembly,
    HeatFlowDirection,
    UnconditionedZone,
    VentilationZone
)
from hvac.charts import LineChart


Q_ = Quantity
Water = Fluid('Water')


class SomeAirCoil:
    """Encapsulates an `AirCoil` (`PlainFinTubeAirToWaterCounterFlowHeatExchanger`)
    instance that models a given air-cooling coil and contains a method that
    returns the cooling capacity of this air-cooling coil for a given zone air
    temperature at the air inlet of the air-cooling coil.
    """
    def __init__(self) -> None:
        """Creates the model of the air-cooling coil."""
        self.air_coil = AirCoil(
            width=Q_(900, 'mm'),
            height=Q_(180, 'mm'),
            num_rows=3,
            pitch_trv=Q_(25.4, 'mm'),
            pitch_lon=Q_(22.0, 'mm'),
            d_i=Q_(8.422, 'mm'),
            d_o=Q_(10.2, 'mm'),
            t_fin=Q_(0.3302, 'mm'),
            fin_density=1 / Q_(3.175, 'mm'),
            num_circuits=2
        )
        self.water_in: FluidState | None = None
        self.m_dot_w: Quantity | None = None
        self.V_dot_a: Quantity | None = None
        # Internal assumption: the relative humidity of the entering air is
        # always 50 %.
        self.RH_a_in: Quantity = Q_(50, 'pct')

    def set_fixed_operating_conditions(
        self,
        T_w_in: Quantity,
        V_dot_w: Quantity,
        V_dot_a: Quantity
    ) -> None:
        """Sets the fixed operating conditions of the air-cooling coil.

        Parameters
        ----------
        T_w_in:
            The entering water temperature.
        V_dot_w:
            The volume flow rate of water through the air-cooling coil.
        V_dot_a:
            The volume flow rate of air through the air-cooling coil.
        """
        self.water_in = Water(T=T_w_in, P=Q_(2, 'bar'))
        self.m_dot_w = V_dot_w * self.water_in.rho
        self.V_dot_a = V_dot_a

    def Q_dot_fun(self, T_a_in: Quantity) -> Quantity:
        """Returns the cooling capacity of the air-cooling coil for a given
        entering air temperature `T_a_in`, all other operating conditions
        being fixed.
        """
        air_in = HumidAir(Tdb=T_a_in, RH=self.RH_a_in)
        m_dot_a = self.V_dot_a * air_in.rho
        self.air_coil.set_operating_conditions(
            air_in=air_in,
            water_in=self.water_in,
            air_m_dot=m_dot_a,
            water_m_dot=self.m_dot_w
        )
        output = self.air_coil.rate()
        return output['Q_dot'].to('W')


class SomeZone:
    """Encapsulates an `UnconditionedZone` object, representing a simple
    single zone conditioned building consisting of 4 exterior walls and a roof
    at a given geographic location on a given day of the year.
    Also encapsulates an instance of `SomeAirCoil`, representing an air-cooling
    coil for cooling the zone.
    As such, this class models a single-zone building equipped with a given
    air-cooling coil for conditioning the zone air.
    Using the climatic design information for the location on the given day
    and the cooling capacity of the air-cooling coil as a function of the zone
    air temperature, the method `solve` of this class, will return a Pandas
    DataFrame table with the heat gains and the resulting zone air temperature
    at each hour of the selected day in this zone.
    """
    class ExteriorBuildingElements:
        """Inner class used to create the exterior building elements of the
        zone.
        """
        def __init__(
            self,
            weather_data: WeatherData,
            T_zone_des: Quantity,
            gross_area_sw: Quantity,
            gross_area_ww: Quantity,
            gross_area_nw: Quantity,
            gross_area_ew: Quantity,
            gross_area_rf: Quantity
        ) -> None:
            """Creates the exterior building elements of the zone:
            - south wall
            - west wall
            - north wall
            - east wall
            - roof

            Parameters
            ----------
            weather_data:
                `WeatherData` object that contains the climatic design
                information needed for instantiating the exterior building
                elements.
            T_zone_des:
                Design value of the zone air temperature needed for
                instantiating the exterior building elements (for determining
                the thermal resistance of the construction assemblies).
            gross_area_sw:
                The total surface area of the south wall.
            gross_area_ww:
                Same for the west wall.
            gross_area_nw:
                Same for the north wall.
            gross_area_ew:
                Same for the east wall.
            gross_area_rf:
                Same for the roof.
            """
            self.weather_data = weather_data
            self.T_zone_des = T_zone_des
            self.gross_area_sw = gross_area_sw
            self.gross_area_ww = gross_area_ww
            self.gross_area_nw = gross_area_nw
            self.gross_area_ew = gross_area_ew
            self.gross_area_rf = gross_area_rf
            # Create the construction assemblies of the exterior building
            # elements:
            T_int = self.T_zone_des.to('K')
            T_ext = self.weather_data.T_db_max.to('K')
            T_asp = (T_int + T_ext) / 2
            self.constr_assem_wall = wtcb.exterior_walls.create_ext_wall_wtcb_F1(
                t_ins=Q_(12, 'cm'),
                T_ext=T_ext,
                T_int=T_int,
                T_asp=T_asp,
                dT_asp=Q_(5, 'K'),
            )
            self.constr_assem_roof = wtcb.roofs.create_roof_wtcb_F1(
                t_ins=Q_(16, 'cm'),
                heat_flow_dir=HeatFlowDirection.DOWNWARDS,
                T_ext=T_ext,
                T_int=T_int,
                T_asp=T_asp,
                dT_asp=Q_(5, 'K')
            )
            self.window_props = wtcb.WindowPropertiesShelf.load(
                ID='window-5a-operable-wood/vinyl'
            )
            # Create the exterior building elements:
            self.south_wall = self._create_south_wall()
            self.west_wall = self._create_west_wall()
            self.north_wall = self._create_north_wall()
            self.east_wall = self._create_east_wall()
            self.roof = self._create_roof()

        def _create_south_wall(self):
            ext_wall = ExteriorBuildingElement.create(
                ID='south_wall',
                T_zone=lambda t_sol_sec: self.T_zone_des,
                # Note: in a `UnconditionedZone` object this parameter is
                # actually of no importance, as in a `UnconditionedZone` object
                # the zone air temperature is an unknown and it will be determined
                # from a energy balance of the zone air node in the thermal model
                # of the zone. However, the method `create` of class
                # `ExteriorBuildingElement` requires a value for this parameter,
                # so we need to give it some value.
                constr_assem=self.constr_assem_wall,
                gross_area=self.gross_area_sw,
                weather_data=self.weather_data,
                gamma=Q_(0, 'deg'),
                beta=Q_(90, 'deg')
            )
            ext_wall.add_window(
                ID='window_south',
                width=Q_(0.5 * self.gross_area_sw.to('m**2').m, 'm'),
                height=Q_(1, 'm'),
                # The window takes up 50 % of the wall's gross area.
                props=self.window_props
            )
            return ext_wall

        def _create_west_wall(self):
            ext_wall = ExteriorBuildingElement.create(
                ID='west_wall',
                T_zone=lambda t_sol_sec: self.T_zone_des,
                constr_assem=self.constr_assem_wall,
                gross_area=self.gross_area_ww,
                weather_data=self.weather_data,
                gamma=Q_(90, 'deg'),
                beta=Q_(90, 'deg')
            )
            ext_wall.add_door(
                ID='ext_door_west',
                width=Q_(1, 'm'),
                height=Q_(2.2, 'm'),
                constr_assem=ConstructionAssembly.create(
                    ID='ext_door_assem',
                    U=Q_(4.0, 'W / (K * m**2)')
                )
            )
            return ext_wall

        def _create_north_wall(self):
            ext_wall = ExteriorBuildingElement.create(
                ID='north_wall',
                T_zone=lambda t_sol_sec: self.T_zone_des,
                constr_assem=self.constr_assem_wall,
                gross_area=self.gross_area_nw,
                weather_data=self.weather_data,
                gamma=Q_(0, 'deg'),
                beta=Q_(90, 'deg')
            )
            return ext_wall

        def _create_east_wall(self):
            ext_wall = ExteriorBuildingElement.create(
                ID='east_wall',
                T_zone=lambda t_sol_sec: self.T_zone_des,
                constr_assem=self.constr_assem_wall,
                gross_area=self.gross_area_ew,
                weather_data=self.weather_data,
                gamma=Q_(-90, 'deg'),
                beta=Q_(90, 'deg')
            )
            window_props = wtcb.WindowPropertiesShelf.load(
                ID='window-5a-operable-wood/vinyl'
            )
            ext_wall.add_window(
                ID='window_east',
                width=Q_(0.5 * self.gross_area_ew.to('m**2').m, 'm'),
                height=Q_(1, 'm'),
                # The window takes up 50 % of the wall's gross area.
                props=window_props
            )
            return ext_wall

        def _create_roof(self):
            roof = ExteriorBuildingElement.create(
                ID='roof',
                T_zone=lambda t_sol_sec: self.T_zone_des,
                constr_assem=self.constr_assem_roof,
                gross_area=self.gross_area_rf,
                weather_data=self.weather_data,
                gamma=Q_(0, 'deg'),
                beta=Q_(0, 'deg')
            )
            return roof

    def __init__(self, weather_data: WeatherData) -> None:
        """Creates the `UnconditionedZone` object and the `SomeAirCoil` object
        to model the single-zone building, equipped with an air-cooling coil.
        """
        width = Q_(10, 'm')   # interior width of the space
        length = Q_(10, 'm')  # interior length of the space
        height = Q_(3, 'm')   # interior height of the space

        floor_area = width * length

        self.zone = UnconditionedZone.create(
            ID='zone',
            weather_data=weather_data,
            floor_area=floor_area,
            height=height,
            ventilation_zone=VentilationZone.create(ID='vez'),
            # characteristics of the interior thermal mass in the space:
            C_tsn=Q_(300, 'kJ / (K * m**2)'),
            A_tsn=floor_area,
            R_tsn=Q_(0.015, 'K * m**2 / W')
            # unit thermal resistance between the interior thermal mass and
            # the zone air (determined by the kind of floor covering and the
            # convective thermal resistance between the floor and zone air)
        )

        h_wall = height + Q_(0.5, 'm')  # exterior height of the space
        w_wall = width + Q_(0.5, 'm')   # exterior width of the space
        l_wall = length + Q_(0.5, 'm')  # exterior length of the space

        # Set the gross area of the exterior building elements:
        gross_area_wall_1 = h_wall * w_wall
        gross_area_wall_2 = h_wall * l_wall
        gross_area_roof = w_wall * l_wall

        # Create the exterior building elements:
        ext_build_elems = SomeZone.ExteriorBuildingElements(
            weather_data=weather_data,
            T_zone_des=Q_(24, 'degC'),
            gross_area_sw=gross_area_wall_1,
            gross_area_ww=gross_area_wall_2,
            gross_area_nw=gross_area_wall_1,
            gross_area_ew=gross_area_wall_2,
            gross_area_rf=gross_area_roof
        )

        # Add the exterior building elements to the zone:
        self.zone.add_ext_build_elem([
            ext_build_elems.south_wall,
            ext_build_elems.west_wall,
            ext_build_elems.north_wall,
            ext_build_elems.east_wall,
            ext_build_elems.roof
        ])

        # Add default space ventilation to the zone:
        self.zone.add_ventilation()

        # Create the air-cooling coil in the zone:
        self.air_coil: SomeAirCoil = SomeAirCoil()

    # noinspection PyUnusedLocal
    def _Q_dot_sys_fun(self, t_sol_sec: float, T_zone: float) -> Quantity:
        """This function is passed to the thermal zone model and will be called
        internally by it. The signature of this function must have the following
        two arguments: (1) solar time in seconds from midnight (0 s)
        (a float), and (2) the zone air temperature in degrees Kelvin (also a
        float). The function must return the cooling capacity of the air-cooling
        coil as a `Quantity` object.
        """
        # noinspection PyBroadException
        try:
            Q_dot = self.air_coil.Q_dot_fun(Q_(T_zone, 'K'))
        except Exception:
            return Q_(0.0, 'W')
        else:
            return Q_dot.to('W')

    def solve(
        self,
        T_w_in: Quantity,
        V_dot_w: Quantity,
        V_dot_a: Quantity
    ) -> pd.DataFrame:
        """Calculates the zone air temperature and the heat gains in the zone at
        each hour of the selected day and returns the results in a Pandas
        DataFrame table.

        Parameters
        ----------
        T_w_in:
            The entering water temperature.
        V_dot_w:
            The volume flow rate of water through the air-cooling coil.
        V_dot_a:
            The volume flow rate of air through the air-cooling coil.
        """
        self.air_coil.set_fixed_operating_conditions(T_w_in, V_dot_w, V_dot_a)
        self.zone.solve(
            F_rad=Q_(0.46, 'frac'),
            Q_dot_sys_fun=self._Q_dot_sys_fun,
            dt_hr=1.0,
            num_cycles=15
        )
        df = self.zone.temperature_heat_gain_table()
        return df


def main_01():
    """Solves the thermal model of the zone with the installed air-cooling
    coil and prints the table with the zone air temperatures and heat gains
    at each hour of the considered day of the year.
    """
    # Define the geographic location of the building:
    location = Location(
        fi=Q_(51.183, 'deg'),
        L_loc=Q_(3.8, 'deg'),
        altitude=Q_(8, 'm'),
        climate_type=ClimateType.MID_LATITUDE_SUMMER,
        timezone='Etc/GMT-1'
    )

    # Create the weather data using the climatic design information found in
    # the ASHRAE 2017 climate data tables for the given location:
    weather_data = WeatherData.create_from_climatic_design_data(
        location=location,
        date=ReferenceDates.get_date_for('Jul'),  # take the reference date for July
        T_db_des=Q_(26.7, 'degC'),
        T_db_rng=Q_(11.3, 'K'),
        T_wb_mc=Q_(19.2, 'degC'),
        T_wb_rng=Q_(4.7, 'K')
    )

    # Create the thermal model of our zone:
    zone = SomeZone(weather_data)

    # Solve the thermal model for the zone air temperature and heat gains at
    # each hour of the selected day:
    # Entering water temperature (EWT) 7 °C:
    df_7 = zone.solve(
        T_w_in=Q_(7, 'degC'),
        V_dot_w=Q_(853.02, 'L / hr'),
        V_dot_a=Q_(850, 'm**3 / hr')
    )
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', 800
    ):
        print(df_7)
    # Entering water temperature (EWT) 12 °C:
    df_12 = zone.solve(
        T_w_in=Q_(12, 'degC'),
        V_dot_w=Q_(853.02, 'L / hr'),
        V_dot_a=Q_(850, 'm**3 / hr')
    )
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', 800
    ):
        print(df_12)

    # Plot a line chart with the zone air temperature for EWT 7°C and EWT 12 °C:
    chart = LineChart()
    chart.add_xy_data(
        label='EWT 7°C',
        x1_values=df_7.index,
        y1_values=df_7['T_zone'],
        style_props={'marker': 'o'}
    )
    chart.add_xy_data(
        label='EWT 12°C',
        x1_values=df_12.index,
        y1_values=df_12['T_zone'],
        style_props={'marker': 'o'}
    )
    chart.add_legend()
    chart.x1.add_title('time index')
    chart.y1.add_title('T_zone, °C')
    chart.show()


def main_02():
    """Investigates the cooling capacity of the air-cooling coil as a function
    of the zone air temperature, while the other operating conditions are
    fixed.
    """
    air_coil = SomeAirCoil()
    air_coil.set_fixed_operating_conditions(
        T_w_in=Q_(7, 'degC'),
        V_dot_w=Q_(853.02, 'L / hr'),
        V_dot_a=Q_(850, 'm**3 / hr')
    )
    T_zone_rng = Q_(np.arange(15, 25, 0.25), 'degC')
    Q_dot_cc = [air_coil.Q_dot_fun(T_zone) for T_zone in T_zone_rng]
    Q_dot_cc = Quantity.from_list(Q_dot_cc)

    chart = LineChart()
    chart.add_xy_data(
        label='cooling capacity',
        x1_values=T_zone_rng.to('degC').m,
        y1_values=Q_dot_cc.to('W').m,
        style_props={'marker': 'o'}
    )
    chart.x1.add_title('entering air temperature, °C')
    chart.y1.add_title('cooling capacity, W')
    chart.show()


if __name__ == '__main__':
    main_01()
