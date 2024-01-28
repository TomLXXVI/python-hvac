"""EXAMPLE 02 (part 2)
----------------------
SIMULATION OF AN "UNCONDITIONED" ZONE HEATED BY PANEL RADIATORS CONTROLLED
BY AN ON/OFF CONTROLLER.

In this example a simulation is demonstrated of a space with a heating system
consisting of panel radiators. The heating load calculation of this space under
design conditions was done in part 1 of this example. The supply water volume
flow rate to the panel radiators is controlled by an on/off controller.

The weather data for this simulation is read from a TMY-datafile (csv-format)
that we have imported from https://re.jrc.ec.europa.eu/pvg_tools/en/#api_5.2 (see
also module hvac sun.tmy.py). Note that the code expects to see the same column
titles as in the csv-file imported from this web application. This means that:
- 'G(h)' must be the column title for the solar irradiance on the horizontal
   surface;
- 'T2m' must be the column title for the dry-bulb outdoor air temperature;
- 'RH' must be the column title for the outdoor air relative humidity.
"""
import warnings
from hvac.fluids import CoolPropWarning

warnings.filterwarnings('ignore', category=CoolPropWarning)

from datetime import date as Date
import pandas as pd
from hvac import Quantity
from hvac.sun import Location, ClimateType
from hvac.heating_load_calc import ClimateDesignData
from hvac.cooling_load_calc import (
    WeatherData,
    wtcb,
    HeatFlowDirection,
    ExteriorBuildingElement,
    InteriorBuildingElement,
    ConstructionAssembly,
    VentilationZone,
    UnconditionedZone
)
from hvac.radiant_emitter import PanelRadiator
from hvac.control.controller import OnOffController
from hvac.charts import LineChart, BarChart


Q_ = Quantity


def main():
    # Set the climatic design data and set the indoor design temperature needed
    # to create and to configure the construction assemblies based on the peak-
    # winter design conditions:
    climatic_design_data = ClimateDesignData(
        # design value of outdoor air temperature (in the coldest month):
        T_ext_d=Q_(-5.7, 'degC'),
        # annual mean outdoor air temperature:
        T_ext_an=Q_(11.5, 'degC'),
        # average of the daily minimum outdoor air temperatures during the
        # coldest month:
        T_ext_min=Q_(4.4, 'degC')
    )
    # design value for the indoor air temperature:
    T_int_des = Q_(20, 'degC')

    # Get the weather data from a TMY-datafile for doing the simulation:
    weather_data = WeatherData.create_from_tmy_data(
        location=Location(
            fi=Q_(50.911, 'deg'),
            L_loc=Q_(3.192, 'deg'),
            altitude=Q_(8, 'm'),
            climate_type=ClimateType.MID_LATITUDE_SUMMER,
            timezone='Etc/GMT-1'
        ),
        date=Date(2001, 1, 15),  # <-- select the day
        tmy_file='tmy_50.911_3.192_2005_2020.csv'
    )

    # Create the model of the building with the heating system:
    building = SomeSingleZoneBuilding(climatic_design_data, weather_data, T_int_des)

    # Solve this model for the zone air temperature and the heat losses/gains in
    # the building during the selected day:
    df = building.solve()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', 800
    ):
        print(df)

    # Print the total amount of thermal energy delivered by the heating system
    # to the zone air:
    _, Q_sys = building.zone.get_system_heat_transfer()
    print(
        f"Thermal energy delivered to the zone by the heating system "
        f"= {Q_sys.to('kWh'):~P.3f}"
    )

    # Plot the zone air temperature in a line chart:
    chart1 = LineChart(size=(8, 6))
    chart1.add_xy_data(
        label='T_zone',
        x1_values=df.index,
        y1_values=df['T_zone']
    )
    chart1.x1.add_title('time index')
    chart1.y1.add_title('temperature, °C')
    chart1.show()

    # Plot the heat rate supplied by the system:
    chart2 = BarChart(size=(8, 6))
    chart2.add_xy_data(
        label='Q_dot_sys',
        x1_values=df.index,
        y1_values=df['Q_dot_sys'],
        style_props={
            'width': 0.8,
            'align': 'edge'
        }
    )
    chart2.x1.add_title('time index')
    chart2.y1.add_title('heat rate, W')
    chart2.show()


class ExteriorBuildingElementFactory:
    """This is a separate class for creating the exterior building elements of
    the single-zone building.
    """
    def __init__(
        self,
        climatic_design_data: ClimateDesignData,
        weather_data: WeatherData,
        T_int_des: Quantity
    ) -> None:
        """Instantiate the class with the data that's needed to create the
        construction assemblies and the exterior building elements of the
        building.

        Parameters
        ----------
        climatic_design_data:
            Contains the data for creating the construction assemblies: the
            outdoor air design temperature, the annual average outdoor air
            temperature, and the average of the minimum outdoor air temperatures
            during the coldest month.
        weather_data:
             Contains the weather data (solar radiation data and outdoor air
             temperature) for the selected day for doing the simulation.
        T_int_des:
            The design-value of the indoor air temperature used in the heating
            load calculation of the building. We need this to create the
            construction assemblies.
        """
        # Input data needed to create the construction assemblies and the
        # exterior building elements:
        self.climatic_design_data = climatic_design_data
        self.weather_data = weather_data
        self.T_int_des = T_int_des

        # Create the construction assemblies of the exterior building elements:
        self._ca_ew = self._create_constr_assem_ext_wall()
        self._ca_rf = self._create_constr_assem_roof()
        self._ca_fl = self._create_constr_assem_floor()

    def _create_constr_assem_ext_wall(self) -> ConstructionAssembly:
        ca = wtcb.exterior_walls.create_ext_wall_wtcb_F1(
            t_ins=Q_(12, 'cm'),
            T_ext=self.climatic_design_data.T_ext_d,
            T_int=self.T_int_des
        )
        return ca

    def _create_constr_assem_roof(self) -> ConstructionAssembly:
        ca = wtcb.roofs.create_roof_wtcb_F1(
            t_ins=Q_(18, 'cm'),
            heat_flow_dir=HeatFlowDirection.UPWARDS,
            T_ext=self.climatic_design_data.T_ext_d,
            T_int=self.T_int_des
        )
        return ca

    def _create_constr_assem_floor(self) -> ConstructionAssembly:
        ca = wtcb.floors.create_floor_wtcb_F3(
            t_ins=Q_(8, 'cm'),
            heat_flow_dir=HeatFlowDirection.DOWNWARDS,
            T_ext=self.climatic_design_data.T_ext_min,
            T_int=self.T_int_des
        )
        return ca

    def _create_exterior_wall(
        self,
        name: str,
        L_ext: Quantity,
        H_ext: Quantity,
        gamma: Quantity,
        ca: ConstructionAssembly
    ) -> ExteriorBuildingElement:
        ew = ExteriorBuildingElement.create(
            ID=f'exterior-wall-{name}',
            T_zone=lambda t_sol_sec: self.T_int_des,
            # Note: in a `UnconditionedZone` object this parameter is
            # actually of no importance, as in a `UnconditionedZone` object
            # the zone air temperature is an unknown and it will be determined
            # from a energy balance of the zone air node in the thermal model
            # of the zone. However, the method `create` of class
            # `ExteriorBuildingElement` requires a value for this parameter,
            # so we need to give it some value.
            constr_assem=ca,
            gross_area=L_ext * H_ext,
            weather_data=self.weather_data,
            gamma=gamma,
            beta=Q_(90, 'deg')
        )
        return ew

    def _create_roof(
        self,
        W_ext: Quantity,
        L_ext: Quantity,
        ca: ConstructionAssembly
    ) -> ExteriorBuildingElement:
        rf = ExteriorBuildingElement.create(
            ID='roof',
            T_zone=lambda t_sol_sec: self.T_int_des,
            constr_assem=ca,
            gross_area=L_ext * W_ext,
            weather_data=self.weather_data,
            gamma=Q_(0, 'deg'),
            beta=Q_(0, 'deg')
        )
        return rf

    def _create_floor(
        self,
        W_ext: Quantity,
        L_ext: Quantity,
        ca: ConstructionAssembly
    ) -> InteriorBuildingElement:
        # We use an `InteriorBuildingElement` object to represent the floor
        # and assign the ground temperature to its `T_adj` attribute.
        # For the ground temperature, the average minimum outdoor temperature
        # during the coldest month is used.
        fl = InteriorBuildingElement.create(
            ID='floor',
            T_zone=lambda t_sol_sec: self.T_int_des,
            T_adj=lambda t_sol_sec: self.climatic_design_data.T_ext_min,
            constr_assem=ca,
            gross_area=L_ext * W_ext
        )
        return fl

    def get_south_wall(self, L_int: Quantity, H_int: Quantity) -> ExteriorBuildingElement:
        L_ext = L_int + 2 * self._ca_ew.thickness
        H_ext = H_int + self._ca_rf.thickness
        sew = self._create_exterior_wall(
            'south', L_ext, H_ext,
            Q_(0, 'deg'), self._ca_ew
        )
        return sew

    def get_west_wall(self, L_int: Quantity, H_int: Quantity) -> ExteriorBuildingElement:
        L_ext = L_int + 2 * self._ca_ew.thickness
        H_ext = H_int + self._ca_rf.thickness
        wew = self._create_exterior_wall(
            'west', L_ext, H_ext,
            Q_(90, 'deg'), self._ca_ew
        )
        return wew

    def get_north_wall(self, L_int: Quantity, H_int: Quantity) -> ExteriorBuildingElement:
        L_ext = L_int + 2 * self._ca_ew.thickness
        H_ext = H_int + self._ca_rf.thickness
        new = self._create_exterior_wall(
            'north', L_ext, H_ext,
            Q_(180, 'deg'), self._ca_ew
        )
        return new

    def get_east_wall(self, L_int: Quantity, H_int: Quantity) -> ExteriorBuildingElement:
        L_ext = L_int + 2 * self._ca_ew.thickness
        H_ext = H_int + self._ca_rf.thickness
        eew = self._create_exterior_wall(
            'east', L_ext, H_ext,
            Q_(270, 'deg'), self._ca_ew
        )
        return eew

    def get_roof(self, W_int: Quantity, L_int: Quantity) -> ExteriorBuildingElement:
        W_ext = W_int + 2 * self._ca_ew.thickness
        L_ext = L_int + 2 * self._ca_ew.thickness
        rf = self._create_roof(W_ext, L_ext, self._ca_rf)
        # Add skylight:
        rf.add_window(
            ID='skylight',
            width=Q_(5, 'm'),
            height=Q_(5, 'm'),
            props=wtcb.WindowPropertiesShelf.load('window-5a-operable-wood/vinyl')
        )
        return rf

    def get_floor(self, W_int: Quantity, L_int: Quantity) -> InteriorBuildingElement:
        W_ext = W_int + 2 * self._ca_ew.thickness
        L_ext = L_int + 2 * self._ca_ew.thickness
        fl = self._create_floor(W_ext, L_ext, self._ca_fl)
        return fl


class HeatingSystem:
    """Represents the heating system in our single-zone building.

    The single-zone building has a design heating load around 8 kW. From a
    catalog we select 4 identical panel radiators with a height of 600 mm and a
    length of 1000 mm, having a nominal heat output of 1.861 kW when the nominal
    supply water temperature is 70 °C, the nominal return water temperature is
    50 °C, and the nominal indoor air temperature is 20 °C. The radiator
    exponent is 1.35.

    The setpoint of the on/off controller, which will control the zone air
    temperature, is 22 °C. The controller has a symmetrical dead band of 1 K
    around the setpoint (-0.5 K / +0.5 K).
    """
    def __init__(self):
        # Configure the panel radiator:
        self.radiator = PanelRadiator(
            Qe_dot_nom=Q_(1.861, 'kW'),
            Tw_sup_nom=Q_(70, 'degC'),
            Tw_ret_nom=Q_(50, 'degC'),
            Ti_nom=Q_(20, 'degC'),
            n_exp=1.35
        )
        # Configure the on/off controller:
        self.controller = OnOffController(
            SP=Q_(22, 'degC'),
            HL_offset=Q_(0.5, 'K'),
            LL_offset=Q_(0.5, 'K'),
            PV_range=(Q_(0, 'degC'), Q_(50, 'degC')),
            dt=Q_(1/6, 'hr'),  # --> time step of the control action
            ctrl_dir='inverse'
        )
        # Set the number of panel radiators in the zone:
        self.num_rad: int = 4
        # Set the maximum water volume flow rate through one radiator:
        self.Vw_dot_max = 1.00 * self.radiator.Vw_dot_nom
        # Set the supply water temperature of the radiators:
        self.Tw_sup = self.radiator.Tw_sup_nom

    def Qe_dot(self, t_sol_sec: float, T_zone: float) -> Quantity:
        """Returns the heat rate emitted by the panel radiators given the zone
        air temperature `T_zone` in Kelvins.

        Notes
        -----
        This function will be coupled to the nodal thermal zone model of our
        single-zone building (see class `SomeSingleZoneBuilding` below, in its
        method `solve`). This function must therefore satisfy the call signature
        imposed by the nodal thermal zone model (see docstring of method `solve`
        of class `UnconditionedZone`).

        Also note that the heating system supplies heat directly to the zone air.
        The implemented sign convention demands that heat supplied by the heating
        system to the zone air must be associated with a negative sign (while
        heat extracted from the zone air by a cooling system must be associated
        with a positive sign).
        """
        T_zone = Q_(T_zone, 'K')
        # Get the output signal (being a percentage) from the controller:
        out = self.controller(t_sol_sec, T_zone)
        # The water volume flow rate through one radiator:
        Vw_dot = out * self.Vw_dot_max
        # The heat rate output of one radiator...
        Qe_dot = self.radiator(
            T_zone,
            Vw_dot=Vw_dot,
            Tw_sup=self.Tw_sup
        )
        # ...is multiplied with the number of radiators in the zone:
        Qe_dot *= self.num_rad
        # print(
        #     f"{out.to('pct'):~P.0f}, "
        #     f"{Vw_dot.to('L / min'):~P.0f}, "
        #     f"{Qe_dot.to('kW'):~P.3f}"
        # )
        return -Qe_dot


class SomeSingleZoneBuilding:
    """Models a very simple single-zone building with an interior floor area
    of 100 m² and a ceiling height of 3 m. The building has four exterior walls,
    a roof, and a floor. All the walls have the same interior dimensions:
    their length is 10 m and their height is equal to the ceiling height. On the
    roof, there is a flat skylight with a window area of 25 m².
    """
    def __init__(
        self,
        climatic_design_data: ClimateDesignData,
        weather_data: WeatherData,
        T_int_des: Quantity
    ) -> None:
        """Instantiates the building model.

        Parameters
        ----------
        climatic_design_data:
            Instance of class `ClimateDesignData` in subpackage `heating_load_calc`.
            Contains the climatic winter design data needed to create the
            construction assemblies.
        weather_data:
            Instance of class `WeatherData` in subpackage `cooling_load_calc`.
            Contains the climatic data needed to determine the sol-air
            temperature at the exterior surface of the exterior building
            elements, which are instances of the class `ExteriorBuildingElement`
            in the subpackage `cooling_load_calc`. Note that this class is
            completely different from the identically named class in subpackage
            `heating_load_calc`, which cannot be used for simulation, but only
            for the calculation of the conduction heat loss through the exterior
            building envelope under design conditions according to the
            implemented heat load calculation method (standard EN 12831-1).
        T_int_des:
            The indoor air temperature used in the heating load design
            calculation of the building.
        """
        # Create the `UnconditionedZone` object:
        self.zone = UnconditionedZone.create(
            ID='my-test-zone',
            weather_data=weather_data,
            floor_area=Q_(100, 'm**2'),
            height=Q_(3, 'm'),
            ventilation_zone=VentilationZone.create(ID='my-test-ventilation-zone'),
            C_tsn=Q_(500, 'kJ / (K * m**2)'),
            A_tsn=Q_(100, 'm**2'),
            R_tsn=Q_(0.15, 'K * m**2 / W')
        )

        # Create the exterior building elements of the zone:
        ebf = ExteriorBuildingElementFactory(
            climatic_design_data,
            weather_data,
            T_int_des
        )
        L_int, H_int = Q_(10, 'm'), Q_(3, 'm')
        ext_build_elems = [
            ebf.get_south_wall(L_int, H_int),
            ebf.get_west_wall(L_int, H_int),
            ebf.get_north_wall(L_int, H_int),
            ebf.get_east_wall(L_int, H_int),
            ebf.get_roof(L_int, L_int)
        ]

        # Add the exterior building elements to the zone:
        self.zone.add_ext_build_elem(ext_build_elems)

        # Add the floor to the zone:
        # --> Note that the floor is represented as an `InteriorBuildingElement`
        # object.
        self.zone.add_int_build_elem(ebf.get_floor(L_int, L_int))

        # Add ventilation/infiltration to the zone:
        self.zone.add_ventilation()

        # Add the heating system to the zone:
        self.heating_system = HeatingSystem()

    def solve(self) -> pd.DataFrame:
        """Solve for the zone air temperature and the heat losses/gains in the
        zone. Returns a Pandas DataFrame with these results for each "calculation
        time moment" during the selected day. The calculations always start at
        midnight of the selected day (0 s). The number of "calculation time
        moments" will depend on the time step we set.
        """
        self.zone.solve(
            Q_dot_sys_fun=self.heating_system.Qe_dot,
            # --> Here we couple the heating system with the thermal zone model.
            num_cycles=15,
            dt_hr=1/6  # calculation time step = 10 min
        )
        return self.zone.temperature_heat_gain_table()


if __name__ == '__main__':
    main()
