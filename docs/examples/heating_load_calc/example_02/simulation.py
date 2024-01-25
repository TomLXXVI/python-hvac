"""EXAMPLE 02 (part 2)
----------------------
SIMULATION OF AN "UNCONDITIONED" ZONE HEATED BY PANEL RADIATORS.

This example is similar to example_09 in docs/examples/cooling_load_calc, where
we simulate a space with an air cooling coil in cooling mode operation.

In this example we simulate a space during winter time. This space is now
equipped with panel radiators to heat the space.

The weather data for the simulation is taken from a TMY-datafile (csv-format)
that we have imported from https://re.jrc.ec.europa.eu/pvg_tools/en/#api_5.2 (see
also module hvac sun.tmy.py). The program expects to see the same column titles
as in the csv-file downloaded from this web application. That is:
- 'G(h)' is the column title for solar irradiance on the horizontal surface
- 'T2m' is the column title for dry-bulb outdoor air temperature
- 'RH' is the column title for outdoor air relative humidity.
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

Q_ = Quantity


def main():
    # Get climatic design data for creating and configuring the construction
    # assemblies based on peak-winter design conditions:
    climatic_design_data = ClimateDesignData(
        T_ext_d=Q_(-5.7, 'degC'),
        T_ext_an=Q_(11.5, 'degC'),
        T_ext_min=Q_(4.4, 'degC')
    )

    # Get weather data from a TMY-datafile for the simulation:
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

    # Set the indoor design temperature for creating and configuring the
    # construction assemblies based on peak-winter design conditions::
    T_int_des = Q_(20, 'degC')

    # Create our single-zone building model:
    building = SomeSingleZoneBuilding(climatic_design_data, weather_data, T_int_des)

    # Set the fixed operating conditions of the heating system in the building:
    # --> We can play with these parameters to see their effect on the zone
    # air temperature.
    building.heating_system.set_fixed_operating_conditions(
        Vw_dot_frac=Q_(75, 'pct'),
        Tw_sup=Q_(55, 'degC')
    )

    # Solve for the zone air temperature in the building:
    df = building.solve()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', 800
    ):
        print(df)


class ExteriorBuildingElementFactory:
    """Creates and returns the exterior building elements of the single-zone
    building.
    """
    def __init__(
        self,
        climatic_design_data: ClimateDesignData,
        weather_data: WeatherData,
        T_int_des: Quantity
    ) -> None:
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
    """Configures the heating system of our single-zone building.

    The single-zone building has a design heating load around 8 kW. We will
    place 4 panel radiators in the zone. From a catalog we select 4 identical
    panel radiators with a height of 600 mm and a length of 1000 mm, having a
    nominal heat output of 1.861 kW when the nominal supply water temperature is
    70 °C, the nominal return water temperature is 50 °C, and the nominal indoor
    air temperature is 20 °C. The radiator exponent is 1.35.
    """
    def __init__(self):
        self.radiator = PanelRadiator(
            Qe_dot_nom=Q_(1.861, 'kW'),
            Tw_sup_nom=Q_(70, 'degC'),
            Tw_ret_nom=Q_(50, 'degC'),
            Ti_nom=Q_(20, 'degC'),
            n_exp=1.35
        )
        self.num_rad: int = 4
        self.Vw_dot: Quantity | None = None
        self.Tw_sup: Quantity | None = None

    def set_fixed_operating_conditions(self, Vw_dot_frac: Quantity, Tw_sup: Quantity):
        """Set the fixed operating conditions of the heating system.

        Parameters
        ----------
        Vw_dot_frac:
            Fraction of the nominal water volume flow rate through the radiators.
        Tw_sup:
            The supply water temperature of the radiators.
        """
        self.Vw_dot = Vw_dot_frac.to('frac').m * self.radiator.Vw_dot_nom
        self.Tw_sup = Tw_sup

    def Qe_dot(self, _, T_zone: float) -> Quantity:
        """Returns the heat rate emitted by the panel radiators given the zone
        air temperature `T_zone` in Kelvins.

        Notes
        -----
        This function will be coupled to the nodal thermal zone model of our
        single-zone building (see class `SomeSingleZoneBuilding` below, in its
        method `solve`). This function must satisfy the call signature imposed
        by the nodal thermal zone model (see docstring of method `solve` of
        class `UnconditionedZone`).

        Also note that the heating system supplies heat to the zone air. The
        implemented sign convention demands that heat supplied by the heating
        system to the zone air must be associated with a negative sign (while
        heat extracted from the zone air by a cooling system must be associated
        with a positive sign).
        """
        T_zone = Q_(T_zone, 'K')
        Qe_dot = self.radiator(T_zone, Vw_dot=self.Vw_dot, Tw_sup=self.Tw_sup)
        Qe_dot *= self.num_rad
        return -Qe_dot


class SomeSingleZoneBuilding:
    """Models a very simple single-zone building with an interior floor area
    of 100 m² and a ceiling height of 3 m. The building has four exterior walls,
    a roof, and a floor. All the walls have the same interior dimensions:
    length is 10 m and height is equal to the ceiling height. The building has
    no windows.
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
            Contains the climatic data needed to create the exterior surface of
            exterior building elements, which are here instances of the class
            `ExteriorBuildingElement` in subpackage `cooling_load_calc`. Note
            that this class is completely different from the identically named
            class in subpackage `heating_load_calc`, which cannot be used for
            simulation, but only for the calculation of the conduction heat loss
            through the exterior building envelope according to the implemented
            heat load calculation method (see standard EN 12831-1).
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
            C_tsn=Q_(300, 'kJ / (K * m**2)'),
            A_tsn=Q_(100, 'm**2'),
            R_tsn=Q_(0.015, 'K * m**2 / W')
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

        # Add the ventilation/infiltration to the zone:
        self.zone.add_ventilation()

        # Add the heating system to the zone:
        self.heating_system = HeatingSystem()

    def solve(self) -> pd.DataFrame:
        self.zone.solve(
            Q_dot_sys_fun=self.heating_system.Qe_dot,
            num_cycles=15,
        )
        return self.zone.temperature_heat_gain_table()


if __name__ == '__main__':
    main()
