from __future__ import annotations

import typing
from dataclasses import dataclass
import numpy as np
import pandas as pd
import scipy.integrate as integrate
import control as ct
from hvac import Quantity
from hvac.fluids.fluid_experimental import PureFluid
from .weather_data import WeatherData
from .building_elements import ExteriorBuildingElement, InteriorBuildingElement
from .fenestration import Window
from ..internal_heat_gains import InternalHeatGain
from ..ventilation import ThermalZoneVentilation, VentilationZone
from ..thermal_models import ThermalZoneModel, Units, TemperatureNode


Q_ = Quantity
Air = PureFluid('Air')

class ThermalZone:

    def __init__(self) -> None:
        self.name: str = ''
        self.floor_area: Quantity | None = None
        self.ceiling_height: Quantity | None = None
        self.unit_R_im: Quantity | None = None
        self.unit_C_im: Quantity | None = None
        self.factor_surf_im: float | None = None
        self.weather_data: WeatherData | None = None
        self.air_temperature: Quantity | None = None
        self.air_pressure: Quantity | None = None
        self.ext_build_elems: dict[str, ExteriorBuildingElement] = {}
        self.windows: dict[str, dict[str, Window]] = {}
        self.ext_doors: dict[str, dict[str, ExteriorBuildingElement]] = {}
        self.int_build_elems: dict[str, InteriorBuildingElement] = {}
        self.int_doors: dict[str, dict[str, InteriorBuildingElement]] = {}
        self.int_heat_gains: dict[str, InternalHeatGain] = {}
        self.ventilation: ThermalZoneVentilation | None = None
        self.ventilation_zone: VentilationZone | None = None
        self.RH_zone: Quantity | None = None
        self._model_factory: typing.Type[ThermalZoneModel] = ThermalZoneModel
        
        self.time: Quantity | None = None
        self.Q_gain_cond: dict[str, Quantity] | None = None  # conduction heat gains
        self.Q_gain_sol: dict[str, Quantity] | None = None  # solar heat gains
        self.Q_gain_ihg: dict[str, Quantity] | None = None  # internal heat gains
        self.Q_gain_vent: dict[str, Quantity] | None = None  # ventilation heat gains
        self.Q_gain_im: Quantity | None = None  # heat transfer between interior mass and zone air node
        
    @classmethod
    def create(
        cls,
        name: str,
        floor_area: Quantity,
        ceiling_height: Quantity,
        unit_R_im: Quantity,
        unit_C_im: Quantity,
        weather_data: WeatherData,
        factor_surf_im: float = 1.0,
        **kwargs
    ) -> ThermalZone:
        """Creates a `Zone` object.

        Parameters
        ----------
        name:
            Name of the temperature zone.
        floor_area:
            Floor area of the zone.
        ceiling_height:
            Ceiling height in the zone.
        unit_R_im:
            Unit thermal resistance of the interior mass surface. 
        unit_C_im:
            Unit thermal capacity of the interior mass, i.e. per unit area.
        weather_data:
            Instance of class `WeatherData` encapsulating the climatic design
            data for the selected design day, which was given when the 
            `WeatherData` object was instantiated.
        factor_surf_im:
            Multiplication factor for the interior mass surface area.
        
        Keyword Arguments
        -----------------
        air_temperature:
            Reference zone air temperature used for calculating the thermal
            capacity of the zone air. Default value is 20 °C.
        air_pressure:
            Reference zone air pressure used for calculating the thermal 
            capacity of the zone air. Default value is 101,325 Pa.
        RH_zone:
            Setpoint or desired value of the relative air humidity in the zone.
            Default value is 50 %.
        """
        zone = cls()
        zone.name = name
        zone.floor_area = floor_area
        zone.ceiling_height = ceiling_height
        zone.unit_R_im = unit_R_im
        zone.unit_C_im = unit_C_im
        zone.factor_surf_im = factor_surf_im
        zone.weather_data = weather_data
        zone.air_temperature = kwargs.get('air_temperature', Q_(20, 'degC'))
        zone.air_pressure = kwargs.get('air_pressure', Q_(101_325, 'Pa'))
        zone.RH_zone = kwargs.get('RH_zone', Q_(50, 'pct'))
        return zone

    @property
    def envelope_area(self) -> Quantity:
        """Returns the total area of the exterior building elements that 
        surround the thermal zone.
        """
        A_env = sum(ebe.gross_area for ebe in self.ext_build_elems.values())
        return A_env

    @property
    def volume(self) -> Quantity:
        """Returns the volume of the thermal zone."""
        V = self.floor_area * self.ceiling_height
        return V
    
    def get_input_name(self, key: str) -> str:
        """Returns the input label for the zone system. Only for internal use."""
        return f"{ThermalZoneModel.input_names[key]}@{self.name}"
    
    def add_ext_build_elems(self, *ext_build_elems: ExteriorBuildingElement) -> None:
        """Adds exterior building elements to the zone."""
        for ext_build_elem in ext_build_elems:
            self.ext_build_elems[ext_build_elem.name] = ext_build_elem
            # If the exterior building element has windows, add them to the
            # `windows` dictionary of the zone.
            if ext_build_elem.windows:
                self.windows[ext_build_elem.name] = ext_build_elem.windows
            # If the exterior building element has doors, add them to the
            # `ext_doors` dictionary of the zone.
            if ext_build_elem.doors:
                self.ext_doors[ext_build_elem.name] = ext_build_elem.doors
    
    def add_int_build_elems(self, *int_build_elems: InteriorBuildingElement) -> None:
        """Adds interior building elements to the zone."""
        for int_build_elem in int_build_elems:
            self.int_build_elems[int_build_elem.name] = int_build_elem
            # If the interior building element has doors, add them to the
            # `int_doors` dictionary of the zone.
            if int_build_elem.doors:
                self.int_doors[int_build_elem.name] = int_build_elem.doors
    
    def add_internal_heat_gains(self, *ihg: InternalHeatGain) -> None:
        """Adds an internal heat gain or a sequence of internal heat gains to
        the zone (see cooling_load_calc.core.internal_heat_gains).
        """
        self.int_heat_gains.update({ihg_.name: ihg_ for ihg_ in ihg})
    
    def add_ventilation(
        self,
        ventilation_zone: VentilationZone,
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_dot_open: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_ATD_d: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_sup: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_exh: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_comb: Quantity = Q_(0.0, 'm ** 3 / hr'),
        T_sup: typing.Callable[[float], Quantity] | None = None
    ) -> None:
        """Adds ventilation to the thermal zone.
        
        Parameters
        ----------
        ventilation_zone:
            The ventilation zone to which the thermal zone belongs. 
        n_min:
            Minimum air change rate required for the space for reasons of air
            quality/hygiene and comfort (EN 12831-1, B.2.10 - Table B.7).
            The default value applies to permanent dwelling areas (living rooms,
            offices) and a ceiling height less than 3 m.
        V_dot_open:
            External air volume flow rate into the space through large openings
            (EN 12831-1, Annex G).
        V_dot_ATD_d:
            Design air volume flow rate of the ATDs in the room (EN 12831-1,
            B.2.12). Only relevant if ATDs (Air Terminal Devices) are used for
            ventilation (i.e., passive devices that allow air flow through a
            building element; it does not include the air out- or inlets of 
            fan-assisted ventilation systems).
        V_dot_sup:
            Supply air volume flow rate from the ventilation system into the
            space.
        V_dot_exh:
            Exhaust ventilation air volume flow rate from the space.
        V_dot_comb:
            Air volume flow rate exhausted from the space that has not been
            included in the exhaust air volume flow of the ventilation system
            (typically, but not necessarily, combustion air if an open flue
            heater is present in the heated space).
        T_sup: optional
            Temperature of ventilation air supplied to the zone. If `None`, but
            `V_dot_sup` is not zero, the outdoor dry-bulb temperature will be
            assigned to `T_sup`. 
        """
        self.ventilation_zone = ventilation_zone
        self.ventilation_zone.add_thermal_zone(self)
        V_dot_trf_dict = {
            ibe.adj_zone_name: ibe.V_dot_trf 
            for ibe in self.int_build_elems.values()
            if ibe.V_dot_trf is not None
        }
        if not V_dot_trf_dict: V_dot_trf_dict = None
        if T_sup is None and V_dot_sup > 0.0:
            T_sup = self.weather_data.T_db
        self.ventilation = ThermalZoneVentilation.create(
            self, n_min, V_dot_open, V_dot_ATD_d, 
            V_dot_sup, V_dot_trf_dict, V_dot_exh, 
            V_dot_comb, T_sup
        )
    
    def create_model(self) -> ThermalZoneModel:
        """Creates the linear thermal network model of the zone."""
        # First create the "bare" zone model (only an interior mass node and 
        # a zone-air node).
        bare_model = self._model_factory(
            self.name,
            self.floor_area,
            self.ceiling_height,
            self.unit_R_im,
            self.unit_C_im,
            self.factor_surf_im,
            air_temperature=self.air_temperature,
            air_pressure=self.air_pressure,
        )
        # Collect all the building elements that surround the zone.
        ext_build_elem_models_ = None
        windows_ = None
        ext_doors_ = None
        int_build_elems_ = None
        int_doors_ = None
        if self.ext_build_elems:
            # Collect the exterior building elements surrounding the zone in a 
            # list and create their linear thermal network model.
            ext_build_elems = list(self.ext_build_elems.values())
            ext_build_elem_models_ = tuple(
                ext_build_elem.create_model(closed=True) 
                for ext_build_elem in ext_build_elems
            )
            # Collect the windows surrounding the zone in a single list.
            if self.windows:
                windows_ = []
                for key in self.windows.keys():
                    w_list = list(self.windows[key].values())
                    windows_.extend(w_list)
                windows_ = tuple(windows_)
            # Collect the exterior doors surrounding the zone in a single list.
            if self.ext_doors:
                ext_doors_ = []
                for key in self.ext_doors.keys():
                    ed_list = list(self.ext_doors[key].values())
                    ext_doors_.extend(ed_list)
                ext_doors_ = tuple(ext_doors_)
        if self.int_build_elems:
            # Collect the interior building elements surrounding the zone in
            # a single list.
            int_build_elems_ = tuple(self.int_build_elems.values())
            # Collect the interior doors surrounding the zone in a single list.
            if self.int_doors:
                int_doors_ = []
                for key in self.int_doors.keys():
                    id_list = list(self.int_doors[key].values())
                    int_doors_.extend(id_list)
                int_doors_ = tuple(int_doors_)
        # Assemble the zone model. This creates the final thermal network of the
        # zone.
        model = bare_model.assemble(
            ext_build_elem_models=ext_build_elem_models_,
            int_build_elems=int_build_elems_,
            windows=windows_, 
            exterior_doors=ext_doors_,
            interior_doors=int_doors_,
            ventilation=self.ventilation
        )
        return model

    def collect_input_data(
        self,
        num_cycles: int,
        setpoint: Quantity | None = None
    ) -> tuple[np.ndarray, dict[str, Quantity]]:
        """Returns all the input data that is needed for running a simulation of
        the zone system.

        Parameters
        ----------
        num_cycles:
            The number of days the selected (design) day needs to be repeated 
            before the simulation ends. When the simulation starts, the initial 
            values of all system state variables are set to zero, and the daily,
            24h cycle of system input values needs to be repeated a number of 
            days until the system shows a harmonic periodic behavior in time or 
            the zone-air temperature 'T_za' has settled around the setpoint. 
        setpoint:
            Setpoint of the zone-air temperature (a fixed value). If set to 
            `None`, the setpoint is ignored and the HVAC-system is considered to
            be turned off.

        Returns
        -------
        time:
            Numpy 1D-array holding the timeline for the simulation. The timeline
            counts the hours of `num_cycles` successive days (the selected
            design day repeats itself `num_cycles` times). 
        input_data: dict[str, Quantity]
            The keys of the dictionary are the input names of the zone system.
            The corresponding values of the dictionary are holding the input 
            values to the zone system at each time moment (hour) in the 
            timeline.
            The inputs to the zone system can be:
            1.  The setpoint temperature for the zone air if the setpoint has
                been set.
            2.  The sol-air temperatures at each exterior building element 
                surrounding the zone (these are multiple inputs: one input for
                each exterior building element, including exterior doors).
            3.  The dry-bulb outdoor air temperature.
            4.  The radiative heat gains in the zone.
            5.  The convective heat gains in the zone.
            6.  The zone-air temperatures of any adjacent zones (this is the
                case when the zone has interior building elements).
        """
        # Set up the time axis. The time axis is an array with `num_cycles` of
        # daily cycles that count each hour of the day (it is as if the weather 
        # and sun trajectory on the selected design day is repeated for 
        # `num_cycles` days in a row). 
        time_list = []
        for i in range(num_cycles):
            time = np.array([
                Q_(h, 'hr').to(Units.unit_t).m
                for h in range(i * 24, (i + 1) * 24)
            ])
            time_list.append(time)
        time = np.concatenate(time_list)
        # Collect temperature independent heat gains in the zone.
        self._collect_temp_indep_heat_gains()
        # Gather all input data needed for a simulation in a dict.
        input_data = {}
        input_data.update(self._collect_radiative_heat_gains(num_cycles))
        input_data.update(self._collect_convective_heat_gains(num_cycles))
        input_data.update(self._collect_sol_air_temperatures(num_cycles))
        if setpoint is not None:
            input_data.update(self._collect_setpoint_values(setpoint, num_cycles))
        if self.windows:
            input_data.update(self._collect_outdoor_air_temperatures(num_cycles))
        if self.int_build_elems:
            input_data.update(self._collect_adjacent_zone_temperatures(num_cycles))
        if self.ventilation and self.ventilation.V_dot_sup > 0.0:
            input_data.update(self._collect_supply_air_temperatures(num_cycles))
        return time, input_data

    def _collect_temp_indep_heat_gains(self) -> None:
        """Collects the solar and internal heat gains in the zone at each hour
        of the design day and assigns them to object data members.

        Attributes
        ----------
        self.Q_sol_conv: `Quantity`-array
            The convective solar heat gains at each hour of the selected design 
            day.
        self.Q_sol_rad: `Quantity`-array
            The radiative solar heat gains at each hour of the selected design 
            day.
        self.Q_ihg_conv: `Quantity`-array
            The convective sensible internal heat gains at each hour of the
            selected design day.
        self.Q_ihg_rad: `Quantity`-array
            The radiative sensible internal heat gains at each hour of the
            selected design day.
        self.Q_ihg_lat: `Quantity`-array
            The latent internal heat gains at each hour of the selected design 
            day.
        self.Q_vent_lat: `Quantity`-array
            The latent heat gains due to outdoor air infiltration / ventilation 
            at each hour of the selected design day.
        """
        # Collect solar heat gains through the windows of the zone.
        Q_sol_conv, Q_sol_rad = self._collect_solar_heat_gains()
        self.Q_gain_sol = {'conv': Q_sol_conv, 'rad': Q_sol_rad}
        # Collect internal heat gains in the zone.
        Q_ihg_conv, Q_ihg_rad, Q_ihg_lat = self._collect_internal_heat_gains()
        self.Q_gain_ihg = {'conv': Q_ihg_conv, 'rad': Q_ihg_rad, 'lat': Q_ihg_lat}
        # Collect latent ventilation heat gain in the zone.
        Q_vent_lat = self._collect_ventilation_latent_heat_gains()
        self.Q_gain_vent = {'lat': Q_vent_lat}
        
    def _collect_solar_heat_gains(self) -> tuple[Quantity, Quantity]:
        """Get the solar heat gains through all windows of the zone at each hour
        of the selected design day. The solar heat gains have a convective 
        fraction, which is heat gain directly transferred to the zone air, and a
        radiative fraction, which is heat gain transferred to the interior mass 
        of the zone.
        
        Returns
        -------
        Q_sol_conv_sum: Quantity-array
            The total convective solar heat gain through all windows of the zone
            at each hour of the selected design day.
        Q_sol_rad_sum: Quantity-array
            The total radiative solar heat gain through all windows of the zone
            at each hour of the selected design day.
        """
        time = Q_([h * 3600 for h in range(24)], 's')
        if self.windows:
            # Collect all windows of the zone in a single list.
            windows = []
            for windows_dict in self.windows.values():
                windows.extend(list(windows_dict.values()))
            # Initialize a 2D-Numpy array to hold the convective solar heat 
            # gains and a Numpy 2D-array to hold the radiative solar heat gains.
            Q_sol_conv_arr = Q_(np.zeros((len(windows), len(time))), Units.unit_Q)
            Q_sol_rad_arr = Q_(np.zeros((len(windows), len(time))), Units.unit_Q)
            # Iterate over all windows in the zone. For each hour of the design 
            # day calculate the convective and radiative solar heat gain through
            # each window.
            for i, window in enumerate(windows):
                for j, t in enumerate(time):
                    _, Q_sol_conv, Q_sol_rad = window.solar_heat_gain(t)
                    Q_sol_conv_arr[i, j] = Q_sol_conv.to(Units.unit_Q)
                    Q_sol_rad_arr[i, j] = Q_sol_rad.to(Units.unit_Q)
            # At each hour, get the total convective and radiative solar heat
            # gain through all windows.
            Q_sol_conv_sum = np.sum(Q_sol_conv_arr, axis=0)
            Q_sol_rad_sum = np.sum(Q_sol_rad_arr, axis=0)
        else:
            Q_sol_conv_sum = Q_(np.zeros_like(time), Units.unit_Q)
            Q_sol_rad_sum = Q_(np.zeros_like(time), Units.unit_Q)
        return Q_sol_conv_sum, Q_sol_rad_sum
    
    def _collect_internal_heat_gains(self):
        time = Q_([h * 3600 for h in range(24)], 's')
        if self.int_heat_gains:
            # Collect all internal heat gains of the zone in a single list.
            int_heat_gains = list(self.int_heat_gains.values())
            # Initialize a 2D-Numpy array to hold the convective internal heat 
            # gains and a Numpy 2D-array to hold the radiative internal heat 
            # gains.
            Q_ihg_conv_arr = Q_(np.zeros((len(int_heat_gains), len(time))), Units.unit_Q)
            Q_ihg_rad_arr = Q_(np.zeros((len(int_heat_gains), len(time))), Units.unit_Q)
            Q_ihg_lat_arr = Q_(np.zeros((len(int_heat_gains), len(time))), Units.unit_Q)
            # Iterate over all internal heat gains in the zone. For each hour of
            # the design day calculate the convective and radiative sensible 
            # heat gain and latent heat gain.
            for i, ihg in enumerate(int_heat_gains):
                for j, t in enumerate(time):
                    Q_ihg_conv, Q_ihg_rad, Q_ihg_lat = ihg.Q_dot(t.m)
                    Q_ihg_conv_arr[i, j] = Q_ihg_conv.to(Units.unit_Q)
                    Q_ihg_rad_arr[i, j] = Q_ihg_rad.to(Units.unit_Q)
                    Q_ihg_lat_arr[i, j] = Q_ihg_lat.to(Units.unit_Q)
            # At each hour, get the total convective and radiative solar heat
            # gain through all windows.
            Q_ihg_conv_sum = np.sum(Q_ihg_conv_arr, axis=0)
            Q_ihg_rad_sum = np.sum(Q_ihg_rad_arr, axis=0)
            Q_ihg_lat_sum = np.sum(Q_ihg_lat_arr, axis=0)
        else:
            Q_ihg_conv_sum = Q_(np.zeros_like(time), Units.unit_Q)
            Q_ihg_rad_sum = Q_(np.zeros_like(time), Units.unit_Q)
            Q_ihg_lat_sum = Q_(np.zeros_like(time), Units.unit_Q)
        return Q_ihg_conv_sum, Q_ihg_rad_sum, Q_ihg_lat_sum
    
    def _collect_ventilation_latent_heat_gains(self):
        time = Q_([h * 3600 for h in range(24)], 's')
        if self.ventilation:
            Q_vent_lat = Quantity.from_list([self.ventilation.Q_dot_lat(t.m) for t in time])
            Q_vent_lat = Q_vent_lat.to(Units.unit_Q)
        else:
            Q_vent_lat = Q_(np.zeros_like(time), Units.unit_Q)
        return Q_vent_lat
    
    def _collect_radiative_heat_gains(self, num_cycles: int) -> dict[str, Quantity]:
        """Collects all the radiative heat gains in the zone for each hour of
        the selected design day. This daily 24h cycle is then repeated 
        `num_cycles` times.

        Returns
        -------
        Q_gain_rad_dict: dict[str, Quantity]
            Dictionary with the input name of the radiative heat gain which
            maps to a `Quantity` array with the hourly radiative heat gain
            values over `num_cycles` days.
        """
        input_name = self.get_input_name('Q_gain_rad')
        Q_gain_rad = self.Q_gain_sol['rad'] + self.Q_gain_ihg['rad']
        Q_gain_rad = np.tile(Q_gain_rad, num_cycles)
        return {input_name: Q_gain_rad}

    def _collect_convective_heat_gains(self, num_cycles: int) -> dict[str, Quantity]:
        """Collects all the convective heat gains in the zone for each hour of
        the selected design day. This daily 24h cycle is then repeated 
        `num_cycles` times.

        Returns
        -------
        Q_gain_conv_dict: dict[str, Quantity]
            Dictionary with the input name of the convective heat gain which
            maps to a `Quantity` array with the hourly convective heat gain
            values over `num_cycles` days.
        """
        input_name = self.get_input_name('Q_gain_conv')
        Q_gain_conv = self.Q_gain_sol['conv'] + self.Q_gain_ihg['conv']
        Q_gain_conv = np.tile(Q_gain_conv, num_cycles)
        return {input_name: Q_gain_conv}
  
    def _collect_sol_air_temperatures(self, num_cycles: int) -> dict[str, Quantity]:
        """Get the sol-air temperature at each exterior building element 
        surrounding the zone at each hour of the selected design day. This 24h 
        cycle is then repeated `num_cycles` times.
        
        Returns
        -------
        T_sa_dict: dict[str, Quantity]
            Dictionary with the input names for the sol-air temperatures which
            map to `Quantity` arrays with the hourly temperature values repeated
            for `num_cycles` days.             
        """
        time = [h * 3600 for h in range(24)]
        d = {}
        for ext_build_elem in self.ext_build_elems.values():
            input_name = ext_build_elem.get_input_name()
            T_sa = Quantity.from_list([ext_build_elem.ext_surf.T_sa(t) for t in time])
            T_sa = T_sa.to(Units.unit_T)
            T_sa = np.tile(T_sa, num_cycles)
            d[input_name] = T_sa
            if ext_build_elem.doors:
                ext_doors = list(ext_build_elem.doors.values())
                for ext_door in ext_doors:
                    input_name = ext_door.get_input_name()
                    T_sa = Quantity.from_list([ext_door.ext_surf.T_sa(t) for t in time])
                    T_sa = T_sa.to(Units.unit_T)
                    T_sa = np.tile(T_sa, num_cycles)
                    d[input_name] = T_sa
        return d
    
    def _collect_adjacent_zone_temperatures(self, num_cycles: int) -> dict[str, Quantity]:
        time = [h * 3600 for h in range(24)]
        d = {}
        for int_build_elem in self.int_build_elems.values():
            input_name = int_build_elem.get_input_name()
            T_za_adj = Quantity.from_list([int_build_elem.T_za_adj for _ in time])
            T_za_adj = T_za_adj.to(Units.unit_T)
            T_za_adj = np.tile(T_za_adj, num_cycles)
            d[input_name] = T_za_adj
            # interior doors share the same input with the interior building
            # element they belong to.
        return d
    
    def _collect_outdoor_air_temperatures(self, num_cycles: int) -> dict[str, Quantity]:
        """Get the dry-bulb outdoor air temperature at each hour of the selected
        design day. This daily 24h cycle is then repeated `num_cycles` times. 
        
        Returns
        -------
        T_ext_db_dict: dict[str, Quantity]
            Dictionary with the input name of the dry-bulb outdoor air 
            temperature which maps to a `Quantity` array with the hourly 
            temperature values over `num_cycles` days. 
        """
        input_name = self.get_input_name('T_ext_db')
        time = [h * 3600 for h in range(24)]
        T_ext_db = Quantity.from_list([self.weather_data.T_db(t) for t in time])
        T_ext_db = T_ext_db.to(Units.unit_T)
        T_ext_db = np.tile(T_ext_db, num_cycles)
        return {input_name: T_ext_db}
    
    def _collect_supply_air_temperatures(self, num_cycles: int) -> dict[str, Quantity]:
        """Get the ventilation supply air temperature at each hour of the 
        selected design day. This daily 24h cycle is then repeated `num_cycles` 
        times.
        
        Returns
        -------
        T_sup_dict: dict[str, Quantity]
            Dictionary with the input name of the ventilation supply air 
            temperature which maps to a `Quantity` array with the hourly 
            temperature values over `num_cycles` days.
        """
        input_name = self.get_input_name('T_sup')
        time = [h * 3600 for h in range(24)]
        T_sup = Quantity.from_list([self.ventilation.T_sup(t) for t in time])
        T_sup = T_sup.to(Units.unit_T)
        T_sup = np.tile(T_sup, num_cycles)
        return {input_name: T_sup}
    
    def _collect_setpoint_values(self, setpoint: Quantity, num_cycles: int) -> dict[str, Quantity]:
        input_name = self.get_input_name('T_sp')
        time = [h * 3600 for h in range(24)]
        T_sp = Quantity.from_list([setpoint for _ in time])
        T_sp = T_sp.to(Units.unit_T)
        T_sp = np.tile(T_sp, num_cycles)
        return {input_name: T_sp}
        
    def simulate(
        self,
        setpoint: Quantity,
        K_p: Quantity,
        K_i: Quantity,
        num_cycles: int = 10, 
        units: dict[str, str] | None = None,
        reduced_order: int | None = None,
        init_values: np.ndarray | None = None
    ) -> SimulationResults:
        """Simulates the temperature feedback control of the zone-air 
        temperature on the selected design day. 
                
        Parameters
        ----------
        setpoint:
            The setpoint for the zone-air temperature.
        K_p:
            The proportional gain of the temperature controller.
        K_i:
            The integral gain of the temperature controller.
        num_cycles:
            The number of days the selected design day needs to be repeated 
            before the simulation ends. When the simulation starts, the initial 
            values of all system state variables are set to zero, and the daily,
            24h cycle of system input values needs to be repeated a number of 
            days until the system shows a harmonic periodic behavior in time and 
            the zone-air temperature 'T_za' has settled around the setpoint. 
            However, after the simulation has finished, only the last 24 hours 
            of the simulation will be returned.
        units:
            Dictionary with the principal measuring units of the system 
            parameters. There are only three principal measuring units:
            - 'time': default unit is seconds ('s').
            - 'temperature': default unit is Kelvin ('K').
            - 'energy': default unit is Joule ('J').
            All system parameters, e.g. thermal resistance, thermal capacity, 
            heat flow, etc. have measuring units which are derived from these 
            principal measuring units. For example, when the timescale of the 
            simulation is set to hours (`units={'time': 'hr'}`), all measuring
            units that contain time, will be automatically adjusted, e.g, heat
            flow will then be expressed in units of J/hr. 
        reduced_order:
            The order to which the feedback system should be reduced. If `None`, 
            the order of the feedback system won't be reduced.
            As the system is build from a big number of temperature nodes, the
            system order will also be high. However, without noticeable loss of 
            accuracy, the order can normally be reduced to an order of 8 or 6
            (but not below the number of system outputs; should `reduced_order`
            be less than the number of system outputs, it will be limited to the
            number of system outputs). 
            Note that without system reduction, it may also happen that some 
            internal mathematical operations raise an exception (e.g. 
            `FloatingPointError`).
        init_values: optional
            Initial values of the feedback system's state variables. By default,
            `init_values` is `None`, meaning that all initial values will be set
            to zero at the start of the simulation (t = 0).
            When the simulation has finished, the values of the state variables 
            at the start of the last daily cycle are returned in an attribute
            `init_values` of the `SimulationResults` object. These values can 
            then be used as the initial values for a next simulation. 

        Returns
        -------
        SimulationResults
            - the time values (attribute `time`), 
            - the values of all inputs of the zone's feedback system, 
            - the zone-air temperatures (attribute `T_za`),
            - the interior mass node temperatures (attribute `T_im`), 
            - the values of the error signal to the controller (attribute `error`), 
            - the values of the output command of the controller, i.e. the 
              (sensible) heat rate released by the HVAC-system into the zone 
              (positive values) or extracted from the zone (negative values) 
              (attribute `Q_hvac`),
            - the daily amount of heating energy delivered by the HVAC-system 
              (attribute `E_hvac_heat`),
            - the daily amount of cooling energy delivered by the HVAC-system 
              (attribute `E_hvac_cool`),
            - the peak heating rate of the HVAC-system (attribute 
              `Q_hvac_heat_pk`), 
            - the peak cooling rate of the HVAC-system (attribute 
              `Q_hvac_cool_pk`), and
            - the zone's thermal network model (`ZoneModel` object with all the 
              temperature nodes) (attribute `zone_model`).
            - the zone's feedback system (`control.StateSpace` object) 
              (attribute `feedback_system`).
        """
        # Set the measuring units of the quantities used in the system.
        if units is not None:
            if 'time' in units.keys():
                Units.set_time_unit(units['time'])
            if 'temperature' in units.keys():
                Units.set_temperature_unit(units['temperature'])
            if 'energy' in units.keys():
                Units.set_energy_unit(units['energy'])
        # Create the zone model. Get the time axis and collect the input data to
        # drive the system. Get the control feedback system of the zone.
        zone_model = self.create_model()
        time, input_data = self.collect_input_data(num_cycles, setpoint)
        feedback_system = zone_model.create_feedback_system(K_p, K_i, reduced_order)
        # Arrange the system inputs in a Numpy array in the correct order for 
        # passing them to `ct.forced_response()`
        inputs = np.zeros((len(feedback_system.input_labels), len(time)))
        for input_name, input_values in input_data.items():
            try:
                input_index = feedback_system.input_index[input_name]
            except KeyError as err:
                raise err
            inputs[input_index, :] = input_values.magnitude
        # Run the simulation (solve the system).
        if init_values is None:
            resp = ct.forced_response(feedback_system, time, inputs)
        else:
            resp = ct.forced_response(feedback_system, time, inputs, init_values)
        # Collect time and system outputs (zone air temperature, interior mass
        # node temperature, interior surface temperatures, HVAC-system heat
        # transfer).
        self.time = Q_([h for h in range(24)], 'hr').to(Units.unit_t)
        T_za = Q_(resp.outputs[0, -24:].flatten(), Units.unit_T)
        T_im = Q_(resp.outputs[1, -24:].flatten(), Units.unit_T)
        T_is_dict = {}
        for i, output_label in enumerate(resp.output_labels):
            if output_label.startswith('T_is'):
                T_is_dict[output_label] = Q_(resp.outputs[i, -24:].flatten(), Units.unit_T)
        Q_hvac = Q_(resp.outputs[-2, -24:].flatten(), Units.unit_Q)
        # Calculate the thermal energy delivered by the HVAC-system during the
        # day and the peak rates of heating and/or cooling.
        res = self._analyze_hvac_energy(self.time, Q_hvac)
        E_hvac_heat = res[0]
        E_hvac_cool = res[1]
        Q_hvac_heat_pk = res[2]
        Q_hvac_cool_pk = res[3]
        # Only keep the last 24 hours of the input values.
        input_data = {k: v[-24:] for k, v in input_data.items()}
        # Collect conduction heat gains and convective ventilation heat gains.
        self.Q_gain_cond = self.calc_conduction_gains(self.time, T_za, T_im, T_is_dict, zone_model)
        Q_inf, Q_sup, Q_trf = self.calc_vent_heat_gains(self.time, T_za)
        Q_vent_conv = np.sum(np.array([
            Q_inf.to(Units.unit_Q).magnitude,
            Q_sup.to(Units.unit_Q).magnitude,
            Q_trf.to(Units.unit_Q).magnitude
        ]), axis=0)
        self.Q_gain_vent['conv'] = Q_(Q_vent_conv, Units.unit_Q)
        # Calculate heat flow between interior mass node and zone-air node.
        self.Q_gain_im = self._calc_interior_heat_transfer(T_za, T_im)
        # Return simulation results.
        sim_result = SimulationResults(
            time=self.time,
            T_za=T_za,
            T_im=T_im,
            T_is_dict=T_is_dict,
            Q_hvac=Q_hvac,
            error=Q_(resp.outputs[-1, -24:].flatten(), Units.unit_T),
            E_hvac_heat=E_hvac_heat,
            E_hvac_cool=E_hvac_cool,
            Q_hvac_heat_pk=Q_hvac_heat_pk,
            Q_hvac_cool_pk=Q_hvac_cool_pk,
            init_values=resp.states[:, -24].flatten(),
            input_data=input_data,
            zone_model=zone_model,
            feedback_system=feedback_system
        )
        return sim_result
    
    @staticmethod
    def _analyze_hvac_energy(time: Quantity, Q_hvac: Quantity) -> tuple[Quantity, ...]:
        """Calculates peak heat rate delivered by the HVAC-system while heating 
        and/or cooling, and the daily amount of thermal energy delivered for 
        heating and/or for cooling. 
        """
        pos_mask = Q_hvac > 0
        neg_mask = Q_hvac < 0
        if np.any(pos_mask):
            Q_hvac_heat = Q_hvac[pos_mask]
            Q_hvac_heat_pk = np.max(Q_hvac_heat)
            E_hvac_heat = integrate.simpson(Q_hvac_heat.m, time[pos_mask].m)
            E_hvac_heat = Q_(E_hvac_heat, Units.unit_E)
        else:
            Q_hvac_heat_pk = Q_(0.0, Units.unit_Q)
            E_hvac_heat = Q_(0.0, Units.unit_E)
        if np.any(neg_mask):
            Q_hvac_cool = Q_hvac[neg_mask]
            Q_hvac_cool_pk = np.min(Q_hvac_cool)
            E_hvac_cool = abs(integrate.simpson(Q_hvac_cool.m, time[neg_mask].m))
            E_hvac_cool = Q_(E_hvac_cool, Units.unit_E)
        else:
            Q_hvac_cool_pk = Q_(0.0, Units.unit_Q)
            E_hvac_cool = Q_(0.0, Units.unit_E)
        return E_hvac_heat, E_hvac_cool, Q_hvac_heat_pk, Q_hvac_cool_pk
    
    def calc_conduction_gains(
        self,
        time: Quantity,
        T_za: Quantity,
        T_im: Quantity,
        T_is_dict: Quantity,
        zone_model: ThermalZoneModel
    ) -> dict[str, Quantity]:
        """Calculates the hourly values of convective and radiative conduction
        heat gain. The returned conduction heat gains are composed of:
        - conduction heat gain through exterior building elements,
        - conduction heat gain through windows,
        - conduction heat gain through exterior doors,
        - conduction heat gain through interior building elements, and
        - conduction heat gain through interior doors.
        
        Parameters
        ----------
        time:
            Time moments (the hours of the day) for which the heat gains are
            calculated. Available through the returned `SimulationResults` 
            data-object when method `simulate(...)` has been called.
        T_za:
            The hourly values of the zone air temperature, which are available
            after a call to method `simulate(...)` through the returned 
            `SimulationResults` data-object. 
        T_im:
            The hourly values of the interior mass temperature, which are 
            available after a call to method `simulate(...)` through the 
            returned `SimulationResults` data-object.
        T_is_dict:
            Dictionary of which the keys are the names of the interior surface
            nodes of the exterior building elements. These keys map to the 
            hourly values of the interior surface temperatures. This dictionary
            is available after a call to method `simulate(...)` through the 
            returned `SimulationResults` data-object.
        zone_model:
            The thermal network model of the zone, which is available after a 
            call to method `simulate(...)` through the returned 
            `SimulationResults` data-object.
        
        Returns
        -------
        dict[str, Quantity]
            The dictionary has two keys:
            - Key 'conv' maps to the hourly values of convective conduction heat
              gain to the zone.
            - Key 'rad' maps to the hourly values of radiative conduction heat
              gain to the zone.
        """
        ebe_cond_gains = self.calc_ebe_cond_gains(T_za, T_im, T_is_dict, zone_model)
        wnd_cond_gains = self.calc_wnd_cond_gains(time, T_za, T_im)
        edr_cond_gains = self.calc_ext_door_cond_gains(time, T_za, T_im)
        ibe_cond_gains = self.calc_ibe_cond_gains(time, T_za, T_im)
        idr_cond_gains = self.calc_int_door_cond_gains(time, T_za, T_im)
        Q_cnd_conv, Q_cnd_rad = [], []
        for rec in ebe_cond_gains.values():
            Q_cnd_conv.append(rec['conv'])
            Q_cnd_rad.append(rec['rad'])
        for wall_key in wnd_cond_gains.keys():
            for rec in wnd_cond_gains[wall_key].values():
                Q_cnd_conv.append(rec['conv'])
                Q_cnd_rad.append(rec['rad'])
        for wall_key in edr_cond_gains.keys():
            for rec in edr_cond_gains[wall_key].values():
                Q_cnd_conv.append(rec['conv'])
                Q_cnd_rad.append(rec['rad'])
        for rec in ibe_cond_gains.values():
            Q_cnd_conv.append(rec['conv'])
            Q_cnd_rad.append(rec['rad'])
        for wall_key in idr_cond_gains.keys():
            for rec in idr_cond_gains[wall_key].values():
                Q_cnd_conv.append(rec['conv'])
                Q_cnd_rad.append(rec['rad'])
        Q_cnd_conv = np.array([q.to(Units.unit_Q).m for q in Q_cnd_conv])
        Q_cnd_rad = np.array([q.to(Units.unit_Q).m for q in Q_cnd_rad])
        Q_cnd_conv_sum = np.sum(Q_cnd_conv, axis=0)
        Q_cnd_rad_sum = np.sum(Q_cnd_rad, axis=0)
        Q_gain_cond = {
            'conv': Q_(Q_cnd_conv_sum, Units.unit_Q), 
            'rad': Q_(Q_cnd_rad_sum, Units.unit_Q)
        }
        return Q_gain_cond
    
    @staticmethod
    def calc_ebe_cond_gains(
        T_za: Quantity,
        T_im: Quantity,
        T_is_dict: Quantity,
        zone_model: ThermalZoneModel
    ) -> dict[str, dict[str, Quantity]]:
        """Calculates the hourly values of convective and radiative conduction
        heat gains through exterior building elements.
        
        Parameters
        ----------
        T_za:
            The hourly values of the zone air temperature, which are available
            after a call to method `simulate(...)` through the returned 
            `SimulationResults` data-object. 
        T_im:
            The hourly values of the interior mass temperature, which are 
            available after a call to method `simulate(...)` through the 
            returned `SimulationResults` data-object.
        T_is_dict:
            Dictionary of which the keys are the names of the interior surface
            nodes of the exterior building elements. These keys map to the 
            hourly values of the interior surface temperatures. This dictionary
            is available after a call to method `simulate(...)` through the 
            returned `SimulationResults` data-object.
        zone_model:
            The thermal network model of the zone, which is available after a 
            call to method `simulate(...)` through the returned 
            `SimulationResults` data-object.
        
        Returns
        -------
        dict[str, dict[str, Quantity]]
            The keys of the dictionary are the names of the exterior building
            elements surrounding the zone. The values of these keys are a 
            dictionary also, having two keys:
            - Key 'conv' maps to the hourly values of convective conduction heat
              gain through the exterior building element.
            - Key 'rad' maps to the hourly values of radiative conduction heat
              gain through the exterior building element.  
        """
        ext_build_elem_gains = {}
        zan_dict = {
            tup[0].name: tup[1] 
            for tup in zone_model.zone_air_node.resistors_in
            if isinstance(tup[0], TemperatureNode)
        }
        imn_dict = {
            tup[0].name: tup[1] 
            for tup in zone_model.int_mass_node.resistors_in
            if isinstance(tup[0], TemperatureNode)
        }
        for key, T_is in T_is_dict.items():
            R_conv = zan_dict[key]
            Q_cnd_conv = (T_is - T_za) / R_conv.value
            R_rad = imn_dict[key]
            Q_cnd_rad = (T_is - T_im) / R_rad.value
            record = {'conv': Q_cnd_conv, 'rad': Q_cnd_rad}
            _, wall_name = tuple(key.split('@'))
            ext_build_elem_gains[wall_name] = record
        return ext_build_elem_gains
        
    def calc_wnd_cond_gains(
        self, 
        time: Quantity,
        T_za: Quantity,
        T_im: Quantity
    ) -> dict[str, dict[str, dict[str, Quantity]]]:
        """Calculates the hourly values of convective and radiative conduction
        heat gains through windows.
        
        Parameters
        ----------
        time:
            Time moments (the hours of the day) for which the heat gains are
            calculated. Available through the returned `SimulationResults` 
            data-object when method `simulate(...)` has been called.
        T_za:
            The hourly values of the zone air temperature, which are available
            after a call to method `simulate(...)` through the returned 
            `SimulationResults` data-object. 
        T_im:
            The hourly values of the interior mass temperature, which are 
            available after a call to method `simulate(...)` through the 
            returned `SimulationResults` data-object.
        
        Returns
        -------
        dict[str, dict[str, dict[str, Quantity]]]
            The keys of the main dictionary are the names of exterior building
            elements that have windows. These keys map to another dictionary of
            which the keys are the names of the windows. The values of these
            keys are dictionaries that have two keys:
            - Key 'conv' maps to the hourly values of the convective conduction
              heat gain through the window.
            - Key 'rad' maps to the hourly values of the radiative conduction
              heat gain through the window.
        """
        T_ext_db = Quantity.from_list([
            self.weather_data.T_db(t.to('s').m) 
            for t in time
        ]).to('K')
        wnd_cond_gains = {}
        for wall_name, windows_dict in self.windows.items():
            wnd_cond_gains[wall_name] = {}
            for window_name, window in windows_dict.items():
                Q_cnd_conv = (1 - window.F_rad) * window.UA * (T_ext_db - T_za)
                Q_cnd_rad = window.F_rad * window.UA * (T_ext_db - T_im)
                record = {'conv': Q_cnd_conv, 'rad': Q_cnd_rad}
                wnd_cond_gains[wall_name][window_name] = record
        return wnd_cond_gains

    def calc_ext_door_cond_gains(
        self,
        time: Quantity,
        T_za: Quantity,
        T_im: Quantity
    ) -> dict[str, dict[str, dict[str, Quantity]]]:
        """Calculates the hourly values of convective and radiative conduction
        heat gains through exterior doors.
        
        Parameters
        ----------
        time:
            Time moments (the hours of the day) for which the heat gains are
            calculated. Available through the returned `SimulationResults` 
            data-object when method `simulate(...)` has been called.
        T_za:
            The hourly values of the zone air temperature, which are available
            after a call to method `simulate(...)` through the returned 
            `SimulationResults` data-object. 
        T_im:
            The hourly values of the interior mass temperature, which are 
            available after a call to method `simulate(...)` through the 
            returned `SimulationResults` data-object.
            
        Returns
        -------
        dict[str, dict[str, dict[str, Quantity]]]
            The keys of the main dictionary are the names of exterior building
            elements that have exterior doors. These keys map to another 
            dictionary of which the keys are the names of the exterior doors. 
            The values of these keys are dictionaries that have two keys:
            - Key 'conv' maps to the hourly values of the convective conduction
              heat gain through the exterior door.
            - Key 'rad' maps to the hourly values of the radiative conduction
              heat gain through the exterior door.
        """
        ext_door_gains = {}
        for wall_name, ext_door_dict in self.ext_doors.items():
            ext_door_gains[wall_name] = {}
            for ext_door_name, ext_door in ext_door_dict.items():
                T_sa = Quantity.from_list([
                    ext_door.ext_surf.T_sa(t.to('s').m)
                    for t in time
                ]).to('K')
                Q_cnd_conv = (1 - ext_door.F_rad) * ext_door.UA * (T_sa - T_za)
                Q_cnd_rad = ext_door.F_rad * ext_door.UA * (T_sa - T_im)
                record = {'conv': Q_cnd_conv, 'rad': Q_cnd_rad}
                ext_door_gains[wall_name][ext_door_name] = record
        return ext_door_gains
    
    def calc_ibe_cond_gains(
        self,
        time: Quantity,
        T_za: Quantity,
        T_im: Quantity
    ) -> dict[str, dict[str, Quantity]]:
        """Calculates the hourly values of convective and radiative conduction
        heat gains through interior building elements.
        
        Parameters
        ----------
        time:
            Time moments (the hours of the day) for which the heat gains are
            calculated. Available through the returned `SimulationResults` 
            data-object when method `simulate(...)` has been called.
        T_za:
            The hourly values of the zone air temperature, which are available
            after a call to method `simulate(...)` through the returned 
            `SimulationResults` data-object. 
        T_im:
            The hourly values of the interior mass temperature, which are 
            available after a call to method `simulate(...)` through the 
            returned `SimulationResults` data-object.
            
        Returns
        -------
        dict[str, dict[str, Quantity]]
            The keys of the dictionary are the names of the interior building
            elements surrounding the zone. The values of these keys are a 
            dictionary also, having two keys:
            - Key 'conv' maps to the hourly values of convective conduction heat
              gain through the interior building element.
            - Key 'rad' maps to the hourly values of radiative conduction heat
              gain through the interior building element.  
        """
        int_build_elem_gains = {}
        for wall_name, int_build_elem in self.int_build_elems.items():
            T_za_adj = int_build_elem.T_za_adj
            T_za_adj = Quantity.from_list([T_za_adj] * len(time)).to('K')
            Q_cnd_conv = (1 - int_build_elem.F_rad) * int_build_elem.UA * (T_za_adj - T_za)
            Q_cnd_rad = int_build_elem.F_rad * int_build_elem.UA * (T_za_adj - T_im)
            record = {'conv': Q_cnd_conv, 'rad': Q_cnd_rad}
            int_build_elem_gains[wall_name] = record
        return int_build_elem_gains
    
    def calc_int_door_cond_gains(
        self,
        time: Quantity,
        T_za: Quantity,
        T_im: Quantity
    ) -> dict[str, dict[str, dict[str, Quantity]]]:
        """Calculates the hourly values of convective and radiative conduction
        heat gains through interior doors.
        
        Parameters
        ----------
        time:
            Time moments (the hours of the day) for which the heat gains are
            calculated. Available through the returned `SimulationResults` 
            data-object when method `simulate(...)` has been called.
        T_za:
            The hourly values of the zone air temperature, which are available
            after a call to method `simulate(...)` through the returned 
            `SimulationResults` data-object. 
        T_im:
            The hourly values of the interior mass temperature, which are 
            available after a call to method `simulate(...)` through the 
            returned `SimulationResults` data-object.
            
        Returns
        -------
        dict[str, dict[str, dict[str, Quantity]]]
            The keys of the main dictionary are the names of interior building
            elements that have interior doors. These keys map to another 
            dictionary of which the keys are the names of the interior doors. 
            The values of these keys are dictionaries that have two keys:
            - Key 'conv' maps to the hourly values of the convective conduction
              heat gain through the interior door.
            - Key 'rad' maps to the hourly values of the radiative conduction
              heat gain through the interior door.
        """
        int_door_gains = {}
        for wall_name, int_door_dict in self.int_doors.items():
            int_door_gains[wall_name] = {}
            for int_door_name, int_door in int_door_dict.items():
                T_za_adj = int_door.T_za_adj
                T_za_adj = Quantity.from_list([T_za_adj] * len(time)).to('K')
                Q_cnd_conv = (1 - int_door.F_rad) * int_door.UA * (T_za_adj - T_za)
                Q_cnd_rad = int_door.F_rad * int_door.UA * (T_za_adj - T_im)
                record = {'conv': Q_cnd_conv, 'rad': Q_cnd_rad}
                int_door_gains[wall_name][int_door_name] = record
        return int_door_gains
    
    def calc_vent_heat_gains(
        self,
        time: Quantity,
        T_za: Quantity
    ) -> tuple[Quantity, ...]:
        """Calculates the hourly values of the sensible ventilation heat gains 
        to the zone (sensible ventilation heat gains are only convective, i.e.,
        the complete heat gain is directly transferred to the zone air, and no
        fraction of this heat gain is transferred to the interior mass of the 
        zone).
        
        Parameters
        ----------
        time:
            Time moments (the hours of the day) for which the heat gains are
            calculated.
        T_za:
            The hourly values of the zone air temperature, which are available
            after a call to method `simulate(...)` through the returned 
            `SimulationResults` data-object. 
        
        Returns
        -------
        Q_inf:
            Hourly values of ventilation heat gain due to infiltration of 
            outdoor air to the zone.
        Q_sup:
            Hourly values of ventilation heat gain due to supply of ventilation
            air to the zone.
        Q_trf:
            Hourly values of ventilation heat gain due to transfer of air from
            adjacent zones.   
        """
        air = Air(T=self.air_temperature, P=self.air_pressure)
        if self.ventilation.V_dot_ext > 0.0:
            R_inf = 1 / (air.rho * air.cp * self.ventilation.V_dot_ext)
            T_ext_db = Quantity.from_list([
                self.weather_data.T_db(t.to('s').m)
                for t in time
            ]).to('K')
            Q_inf = (T_ext_db - T_za) / R_inf
        else:
            Q_inf = Q_(np.zeros_like(time), Units.unit_Q)
        if self.ventilation.V_dot_sup > 0.0:
            R_sup = 1 / (air.rho * air.cp * self.ventilation.V_dot_sup)
            T_sup = Quantity.from_list([
                self.ventilation.T_sup(t.to('s').m) 
                for t in time
            ]).to('K')
            Q_sup = (T_sup - T_za) / R_sup
        else:
            Q_sup = Q_(np.zeros_like(time), Units.unit_Q)
        Q_trf_list = []        
        for int_build_elem in self.int_build_elems.values():
            if int_build_elem.V_dot_trf is not None and int_build_elem.V_dot_trf > 0.0:
                R_trf = 1 / (air.rho * air.cp * int_build_elem.V_dot_trf)
                T_za_adj = Quantity.from_list([int_build_elem.T_za_adj] * len(time)).to('K')
                Q_trf = (T_za_adj - T_za) / R_trf
                Q_trf_list.append(Q_trf)
        if Q_trf_list:
            Q_trf = np.sum(np.array([Q_trf_list]), axis=0)
        else:
            Q_trf = Q_(np.zeros_like(time), Units.unit_Q)
        return Q_inf, Q_sup, Q_trf
    
    def _calc_interior_heat_transfer(self, T_za: Quantity, T_im: Quantity) -> Quantity:
        """Calculates the heat flow between the interior mass node and the zone
        air node. The heat flow is positive if the interior mass temperature is
        higher than the zone-air temperature.
        """
        R_im = self.unit_R_im / (self.factor_surf_im * self.floor_area)
        Q_gain_im = (T_im - T_za) / R_im
        return Q_gain_im
    
    def heat_gains_overview(self, units: str = 'W') -> pd.DataFrame:
        """Returns a Pandas DataFrame object with a summarized overview of the
        hourly values of heat gains to the zone.
        
        Parameters
        ----------
        units:
            The wanted measuring units of the heat gains. Default value is
            'W' (Watts).
        
        Returns
        -------
        Pandas DataFrame
            The data frame has the following columns:
            - 'Q_cond' ('conv' and 'rad'): 
              The hourly values of convective/radiative conduction heat gain to 
              the zone through exterior building elements, windows, exterior 
              doors, interior building elements, and interior doors.
            - 'Q_sol' ('conv' and 'rad'):
              The hourly values of convective/radiative solar heat gain through 
              windows.
            - 'Q_ihg' ('conv', 'rad', and 'lat'):
              The hourly values of convective/radiative/latent internal heat 
              gains in the zone from equipment, lighting, and people.
            - 'Q_vent' ('conv' and 'lat'):
              The hourly values of convective/latent ventilation heat gains in 
              the zone from outdoor air infiltration, supply of ventilation air 
              to the zone, and transfer of air to the zone from adjacent zones.
            - 'Q_im' ('conv'):
              The hourly values of convective heat transfer from the interior
              mass to the zone air.
        """
        if self.time is not None:
            tuples = [
                ('Q_cond', 'conv'), ('Q_cond', 'rad'),
                ('Q_sol', 'conv'), ('Q_sol', 'rad'),
                ('Q_ihg', 'conv'), ('Q_ihg', 'rad'), ('Q_ihg', 'lat'),
                ('Q_vent', 'conv'), ('Q_vent', 'lat'),
                ('Q_im', 'conv')
            ]
            columns = pd.MultiIndex.from_tuples(tuples, names=('heat gain', 'type'))
            index = self.time.to('hr').magnitude
            data = np.array([
                self.Q_gain_cond['conv'].to(units).magnitude,
                self.Q_gain_cond['rad'].to(units).magnitude,
                self.Q_gain_sol['conv'].to(units).magnitude,
                self.Q_gain_sol['rad'].to(units).magnitude,
                self.Q_gain_ihg['conv'].to(units).magnitude,
                self.Q_gain_ihg['rad'].to(units).magnitude,
                self.Q_gain_ihg['lat'].to(units).magnitude,
                self.Q_gain_vent['conv'].to(units).magnitude,
                self.Q_gain_vent['lat'].to(units).magnitude,
                self.Q_gain_im.to(units).magnitude
            ])        
            data = np.transpose(data)
            df = pd.DataFrame(data, index, columns)
            return df
        else:
            raise RuntimeError("Call method `simulate` first.")


@dataclass
class SimulationResults:
    time: Quantity
    T_za: Quantity
    T_im: Quantity
    T_is_dict: dict[str, Quantity]
    Q_hvac: Quantity
    error: Quantity
    E_hvac_heat: Quantity
    E_hvac_cool: Quantity
    Q_hvac_heat_pk: Quantity
    Q_hvac_cool_pk: Quantity
    init_values: np.ndarray
    input_data: dict[str, Quantity]
    zone_model: ThermalZoneModel
    feedback_system: ct.StateSpace
