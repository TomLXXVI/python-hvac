from __future__ import annotations
from typing import Callable, Sequence
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from hvac import Quantity
from hvac.cooling_load_calc.core import (
    ExteriorBuildingElement,
    InternalHeatGain,
    InteriorBuildingElement,
    WeatherData,
    Window
)
from hvac.cooling_load_calc.building.conditioned_zone import (
    SpaceVentilation
)
from hvac.cooling_load_calc.building.ventilation_zone import (
    VentilationZone
)
from hvac.cooling_load_calc.core.thermal_models import (
    NodalThermalZoneModel
)

Q_ = Quantity


class UnconditionedZone:
    """Represents a space or a group of spaces in a building of which the zone
    air temperature is free-floating. This class can be used to solve for the
    zone air temperature as a function of time if the heat gains in the zone and
    the heat extraction rate of the cooling system can be calculated and/or are
    given.
    """
    def __init__(self):
        self.ID: str = ''
        self.floor_area: Quantity | None = None
        self.height: Quantity | None = None
        self.ext_build_elems: dict[str, ExteriorBuildingElement] = {}
        self.int_build_elems: dict[str, InteriorBuildingElement] = {}
        self.int_heat_gains: dict[str, InternalHeatGain] = {}
        self.ventilation: SpaceVentilation | None = None
        self.vez: VentilationZone | None = None
        self.tsn: tuple[Quantity, ...] | None = None  # specs of thermal storage node
        self.thz_model: NodalThermalZoneModel | None = None
        self.dt_hr: float = 1.0  # the time step used for solving the nodal thermal zone model

    @classmethod
    def create(
        cls,
        ID: str,
        weather_data: WeatherData,
        floor_area: Quantity,
        height: Quantity,
        ventilation_zone: VentilationZone | None = None,
        C_tsn: Quantity = Q_(100, 'kJ / m**2'),
        A_tsn: Quantity = Q_(1.0, 'm**2'),
        R_tsn: Quantity = Q_(float('inf'), 'K * m**2 / W')
    ) -> UnconditionedZone:
        """Creates a `UnconditionedZone` object.

        Parameters
        ----------
        ID:
            Identifies the zone in the building.
        weather_data:
            Encapsulates the climatic design information needed to determine the
            solar radiation incident on the exterior surface of the building
            element and its sol-air temperature during the design day, which was
            specified on instantiation of the `WeatherData` object.
        floor_area:
            The floor area of the zone.
        height:
            The height of the zone.
        ventilation_zone: optional
            The ventilation zone to which the thermal zone belongs.
        C_tsn:
            The unit thermal capacity of the interior thermal mass.
        A_tsn:
            The surface area associated with the interior thermal mass.
        R_tsn:
            The unit thermal resistance between the interior thermal mass and
            the zone air.
        """
        zone = cls()
        zone.ID = ID
        zone.weather_data = weather_data
        zone.floor_area = floor_area
        zone.height = height
        zone.tsn = (A_tsn, C_tsn, R_tsn)
        zone.vez = ventilation_zone
        if ventilation_zone is not None:
            zone.vez.add_thermal_zone(zone)
        return zone

    def add_ext_build_elem(
        self,
        ebe: ExteriorBuildingElement | Sequence[ExteriorBuildingElement]
    ) -> None:
        """Adds an exterior building element or a sequence of exterior building
        elements to the zone.
        """
        if isinstance(ebe, Sequence):
            self.ext_build_elems.update({ebe_.ID: ebe_ for ebe_ in ebe})
        else:
            self.ext_build_elems[ebe.ID] = ebe

    def add_int_build_elem(
        self,
        ibe: InteriorBuildingElement | Sequence[InteriorBuildingElement]
    ) -> None:
        """Adds an interior building element or a sequence of interior building
        elements to the zone."""
        self.int_build_elems[ibe.ID] = ibe

    def add_ventilation(
        self,
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_dot_open: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_ATD_d: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_sup: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_trf: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_exh: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_comb: Quantity = Q_(0.0, 'm ** 3 / hr'),
        T_sup: Callable[[float], Quantity] | None = None,
        T_trf: Callable[[float], Quantity] | None = None
    ) -> None:
        """Configures the ventilation of the zone according to standard
        EN 12831-1 (2017).

        Parameters
        ----------
        n_min:
            Minimum air change rate required for the space for reasons of air
            quality/hygiene and comfort (EN 12831-1, B.2.10 - Table B.7).
            The default value applies to permanent dwelling areas (living rooms,
            offices) and a ceiling height less than 3 m.
        V_dot_open: Quantity, optional
            External air volume flow rate into the space through large openings
            (EN 12831-1, Annex G).
        V_dot_ATD_d:
            Design air volume flow rate of the ATDs in the room
            (EN 12831-1, B.2.12). Only relevant if ATDs are used for
            ventilation.
        V_dot_sup:
            Supply air volume flow rate from the ventilation system into the
            space.
        V_dot_trf:
            Transfer air volume flow rate into the space from adjacent spaces.
        V_dot_exh:
            Exhaust ventilation air volume flow rate from the space.
        V_dot_comb:
            Air volume flow rate exhausted from the space that has not been
            included in the exhaust air volume flow of the ventilation system
            (typically, but not necessarily, combustion air if an open flue
            heater is present in the heated space).
        T_sup:
            Temperature of the ventilation supply air after passing heat recovery
            (see EN 12831-1 §6.3.3.7) or after passive pre-cooling.
            Any temperature decrease by active pre-cooling, which requires power
            from a cooling machine, shall not be taken into account (i.e., the
            ventilation load is not a space load in that case, but a load on
            the air-handling unit).
        T_trf:
            Temperature of the transfer air that flows into the space from
            another space. In case the room height of the other space is less
            than 4 m, it is equal to the internal design temperature of the other
            space; otherwise, it is equal to mean air temperature of the other
            space (see EN 12831-1, 6.3.8.3).
        """
        self.ventilation = SpaceVentilation.create(
            thz=self,
            n_min=n_min,
            V_dot_open=V_dot_open,
            V_dot_ATD_d=V_dot_ATD_d,
            V_dot_sup=V_dot_sup,
            V_dot_trf=V_dot_trf,
            V_dot_exh=V_dot_exh,
            V_dot_comb=V_dot_comb,
            T_sup=T_sup,
            T_trf=T_trf
        )

    def add_internal_heat_gain(
        self,
        ihg: InternalHeatGain | Sequence[InternalHeatGain]
    ) -> None:
        """Adds an internal heat gain or a sequence of internal heat gains to
        the zone (see cooling_load_calc.core.internal_heat_gains).
        """
        if isinstance(ihg, Sequence):
            self.int_heat_gains.update({ihg_.ID: ihg_ for ihg_ in ihg})
        else:
            self.int_heat_gains[ihg.ID] = ihg

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
        V = self.floor_area * self.height
        return V

    def _solar_heat_gain(
        self,
        dt_hr: float = 1.0
    ) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the solar heat gain through all the windows of the zone,
        and also its radiative and convective components, for each time index
        k of the design day.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step between the calculations is one quarter of an
            hour. The default value is 1 hour.

        Returns
        -------
        3-tuple with:
        - a `Quantity`-array with the solar heat gains at each time index
        - a `Quantity`-array with the convective components
        - a `Quantity`-array with the radiative components
        """
        d = {'Q_dot': [], 'Q_dot_conv': [], 'Q_dot_rad': []}
        for ebe in self.ext_build_elems.values():
            for wnd in ebe.windows.values():
                Q_dot, Q_dot_cv, Q_dot_rd = wnd.solar_heat_gain(dt_hr)
                d['Q_dot'].append(Q_dot)
                d['Q_dot_conv'].append(Q_dot_cv)
                d['Q_dot_rad'].append(Q_dot_rd)
        Q_dot = sum(d['Q_dot'])
        Q_dot_conv = sum(d['Q_dot_conv'])
        Q_dot_rad = sum(d['Q_dot_rad'])
        return Q_dot, Q_dot_conv, Q_dot_rad

    def _internal_heat_gain(
        self,
        dt_hr: float
    ) -> tuple[Quantity, Quantity, Quantity, Quantity]:
        """Returns the internal heat gains in the zone for each time index k of
        the design day.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step between the calculations is one quarter of an
            hour. The default value is 1 hour.

        Returns
        -------
        4-tuple with:
        - a `Quantity`-array with the total of the sensible internal heat gains
          at each time index
        - a `Quantity`-array with the total of the convective components of the
          sensible heat gains
        - a `Quantity`-array with the total of the radiative components of the
          sensible heat gains
        - a `Quantity`-array with the total of the latent heat gains at each
          time index
        """
        num_steps = int(round(24 / dt_hr))
        dt_sec = dt_hr * 3600
        d = {
            'Q_dot_sen': [],
            'Q_dot_sen_cv': [],
            'Q_dot_sen_rd': [],
            'Q_dot_lat': []
        }
        for k in range(num_steps):
            d_k = {
                'Q_dot_sen': [],
                'Q_dot_sen_cv': [],
                'Q_dot_sen_rd': [],
                'Q_dot_lat': []
            }
            for ihg in self.int_heat_gains.values():
                Q_dot_sen_cv, Q_dot_sen_rd, Q_dot_lat = ihg.Q_dot(k * dt_sec)
                Q_dot_sen = Q_dot_sen_cv + Q_dot_sen_rd
                d_k['Q_dot_sen'].append(Q_dot_sen)
                d_k['Q_dot_sen_cv'].append(Q_dot_sen_cv)
                d_k['Q_dot_sen_rd'].append(Q_dot_sen_rd)
                d_k['Q_dot_lat'].append(Q_dot_lat)
            d['Q_dot_sen'].append(sum(d_k['Q_dot_sen']))
            d['Q_dot_sen_cv'].append(sum(d_k['Q_dot_sen_cv']))
            d['Q_dot_sen_rd'].append(sum(d_k['Q_dot_sen_rd']))
            d['Q_dot_lat'].append(sum(d_k['Q_dot_lat']))
        Q_dot_sen = Q_(d['Q_dot_sen'], 'W')
        Q_dot_sen_cv = Q_(d['Q_dot_sen_cv'], 'W')
        Q_dot_sen_rd = Q_(d['Q_dot_sen_rd'], 'W')
        Q_dot_lat = Q_(d['Q_dot_lat'], 'W')
        return Q_dot_sen, Q_dot_sen_cv, Q_dot_sen_rd, Q_dot_lat

    def _get_windows_and_ext_doors(self) -> tuple[list[Window], list[ExteriorBuildingElement]]:
        """Returns a list with all the windows present in the exterior building
        elements surrounding the zone, and also a list with all the exterior
        doors.
        """
        windows = [
            window
            for ebe in self.ext_build_elems.values()
            for window in ebe.windows.values()
        ]
        ext_doors = [
            ext_door
            for ebe in self.ext_build_elems.values()
            for ext_door in ebe.doors.values()
        ]
        return windows, ext_doors

    def _get_int_build_elems(self) -> list[InteriorBuildingElement]:
        """Returns a single list with all the interior building elements
        surrounding the zone and all the interior doors.
        """
        ibe_lst = [ibe for ibe in self.int_build_elems.values()]
        ibe_lst.extend([
            int_door
            for ibe in self.int_build_elems.values()
            for int_door in ibe.doors
        ])
        return ibe_lst

    def _init_values(self) -> list[list[Quantity]]:
        """Sets the initial values of all the node temperatures in the nodal
        thermal model of the zone to the sol-air temperature at midnight (0 s)
        of the design day.
        """
        T_init_lst = []
        T_sa_avg = []
        for ebe in self.ext_build_elems.values():
            T_sa_avg.append(ebe._ext_surf.T_sa(12))
            num_nodes = len(ebe.ltn)
            T_init_sub_lst = [ebe._ext_surf.T_sa(12)] * num_nodes
            T_init_lst.extend(T_init_sub_lst)
        T_sa_avg = Quantity.from_list(T_sa_avg)
        # Add an init-value for the zone air node and thermal storage node,
        # being the average value of the sol-air temperatures of the exterior
        # building elements:
        T_init_lst.extend([sum(T_sa_avg.to('K')) / len(T_sa_avg)] * 2)
        return [T_init_lst, T_init_lst]

    def solve(
        self,
        F_rad: Quantity = Q_(0.46, 'frac'),
        Q_dot_sys_fun: Callable[[float], Quantity] | None = None,
        dt_hr: float = 1.0,
        num_cycles: int = 1
    ) -> None:
        """Creates the nodal thermal zone model of the unconditioned space and
        solves it.
        """
        self.dt_hr = dt_hr
        num_steps = int(round(24 / dt_hr))
        windows, ext_doors = self._get_windows_and_ext_doors()
        int_build_elems = self._get_int_build_elems()
        # Solar heat gain and internal heat gains are independent of the zone
        # air temperature. The radiative and convective components of the
        # solar and internal heat gain are linearly interpolated and these
        # functions will be used in the node equations of the zone air and
        # thermal storage node:
        t_axis = 3600 * np.arange(0, num_steps * dt_hr, dt_hr)  # seconds
        _, Q_dot_sol_cv, Q_dot_sol_rd = self._solar_heat_gain(dt_hr)
        Q_dot_sol_cv_interp = interp1d(t_axis, Q_dot_sol_cv.to('W').m)
        Q_dot_sol_rd_interp = interp1d(t_axis, Q_dot_sol_rd.to('W').m)
        _, Q_dot_ihg_cv, Q_dot_ihg_rd, _ = self._internal_heat_gain(dt_hr)
        Q_dot_ihg_cv_interp = interp1d(t_axis, Q_dot_ihg_cv.to('W').m)
        Q_dot_ihg_rd_interp = interp1d(t_axis, Q_dot_ihg_rd.to('W').m)
        # Create the nodal thermal model of the zone:
        self.thz_model = NodalThermalZoneModel.create(
            ext_build_elems=list(self.ext_build_elems.values()),
            F_rad=F_rad,
            A_tsn=self.tsn[0],
            C_tsn=self.tsn[1],
            R_tsn=self.tsn[2],
            windows=windows,
            ext_doors=ext_doors,
            int_build_elems=int_build_elems,
            ventilation=self.ventilation,
            Q_dot_sol_cv=lambda t_sol_sec: Q_(Q_dot_sol_cv_interp(t_sol_sec), 'W'),
            Q_dot_sol_rd=lambda t_sol_sec: Q_(Q_dot_sol_rd_interp(t_sol_sec), 'W'),
            Q_dot_ihg_cv=lambda t_sol_sec: Q_(Q_dot_ihg_cv_interp(t_sol_sec), 'W'),
            Q_dot_ihg_rd=lambda t_sol_sec: Q_(Q_dot_ihg_rd_interp(t_sol_sec), 'W'),
            Q_dot_sys=Q_dot_sys_fun
        )
        # Solve the model for the node temperatures:
        self.thz_model.solve(
            num_steps=num_steps,
            init_values=None,
            dt_hr=dt_hr,
            num_cycles=num_cycles
        )

    def get_zone_air_temperature(self, unit: str = 'degC') -> Quantity:
        """Returns the zone air temperature in the unconditioned space for
        each solar time index of the design day, after the nodal thermal zone
        model of the unconditioned zone has been solved.
        """
        if self.thz_model is not None:
            return self.thz_model.zone_air_temperature(unit)

    def get_ventilation_heat_gain(self) -> Quantity:
        """Returns the sensible ventilation/infiltration heat gain to the zone
        air at each solar time index of the design day, after the nodal thermal
        zone model of the unconditioned zone has been solved.
        """
        num_steps = int(round(24 / self.dt_hr))
        t_sol_axis = 3600 * np.arange(0, num_steps * self.dt_hr, self.dt_hr)  # seconds
        T_zone_arr = self.get_zone_air_temperature()
        Q_dot_vent_lst = []
        for t_sol_sec, T_zone in zip(t_sol_axis, T_zone_arr):
            T_zone = T_zone.to('degC').m
            Q_dot_vent_ext = (
                0.34 * self.ventilation.V_dot_ext
                * (self.ventilation.T_db_ext(t_sol_sec).to('degC').m
                   - T_zone)
            )
            Q_dot_vent_sup = (
                0.34 * self.ventilation.V_dot_sup
                * (self.ventilation.T_sup(t_sol_sec).to('degC').m
                   - T_zone)
            )
            if self.ventilation.T_trf:
                Q_dot_vent_trf = (
                    0.34 * self.ventilation.V_dot_trf
                    * (self.ventilation.T_trf(t_sol_sec).to('degC').m
                       - T_zone)
                )
            else:
                Q_dot_vent_trf = 0.0
            Q_dot_vent = Q_dot_vent_ext + Q_dot_vent_sup + Q_dot_vent_trf
            Q_dot_vent_lst.append(Q_dot_vent)
        return Q_(Q_dot_vent_lst, 'W')

    def get_conduction_heat_gain(self) -> Quantity:
        """Returns the convective part of the conduction heat gain to the zone
        air through opaque exterior building elements, exterior doors, windows,
        and interior building elements, including interior doors, at each solar
        time index of the design day, after the nodal thermal zone model of the
        unconditioned zone has been solved.
        """
        num_steps = int(round(24 / self.dt_hr))
        t_sol_axis = 3600 * np.arange(0, num_steps * self.dt_hr, self.dt_hr)  # seconds
        T_zone_arr = self.get_zone_air_temperature()
        wnd_lst, edr_lst = self._get_windows_and_ext_doors()
        ibe_lst = self._get_int_build_elems()
        # Convective part of conduction heat gain through windows, exterior
        # doors and interior building elements:
        Q_dot_cnd_cv_lst = []
        for t_sol_sec, T_zone in zip(t_sol_axis, T_zone_arr):
            Q_dot_cnd_cv_lst_ = []  # collects conduction heat gains at t_sol_sec
            for window in wnd_lst:
                Q_dot_wnd_cv = (
                    (1 - window.F_rad) * window.UA
                    * (window._ext_surf.T_db(t_sol_sec).to('K')
                       - T_zone.to('K'))
                )
                Q_dot_cnd_cv_lst_.append(Q_dot_wnd_cv.to('W'))
            for door in edr_lst:
                Q_dot_edr = (
                    (1 - door.F_rad) * door.UA
                    * (door._ext_surf.T_db(t_sol_sec).to('K')
                       - T_zone.to('K'))
                )
                Q_dot_cnd_cv_lst_.append(Q_dot_edr.to('W'))
            for ibe in ibe_lst:
                Q_dot_ibe = (
                    (1 - ibe.F_rad) * ibe.UA
                    * (ibe.T_adj(t_sol_sec).to('K')
                       - T_zone.to('K'))
                )
                Q_dot_cnd_cv_lst_.append(Q_dot_ibe.to('W'))
            Q_dot_cnd_cv_lst.append(sum(Q_dot_cnd_cv_lst_))  # add total of conduction heat gains at t_sol_sec to list
        Q_dot_cnd_cv = Quantity.from_list(Q_dot_cnd_cv_lst)  # convert the list to a `Quantity` array
        # Convective part of conduction heat gain through opaque exterior
        # building elements:
        Q_dot_ebe_dict = self.thz_model.conductive_heat_gain()
        Q_dot_ebe_tot_cv = 0
        # --> total convective conduction heat gain from all exterior building
        # elements together at each time index:
        for key in Q_dot_ebe_dict.keys():
            ebe = self.ext_build_elems[key]
            Q_dot = Q_dot_ebe_dict[key]  # cond. heat gain from ext. build. elem. `ebe` at each time index
            Q_dot_cv = (1 - ebe.F_rad) * Q_dot  # convective part of this conduction heat gain
            Q_dot_ebe_tot_cv += Q_dot_cv  # add to total of all exterior building elements together
        Q_dot_cnd_cv += Q_dot_ebe_tot_cv
        return Q_dot_cnd_cv.to('W')

    def get_solar_heat_gain(self) -> Quantity:
        """Returns the convective part of the solar heat gain to the zone air
        through windows at each solar time index of the design day, after the
        nodal thermal zone model of the unconditioned zone has been solved.
        """
        _, Q_dot_conv, _ = self._solar_heat_gain(self.dt_hr)
        return Q_dot_conv.to('W')

    def get_internal_heat_gain(self) -> Quantity:
        """Returns the convective part of the internal heat gains to the zone
        air at each solar time index of the design day, after the nodal thermal
        zone model of the unconditioned zone has been solved.
        """
        _, Q_dot_sen_cv, *_ = self._internal_heat_gain(self.dt_hr)
        return Q_dot_sen_cv.to('W')

    def get_int_thermal_mass_heat_gain(self) -> Quantity:
        """Returns the convective heat gain from the interior thermal mass to
        the zone air at each solar time index of the design day, after the nodal
        thermal zone model of the unconditioned zone has been solved.
        """
        R_itm_z = self.thz_model.tsn.R2
        # --> unit thermal resistance between interior thermal mass and zone air
        A_itm = self.thz_model.tsn.A[-1]
        # --> surface area associated with the interior thermal mass
        T_tsn = self.thz_model.thermal_storage_temperature('K')
        T_zone = self.thz_model.zone_air_temperature('K')
        Q_dot_itm = A_itm * (T_tsn - T_zone) / R_itm_z
        return Q_dot_itm.to('W')

    def temperature_heat_gain_table(
        self,
        units: dict[str, str] | None = None,
        num_decimals: int = 3
    ) -> pd.DataFrame:
        """Returns a Pandas DataFrame object with the zone air temperature and
        the heat gains to the zone air at each solar time index of the
        design day.

        Parameters
        ----------
        units:
            A dictionary with two keys 'T' and 'Q_dot' of which the values
            are the display units.
        num_decimals:
            The displayed number of decimals of values.
        """
        if units is None: units = {'T': 'degC', 'Q_dot': 'W'}
        T_zone = self.get_zone_air_temperature()
        Q_dot_vent = self.get_ventilation_heat_gain()
        Q_dot_cond = self.get_conduction_heat_gain()
        Q_dot_sol = self.get_solar_heat_gain()
        Q_dot_ihg = self.get_internal_heat_gain()
        Q_dot_itm = self.get_int_thermal_mass_heat_gain()
        Q_dot_net = Q_dot_vent + Q_dot_cond + Q_dot_sol + Q_dot_ihg + Q_dot_itm
        d = {
            'T_zone': T_zone.to(units['T']).magnitude,
            'H.G. vent.': Q_dot_vent.to(units['Q_dot']).magnitude,
            'H.G. cond.': Q_dot_cond.to(units['Q_dot']).magnitude,
            'H.G. solar': Q_dot_sol.to(units['Q_dot']).magnitude,
            'H.G. internal': Q_dot_ihg.to(units['Q_dot']).magnitude,
            'H.G. storage': Q_dot_itm.to(units['Q_dot']).magnitude,
            'H.G. net': Q_dot_net.to(units['Q_dot']).magnitude
        }
        df = pd.DataFrame(d).round(num_decimals)
        return df
