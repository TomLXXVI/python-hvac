from __future__ import annotations
from typing import Callable, Sequence, TYPE_CHECKING
import pandas as pd
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.cooling_load_calc.core import (
    WeatherData,
    ExteriorBuildingElement,
    InteriorBuildingElement,
    ThermalStorageNode,
    InternalHeatGain
)
from hvac.cooling_load_calc.core.utils import convert_to_clock_time

if TYPE_CHECKING:
    from .ventilation_zone import VentilationZone


Q_ = Quantity


class ConditionedZone:
    """Represents a space or a group of spaces in a building of which the
    air temperature is controlled by a single thermostat. Aim is to find the
    required rate at which the cooling system must extract heat to keep the
    zone air temperature at its setpoint value.
    """
    def __init__(self):
        self.ID: str = ''
        self.floor_area: Quantity | None = None
        self.height: Quantity | None = None
        self.weather_data: WeatherData | None = None
        self.T_zone: Callable[[float], Quantity] | None = None
        self.RH_zone: Quantity | None = None
        self.ext_build_elems: dict[str, ExteriorBuildingElement] = {}
        self.int_build_elems: dict[str, InteriorBuildingElement] = {}
        self.tsn: ThermalStorageNode | None = None
        self.ventilation: SpaceVentilation | None = None
        self.vez: VentilationZone | None = None
        self.int_heat_gains: dict[str, InternalHeatGain] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        weather_data: WeatherData,
        T_zone: Callable[[float], Quantity],
        floor_area: Quantity,
        height: Quantity,
        RH_zone: Quantity = Q_(50, 'pct'),
        ventilation_zone: VentilationZone | None = None
    ) -> ConditionedZone:
        """Creates a `ConditionedZone` object.

        Parameters
        ----------
        ID:
            Identifies the zone in the building.
        weather_data:
            Encapsulates the climatic design information needed to determine the
            solar radiation incident on the exterior surface of the building
            element and its sol-air temperature during the design day, which was
            specified on instantiation of the `WeatherData` object.
        T_zone:
            The zone air temperature, being a function with signature
            `f(t_sol_sec: float) -> Quantity` which takes the solar time in
            seconds and returns the temperature in the zone as a `Quantity`
            object. This may allow a time-variable setpoint temperature in the
            zone.
        floor_area:
            The floor area of the zone.
        height:
            The height of the zone.
        RH_zone: optional
            Relative air humidity in the zone. The default is 50 %.
        ventilation_zone: optional
            The ventilation zone to which the thermal zone belongs.
        """
        zone = cls()
        zone.ID = ID
        zone.weather_data = weather_data
        zone.T_zone = T_zone
        zone.floor_area = floor_area
        zone.height = height
        zone.RH_zone = RH_zone
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
        if isinstance(ibe, Sequence):
            self.int_build_elems.update({ibe_.ID: ibe_ for ibe_ in ibe})
        else:
            self.int_build_elems[ibe.ID] = ibe

    def add_thermal_storage_node(
        self,
        C: Quantity = Q_(100, 'kJ / m**2'),
        A: Quantity = Q_(1.0, 'm**2'),
        R_tz: Quantity = Q_(float('inf'), 'K * m**2 / W')
    ) -> None:
        """Adds a thermal storage node to the thermal zone that represents the
        interior thermal mass in the zone.

        Parameters
        ----------
        C:
            The thermal capacity of the node per unit area.
        A:
            The surface area associated with the node.
        R_tz:
            The unit thermal resistance between this node and the zone air
            node.
        """
        self.tsn = ThermalStorageNode.create(
            ID='TSN',
            C=C,
            A=A,
            R2=R_tz
        )

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
    
    def _conductive_heat_gain(
        self,
        dt_hr: float = 1.0,
        num_cycles: int = 5
    ) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the total conductive heat gain in the zone, and also its
        convective and radiative components, at each time index k of the design
        day.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step for the calculations is one quarter of an hour.
            The default value is 1 hour.
        num_cycles:
            Number of diurnal calculation cycles before the results of the last
            diurnal cycle are returned.

        Returns
        -------
        3-tuple with:
        - a `Quantity`-array with the conductive heat gains at each time index
        - a `Quantity`-array with the convective components
        - a `Quantity`-array with the radiative components
        """
        d = {'Q_dot': [], 'Q_dot_conv': [], 'Q_dot_rad': []}
        # Get the conductive heat gain from all exterior building elements,
        # windows and doors:
        for ebe in self.ext_build_elems.values():
            Q_dot, Q_dot_cv, Q_dot_rd = ebe.conductive_heat_gain(dt_hr, num_cycles)
            d['Q_dot'].append(Q_dot)
            d['Q_dot_conv'].append(Q_dot_cv)
            d['Q_dot_rad'].append(Q_dot_rd)
            for wnd in ebe.windows.values():
                Q_dot, Q_dot_cv, Q_dot_rd = wnd.conductive_heat_gain(dt_hr)
                d['Q_dot'].append(Q_dot)
                d['Q_dot_conv'].append(Q_dot_cv)
                d['Q_dot_rad'].append(Q_dot_rd)
            for door in ebe.doors.values():
                Q_dot, Q_dot_cv, Q_dot_rd = door.conductive_heat_gain(dt_hr, num_cycles)
                d['Q_dot'].append(Q_dot)
                d['Q_dot_conv'].append(Q_dot_cv)
                d['Q_dot_rad'].append(Q_dot_rd)
        # Get the conductive heat gain from all interior building elements and
        # doors:
        for ibe in self.int_build_elems.values():
            Q_dot, Q_dot_cv, Q_dot_rd = ibe.conductive_heat_gain(dt_hr)
            d['Q_dot'].append(Q_dot)
            d['Q_dot_conv'].append(Q_dot_cv)
            d['Q_dot_rad'].append(Q_dot_rd)
            for door in ibe.doors.values():
                Q_dot, Q_dot_cv, Q_dot_rd = door.conductive_heat_gain(dt_hr)
                d['Q_dot'].append(Q_dot)
                d['Q_dot_conv'].append(Q_dot_cv)
                d['Q_dot_rad'].append(Q_dot_rd)
        # Sum the collected conductive heat gains and their convective and
        # radiative components:
        Q_dot_cond = sum(d['Q_dot'])
        Q_dot_conv = sum(d['Q_dot_conv'])
        Q_dot_rad = sum(d['Q_dot_rad'])
        return Q_dot_cond, Q_dot_conv, Q_dot_rad

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
        Q_dot_cond = sum(d['Q_dot'])
        Q_dot_conv = sum(d['Q_dot_conv'])
        Q_dot_rad = sum(d['Q_dot_rad'])
        return Q_dot_cond, Q_dot_conv, Q_dot_rad

    def _ventilation_heat_gain(
        self,
        dt_hr: float = 1.0
    ) -> tuple[Quantity, Quantity]:
        """Returns the sensible and the latent heat gain due to space
        ventilation for each time index k of the design day.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step between the calculations is one quarter of an
            hour. The default value is 1 hour.

        Returns
        -------
        A 2-tuple:
        - a `Quantity`-array with the sensible heat gains at each time index
        - a `Quantity`-array with the latent heat gains
        """
        num_steps = int(round(24 / dt_hr))
        dt_sec = dt_hr * 3600
        if self.ventilation is not None:
            Q_dot_sen_vent = Q_([
                self.ventilation.Q_dot_sen(k * dt_sec)
                for k in range(num_steps)
            ], 'W')
            Q_dot_lat_vent = Q_([
                self.ventilation.Q_dot_lat(k * dt_sec)
                for k in range(num_steps)
            ], 'W')
            return Q_dot_sen_vent, Q_dot_lat_vent
        return Q_([0.0] * num_steps, 'W'), Q_([0.0] * num_steps, 'W')

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

    def solve(
        self,
        dt_hr: float = 1.0,
        num_cycles: int = 5,
        unit: str = 'W'
    ) -> pd.DataFrame:
        """Calculates the sensible and latent cooling load of the zone.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step for the calculations is one quarter of an hour.
            The default value is 1 hour.
        num_cycles:
            Number of diurnal calculation cycles before the results of the last
            diurnal cycle are returned.
        unit:
            The measuring unit in which heat flows need to be expressed. The
            default unit is Watts (W).

        Returns
        -------
        A Pandas DataFrame with the heat gains and the total zone load for each
        time index k of the design day.
        """
        # CONDUCTION HEAT GAIN
        # Get the total conductive heat gain in the zone, and its convective and
        # radiative components:
        Q_dot_cnd, Q_dot_cnd_cv, Q_dot_cnd_rd = self._conductive_heat_gain(dt_hr, num_cycles)

        # SOLAR HEAT GAIN
        # Get the total solar heat gain in the zone, and its convective and
        # radiative components:
        Q_dot_sol, Q_dot_sol_cv, Q_dot_sol_rd = self._solar_heat_gain(dt_hr)

        # INTERNAL HEAT GAINS
        Q_dot_ihg = self._internal_heat_gain(dt_hr)
        Q_dot_ihg_sen = Q_dot_ihg[0]
        Q_dot_ihg_sen_cv = Q_dot_ihg[1]
        Q_dot_ihg_sen_rd = Q_dot_ihg[2]
        Q_dot_ihg_lat = Q_dot_ihg[3]

        # CONVECTIVE HEAT GAIN FROM INTERIOR THERMAL MASS
        # If the interior thermal mass of the zone has been defined, take the
        # total radiative heat gain which flows into the thermal mass and get
        # the convective heat flow from the thermal mass to the zone air:
        Q_dot_rad = Q_dot_cnd_rd + Q_dot_sol_rd + Q_dot_ihg_sen_rd
        if self.tsn is not None:
            dt_sec = dt_hr * 3600
            self.tsn.solve(
                num_steps=int(round(24 / dt_hr)),
                Q_dot_rad=lambda t_sol_sec: Q_dot_rad[int(t_sol_sec / dt_sec)],
                T_zone=self.T_zone,
                init_values=[self.T_zone(0)] * 2,
                dt_hr=dt_hr,
                num_cycles=num_cycles
            )
            # T_tsn = self.tsn.get_node_temperatures()
            Q_dot_tsn = self.tsn.get_heat_flows()
        else:
            Q_dot_tsn = None

        # VENTILATION AND INFILTRATION HEAT GAIN
        # Get the sensible and latent heat gain due to ventilation:
        Q_dot_sen_vent, Q_dot_lat_vent = self._ventilation_heat_gain(dt_hr)

        # SENSIBLE ZONE COOLING LOAD
        # The sensible cooling load of the zone = total convective heat gain to
        # the zone. If the interior thermal mass of the zone has been defined,
        # add the convective heat flow coming from the thermal mass to the
        # sensible zone load. If the interior thermal mass of the zone has not
        # been defined, add the total radiative heat gain directly to the
        # sensible zone load:
        Q_dot_sen_zone = Q_dot_cnd_cv + Q_dot_sol_cv + Q_dot_ihg_sen_cv
        if Q_dot_tsn is not None:
            Q_dot_sen_zone += Q_dot_tsn
        else:
            Q_dot_sen_zone += Q_dot_rad
        # Add the sensible ventilation heat gain to the sensible cooling load
        # of the zone:
        Q_dot_sen_zone += Q_dot_sen_vent

        # LATENT ZONE COOLING LOAD
        # Latent cooling load of the zone = latent ventilation gain + latent
        # internal heat gains:
        Q_dot_lat_zone = Q_dot_lat_vent + Q_dot_ihg_lat

        # Return the heat gains in a dataframe:
        # - conduction heat gain
        # - solar heat gain
        # - sensible ventilation heat gain
        # - sensible internal heat gain
        # - convective heat gain from the interior thermal mass (if it has been
        #   been defined)
        # - sensible zone load
        # - latent ventilation heat gain
        # - latent internal heat gain
        # - latent zone load
        n = len(Q_dot_cnd)
        loc_time, sol_time = zip(*[
            convert_to_clock_time(
                time_index=k,
                dt_hr=dt_hr,
                date=self.weather_data.date,
                L_loc=self.weather_data.location.L_loc,
                tz_loc=self.weather_data.location.timezone
            ) for k in range(n)
        ])
        d = {
            'solar time': sol_time,
            f'local time {self.weather_data.location.timezone}': loc_time,
            'Q_dot_cnd': Q_dot_cnd.to(unit).m,
            'Q_dot_sol': Q_dot_sol.to(unit).m,
            'Q_dot_sen_vent': Q_dot_sen_vent.m,
            'Q_dot_sen_ihg': Q_dot_ihg_sen.m
        }
        if Q_dot_tsn is not None:
            d['Q_dot_tsn'] = Q_dot_tsn.to(unit).m
        d.update({
            'Q_dot_sen_zone': Q_dot_sen_zone.to(unit).m,
            'Q_dot_lat_vent': Q_dot_lat_vent.m,
            'Q_dot_lat_ihg': Q_dot_ihg_lat.m,
            'Q_dot_lat_zone': Q_dot_lat_zone.m
        })
        df = pd.DataFrame(d)
        return df


class SpaceVentilation:
    """Calculation of the sensible and latent heat gain due to space
    ventilation according to standard EN 12831-1 (2017).
    """
    def __init__(self):
        self.n_min: float = 0.0
        self.V_dot_open: float = 0.0
        self.V_dot_ATD_d: float = 0.0
        self.V_dot_sup: float = 0.0
        self.V_dot_trf: float = 0.0
        self.V_dot_exh: float = 0.0
        self.V_dot_comb: float = 0.0
        self.T_sup: Callable[[float], Quantity] | None = None
        self.T_trf: Callable[[float], Quantity] | None = None
        self.T_db_ext: Callable[[float], Quantity] | None = None
        self.T_wb_ext: Callable[[float], Quantity] | None = None
        self.thz: ConditionedZone | None = None
        self.vez: VentilationZone | None = None
    
    @classmethod
    def create(
        cls,
        thz: ConditionedZone,
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_dot_open: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_ATD_d: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_sup: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_trf: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_exh: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_comb: Quantity = Q_(0.0, 'm ** 3 / hr'),
        T_sup: Callable[[float], Quantity] | None = None,
        T_trf: Callable[[float], Quantity] | None = None
    ) -> SpaceVentilation:
        """Configures the ventilation of the thermal zone according to standard
        EN 12831-1 (2017).

        Parameters
        ----------
        thz:
            The thermal zone to which the space ventilation is added.
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
        obj = cls()
        obj.thz = thz
        obj.vez = thz.vez
        obj.n_min = n_min.to('1 / hr').m
        obj.V_dot_open = V_dot_open.to('m ** 3 / hr').m
        obj.V_dot_ATD_d = V_dot_ATD_d.to('m ** 3 / hr').m
        obj.V_dot_sup = V_dot_sup.to('m ** 3 / hr').m
        obj.V_dot_trf = V_dot_trf.to('m ** 3 / hr').m
        obj.V_dot_exh = V_dot_exh.to('m ** 3 / hr').m
        obj.V_dot_comb = V_dot_comb.to('m ** 3 / hr').m
        obj.T_sup = T_sup or thz.weather_data.T_db
        obj.T_trf = T_trf
        obj.T_db_ext = thz.weather_data.T_db
        obj.T_wb_ext = thz.weather_data.T_wb
        return obj

    @property
    def V_dot_leak_ATD(self) -> float:
        """External airflow rate in m³/h into the space through leakages and ATDs.

        See EN 12831-1 (2017), eq. 19: The leakage airflow rate into the
        ventilation zone is divided among the spaces according to their envelope
        surface area, and the airflow rate through ATDs into the ventilation
        zone is divided among the spaces according to their design airflow rate
        through ATDs.
        """
        try:
            return (
                self.vez.V_dot_leak
                * (self.thz.envelope_area.to('m ** 2').m / self.vez.A_env)
                + self.vez.V_dot_ATD * (self.V_dot_ATD_d / self.vez.V_dot_ATD_d)
            )
        except ZeroDivisionError:
            try:
                return (
                    self.vez.V_dot_leak
                    * (self.thz.envelope_area.to('m ** 2').m / self.vez.A_env)
                )
            except ZeroDivisionError:
                return 0.0

    @property
    def V_dot_env(self) -> float:
        """External airflow rate in m³/h into the space through the envelope
        (see EN 12831-1 (2017), eq. 18).
        """
        try:
            V_dot_env = (
                (self.vez.V_dot_inf_add / self.vez.V_dot_env)
                * min(self.vez.V_dot_env, self.V_dot_leak_ATD * self.vez.f_dir)
            )
            V_dot_env += (
                (self.vez.V_dot_env - self.vez.V_dot_inf_add)
                / self.vez.V_dot_env * self.V_dot_leak_ATD
            )
        except ZeroDivisionError:
            return float('nan')
        else:
            return V_dot_env

    @property
    def V_dot_tech(self) -> float:
        """Technical airflow rate in m³/h into the space
        (see EN 12831-1 (2017), eq. 23).
        """
        return max(
            self.V_dot_sup + self.V_dot_trf,
            self.V_dot_exh + self.V_dot_comb
        )

    @property
    def V_dot_min(self) -> float:
        """The minimum required airflow rate in m³/h of the space that needs to
        be ensured to maintain an appropriate level of air hygiene
        (see EN 12831-1 eq. 33).
        """
        return self.n_min * self.thz.volume.to('m ** 3').m

    @property
    def V_dot_ext(self) -> float:
        """Volume flow rate in m³/h of outdoor air that enters the space."""
        V_dot_ext = max(
            self.V_dot_env + self.V_dot_open,
            self.V_dot_min - self.V_dot_tech
        )
        return V_dot_ext

    def Q_dot_sen(self, t_sol_sec: float) -> float:
        """Sensible infiltration/ventilation load in W of the space 
        (see EN 12831-1 eq. 17).
        """
        T_db_ext = self.T_db_ext(t_sol_sec).to('degC').m
        T_sup = self.T_sup(t_sol_sec).to('degC').m
        T_zone = self.thz.T_zone(t_sol_sec).to('degC').m
        VT_inf = self.V_dot_ext * (T_db_ext - T_zone)
        VT_sup = self.V_dot_sup * (T_sup - T_zone)
        if self.T_trf:
            T_trf = self.T_trf(t_sol_sec).to('degC').m
            VT_trf = self.V_dot_trf * (T_trf - T_zone)
        else:
            VT_trf = 0.0
        Q_dot_sen = 0.34 * (VT_inf + VT_sup + VT_trf)  # see EN 12831-1 B.2.8.
        return Q_dot_sen

    def Q_dot_lat(self, t_sol_sec: float) -> float:
        """Latent infiltration/ventilation load in W of the space (see
        Principles of Heating, Ventilation, and Air Conditioning in Buildings,
        Chapter 10.5).
        """
        h_o = HumidAir(
            Tdb=self.T_db_ext(t_sol_sec),
            Twb=self.T_wb_ext(t_sol_sec)
        ).h.to('J / kg').m
        h_x = HumidAir(
            Tdb=self.T_db_ext(t_sol_sec),
            RH=self.thz.RH_zone
        ).h.to('J / kg').m
        rho_o = HumidAir(
            Tdb=self.T_db_ext(t_sol_sec),
            Twb=self.T_wb_ext(t_sol_sec)
        ).rho.to('kg / m ** 3').m
        V_dot_inf = max(
            self.V_dot_env + self.V_dot_open,
            self.V_dot_min - self.V_dot_tech
        ) / 3600.0  # m³/s
        Q_dot_lat = rho_o * V_dot_inf * (h_o - h_x)
        return Q_dot_lat
