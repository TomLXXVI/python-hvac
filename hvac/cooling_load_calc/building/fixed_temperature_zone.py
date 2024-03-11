from __future__ import annotations

from typing import Callable, TYPE_CHECKING

import pandas as pd

from hvac import Quantity
from hvac.cooling_load_calc.core import (
    WeatherData,
    ExteriorBuildingElement,
    InteriorBuildingElement,
    ThermalStorageNode,
    InternalHeatGain
)
from hvac.fluids import HumidAir

if TYPE_CHECKING:
    from .ventilation_zone import VentilationZone


Q_ = Quantity


class FixedTemperatureZone:
    """Represents a space or a group of spaces in a building of which the
    air temperature is perfectly held at its setpoint by the control action of
    the zone thermostat.
    This class is meant to be used in cooling load calculations, i.e., to find
    the required rate at which the cooling system must extract heat from the
    zone air to keep the zone air temperature fixed at its setpoint value.
    The sensible and latent heat gains and corresponding cooling loads are
    calculated at regular time moments during the design day (see method
    `solve`). By default, the time step is one hour, which means that the
    heat gains and cooling loads are calculated for each hour of the design
    day.
    """
    def __init__(self):
        self.ID: str = ''
        self.floor_area: Quantity | None = None
        self.height: Quantity | None = None
        self.weather_data: WeatherData | None = None
        self.T_zone: Callable[[float], Quantity] | None = None
        self.T_zone_des: Quantity | None = None
        self.RH_zone: Quantity | None = None
        self.ext_build_elems: dict[str, ExteriorBuildingElement] = {}
        self.int_build_elems: dict[str, InteriorBuildingElement] = {}
        self.tsn: ThermalStorageNode | None = None
        self.local_ventilation: LocalVentilation | None = None
        self.ventilation_zone: VentilationZone | None = None
        self.int_heat_gains: dict[str, InternalHeatGain] = {}
        self.thermal_storage_effect: dict[str, Quantity] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        weather_data: WeatherData,
        T_zone: Callable[[float], Quantity],
        floor_area: Quantity,
        height: Quantity,
        RH_zone: Quantity = Q_(50, 'pct'),
        ventilation_zone: VentilationZone | None = None,
        T_zone_des: Quantity = Q_(24, 'degC')
    ) -> FixedTemperatureZone:
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
        T_zone_des: optional
            The design value of the zone air temperature. The default value is
            24 °C. This value can be used to determine the thermal resistance of
            construction assemblies in the exterior building elements that
            surround the conditioned zone.
        """
        zone = cls()
        zone.ID = ID
        zone.weather_data = weather_data
        zone.T_zone = T_zone
        zone.T_zone_des = T_zone_des
        zone.floor_area = floor_area
        zone.height = height
        zone.RH_zone = RH_zone
        zone.ventilation_zone = ventilation_zone
        if ventilation_zone is not None:
            zone.ventilation_zone.add_thermal_zone(zone)
        return zone

    def add_ext_build_elem(self, *ebe: ExteriorBuildingElement) -> None:
        """Adds an exterior building element or a sequence of exterior building
        elements to the zone.
        """
        self.ext_build_elems.update({ebe_.ID: ebe_ for ebe_ in ebe})

    def add_int_build_elem(self, *ibe: InteriorBuildingElement) -> None:
        """Adds an interior building element or a sequence of interior building
        elements to the zone."""
        self.int_build_elems.update({ibe_.ID: ibe_ for ibe_ in ibe})

    def add_thermal_storage_node(
        self,
        C: Quantity = Q_(100, 'kJ / (K * m**2)'),
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
        self.local_ventilation = LocalVentilation.create(
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

    def add_internal_heat_gain(self, *ihg: InternalHeatGain) -> None:
        """Adds an internal heat gain or a sequence of internal heat gains to
        the zone (see cooling_load_calc.core.internal_heat_gains).
        """
        self.int_heat_gains.update({ihg_.ID: ihg_ for ihg_ in ihg})

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
        if d['Q_dot']:
            Q_dot = sum(d['Q_dot'])
            Q_dot_conv = sum(d['Q_dot_conv'])
            Q_dot_rad = sum(d['Q_dot_rad'])
        else:
            num_steps = int(round(24 / dt_hr))
            Q_dot = Q_dot_conv = Q_dot_rad = Q_([0.0] * num_steps, 'W')
        return Q_dot, Q_dot_conv, Q_dot_rad

    def _solar_heat_gain(
        self,
        dt_hr: float = 1.0
    ) -> tuple[Quantity, Quantity, Quantity]:
        """Returns the solar heat gain through all the windows of the zone,
        and also its radiative and convective components, for each time index
        of the design day.

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
        if d['Q_dot']:
            Q_dot = sum(d['Q_dot'])
            Q_dot_conv = sum(d['Q_dot_conv'])
            Q_dot_rad = sum(d['Q_dot_rad'])
        else:
            num_steps = int(round(24 / dt_hr))
            Q_dot = Q_dot_conv = Q_dot_rad = Q_([0.0] * num_steps, 'W')
        return Q_dot, Q_dot_conv, Q_dot_rad

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
        if self.local_ventilation is not None:
            Q_dot_sen_vent = Q_([
                self.local_ventilation.Q_dot_sen(k * dt_sec)
                for k in range(num_steps)
            ], 'W')
            Q_dot_lat_vent = Q_([
                self.local_ventilation.Q_dot_lat(k * dt_sec)
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
        if self.int_heat_gains:
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
        else:
            Q_dot_zero = Q_([0.0] * num_steps, 'W')
            Q_dot_sen = Q_dot_zero
            Q_dot_sen_cv = Q_dot_zero
            Q_dot_sen_rd = Q_dot_zero
            Q_dot_lat = Q_dot_zero
        return Q_dot_sen, Q_dot_sen_cv, Q_dot_sen_rd, Q_dot_lat

    def solve(
        self,
        dt_hr: float = 1.0,
        num_cycles: int = 5,
        unit: str = 'W'
    ) -> pd.DataFrame:
        """Calculates the sensible and latent heat gains in the zone and the
        resulting cooling load (i.e. the heat rate that the cooling system must
        extract from the zone air to keep the zone air temperature constant at
        its setpoint value).

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step for the calculations is one quarter of an hour.
            The default value is 1 hour.
        num_cycles:
            Number of diurnal calculation cycles before the results of the last
            diurnal cycle are returned.
            Each new cycle starts with the node temperatures from the last two
            time indexes k-2 and k-1 of the previous cycle as the initial values.
            Only the last cycle is kept.
        unit:
            The measuring unit in which heat flows need to be expressed. The
            default unit is Watts (W).

        Returns
        -------
        Pandas DataFrame with the heat gains (including both the convective and
        radiative component) and the total zone load at each time index of the
        design day.
        """
        # CONDUCTION HEAT GAIN
        # Get the total conductive heat gain in the zone, its convective and
        # radiative component:
        tup = self._conductive_heat_gain(dt_hr, num_cycles)
        Q_dot_cnd = tup[0]
        Q_dot_cnd_cv = tup[1]
        Q_dot_cnd_rd = tup[2]

        # SOLAR HEAT GAIN
        # Get the total solar heat gain in the zone, and its convective and
        # radiative components:
        tup = self._solar_heat_gain(dt_hr)
        Q_dot_sol = tup[0]
        Q_dot_sol_cv = tup[1]
        Q_dot_sol_rd = tup[2]

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
            # Keep the thermal node temperatures, the heat rate absorbed by
            # the thermal storage node, and the heat rate released to the
            # zone air in a dictionary.
            Q_dot_tsn = self.tsn.get_heat_flows()
            T_tsn = self.tsn.get_node_temperatures()
            self.thermal_storage_effect = {
                'T_tsn': T_tsn,
                'Q_dot_rad': Q_dot_rad,
                'Q_dot_tsn': Q_dot_tsn
            }
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
        d = {
            'Q_dot_cnd': Q_dot_cnd.to(unit).m,
            'Q_dot_sol': Q_dot_sol.to(unit).m,
            'Q_dot_sen_vent': Q_dot_sen_vent.to(unit).m,
            'Q_dot_sen_ihg': Q_dot_ihg_sen.to(unit).m
        }
        if Q_dot_sol is not None:
            d['Q_dot_sol'] = Q_dot_sol.to(unit).m

        if Q_dot_tsn is not None:
            d['Q_dot_tsn'] = Q_dot_tsn.to(unit).m
        d.update({
            'Q_dot_sen_zone': Q_dot_sen_zone.to(unit).m,
            'Q_dot_lat_vent': Q_dot_lat_vent.to(unit).m,
            'Q_dot_lat_ihg': Q_dot_ihg_lat.to(unit).m,
            'Q_dot_lat_zone': Q_dot_lat_zone.to(unit).m
        })
        df = pd.DataFrame(d)
        df['Q_dot_zone'] = df['Q_dot_sen_zone'] + df['Q_dot_lat_zone']
        df.loc['MAX'] = df.max()
        return df


class LocalVentilation:
    """Calculation of the sensible and latent heat gain due to space
    ventilation according to standard EN 12831-1 (2017).
    """
    def __init__(self):
        self.n_min: Quantity = Q_(0.0, '1 / hr')
        self.V_dot_open: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_ATD_d: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_sup: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_trf: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_exh: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_comb: Quantity = Q_(0.0, 'm**3 / hr')
        self.T_sup: Callable[[float], Quantity] | None = None
        self.T_trf: Callable[[float], Quantity] | None = None
        self.T_db_ext: Callable[[float], Quantity] | None = None
        self.T_wb_ext: Callable[[float], Quantity] | None = None
        self.thz: FixedTemperatureZone | None = None
        self.vez: VentilationZone | None = None
    
    @classmethod
    def create(
        cls,
        thz: FixedTemperatureZone,
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_dot_open: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_ATD_d: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_sup: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_trf: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_exh: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_comb: Quantity = Q_(0.0, 'm ** 3 / hr'),
        T_sup: Callable[[float], Quantity] | None = None,
        T_trf: Callable[[float], Quantity] | None = None
    ) -> LocalVentilation:
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
        V_dot_open:
            External air volume flow rate into the space through large openings
            (EN 12831-1, Annex G).
        V_dot_ATD_d:
            Design air volume flow rate of the ATDs in the room (EN 12831-1,
            B.2.12). Only relevant if ATDs (Air Terminal Devices) are used for
            ventilation (i.e., passive devices that allow air flow through a
            building element; it does not include the air out- or inlets of fan-
            assisted ventilation systems).
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
        obj.vez = thz.ventilation_zone
        obj.n_min = n_min.to('1 / hr')
        obj.V_dot_open = V_dot_open.to('m ** 3 / hr')
        obj.V_dot_ATD_d = V_dot_ATD_d.to('m ** 3 / hr')
        obj.V_dot_sup = V_dot_sup.to('m ** 3 / hr')
        obj.V_dot_trf = V_dot_trf.to('m ** 3 / hr')
        obj.V_dot_exh = V_dot_exh.to('m ** 3 / hr')
        obj.V_dot_comb = V_dot_comb.to('m ** 3 / hr')
        obj.T_sup = T_sup or thz.weather_data.T_db
        obj.T_trf = T_trf
        obj.T_db_ext = thz.weather_data.T_db
        obj.T_wb_ext = thz.weather_data.T_wb
        return obj

    @property
    def V_dot_leak_ATD(self) -> Quantity:
        """External airflow rate in m³/h into the space through leakages and ATDs.

        See EN 12831-1 (2017), eq. 19: The leakage airflow rate into the
        ventilation zone is divided among the spaces according to their envelope
        surface area, and the airflow rate through ATDs into the ventilation
        zone is divided among the spaces according to their design airflow rate
        through ATDs.
        """
        try:
            V_dot_leak_ATD = (
                self.vez.V_dot_leak
                * (self.thz.envelope_area.to('m ** 2').m / self.vez.A_env)
                + self.vez.V_dot_ATD * (self.V_dot_ATD_d.m / self.vez.V_dot_ATD_d.m)
            )
            return Q_(V_dot_leak_ATD, 'm**3 / hr')
        except ZeroDivisionError:
            try:
                V_dot_leak_ATD = (
                    self.vez.V_dot_leak
                    * (self.thz.envelope_area.to('m ** 2').m / self.vez.A_env)
                )
                return Q_(V_dot_leak_ATD, 'm**3 / hr')
            except ZeroDivisionError:
                return Q_(0.0, 'm**3 / hr')

    @property
    def V_dot_env(self) -> Quantity:
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
            return Q_(float('nan'), 'm**3 / hr')
        else:
            return V_dot_env.to('m**3 / hr')

    @property
    def V_dot_tech(self) -> Quantity:
        """Technical airflow rate in m³/h into the space
        (see EN 12831-1 (2017), eq. 23).
        """
        V_dot_tech = max(
            self.V_dot_sup + self.V_dot_trf,
            self.V_dot_exh + self.V_dot_comb
        )
        return V_dot_tech.to('m**3 / hr')

    @property
    def V_dot_min(self) -> Quantity:
        """The minimum required airflow rate in m³/h of the space that needs to
        be ensured to maintain an appropriate level of air hygiene
        (see EN 12831-1 eq. 33).
        """
        V_dot_min = self.n_min * self.thz.volume.to('m ** 3')
        return V_dot_min.to('m**3 / hr')

    @property
    def V_dot_ext(self) -> Quantity:
        """Volume flow rate in m³/h of outdoor air that enters the space."""
        V_dot_ext = max(
            self.V_dot_env.m + self.V_dot_open.m,
            self.V_dot_min.m - self.V_dot_tech.m
        )
        return Q_(V_dot_ext, 'm**3 / hr')

    def Q_dot_sen(self, t_sol_sec: float) -> float:
        """Sensible infiltration/ventilation load in W of the space 
        (see EN 12831-1 eq. 17).
        """
        T_db_ext = self.T_db_ext(t_sol_sec).to('degC').m
        T_sup = self.T_sup(t_sol_sec).to('degC').m
        T_zone = self.thz.T_zone(t_sol_sec).to('degC').m
        VT_inf = self.V_dot_ext.m * (T_db_ext - T_zone)
        VT_sup = self.V_dot_sup.m * (T_sup - T_zone)
        if self.T_trf:
            T_trf = self.T_trf(t_sol_sec).to('degC').m
            VT_trf = self.V_dot_trf.m * (T_trf - T_zone)
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
            self.V_dot_env.m + self.V_dot_open.m,
            self.V_dot_min.m - self.V_dot_tech.m
        ) / 3600.0  # m³/s
        Q_dot_lat = rho_o * V_dot_inf * (h_o - h_x)
        return Q_dot_lat
