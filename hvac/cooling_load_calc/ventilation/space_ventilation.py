from __future__ import annotations

import typing
from hvac import Quantity
from hvac.fluids import HumidAir

if typing.TYPE_CHECKING:
    from ..building.thermal_zone import ThermalZone
    from .ventilation_zone import VentilationZone

Q_ = Quantity


class ThermalZoneVentilation:
    """Calculation of the sensible and latent heat gain due to space ventilation
     according to standard EN 12831-1 (2017).
    """
    def __init__(self):
        self.n_min: Quantity = Q_(0.0, '1 / hr')
        self.V_dot_open: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_ATD_d: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_sup: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_trf_dict: dict[str, Quantity] | None = None
        self.V_dot_exh: Quantity = Q_(0.0, 'm**3 / hr')
        self.V_dot_comb: Quantity = Q_(0.0, 'm**3 / hr')
        self.T_sup: typing.Callable[[float], Quantity] | None = None
        self.thz: ThermalZone | None = None
        self.vez: VentilationZone | None = None

    @classmethod
    def create(
        cls,
        thz: ThermalZone,
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_dot_open: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_ATD_d: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_sup: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_trf_dict: dict[str, Quantity] | None = None,
        V_dot_exh: Quantity = Q_(0.0, 'm ** 3 / hr'),
        V_dot_comb: Quantity = Q_(0.0, 'm ** 3 / hr'),
        T_sup: typing.Callable[[float], Quantity] | None = None
    ) -> ThermalZoneVentilation:
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
            building element; it does not include the air out- or inlets of 
            fan-assisted ventilation systems).
        V_dot_sup:
            Supply air volume flow rate from the ventilation system into the
            space.
        V_dot_trf_dict:
            Dictionary with transfer air volume flow rates into the space from 
            adjacent spaces. The keys of the dictionary are the names of the
            adjacent spaces.
        V_dot_exh:
            Exhaust ventilation air volume flow rate from the space.
        V_dot_comb:
            Air volume flow rate exhausted from the space that has not been
            included in the exhaust air volume flow of the ventilation system
            (typically, but not necessarily, combustion air if an open flue
            heater is present in the heated space).
        T_sup:
            Temperature of the ventilation supply air to the zone. A callable
            that takes the solar time in seconds from midnight and returns the
            ventilation supply air temperature as a `Quantity` object.
        """
        obj = cls()
        obj.thz = thz
        obj.vez = thz.ventilation_zone
        obj.n_min = n_min.to('1 / hr')
        obj.V_dot_open = V_dot_open.to('m ** 3 / hr')
        obj.V_dot_ATD_d = V_dot_ATD_d.to('m ** 3 / hr')
        obj.V_dot_sup = V_dot_sup.to('m ** 3 / hr')
        obj.V_dot_trf_dict = V_dot_trf_dict
        obj.V_dot_exh = V_dot_exh.to('m ** 3 / hr')
        obj.V_dot_comb = V_dot_comb.to('m ** 3 / hr')
        obj.T_sup = T_sup
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
        if self.V_dot_trf_dict is not None:
            V_dot_trf = sum(
                V_dot_trf 
                for V_dot_trf in self.V_dot_trf_dict.values()
            )
        else:
            V_dot_trf = Q_(0.0, 'm**3 / hr')
        V_dot_tech = max(
            self.V_dot_sup + V_dot_trf,
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
    
    def Q_dot_lat(self, t_sol_sec: float) -> Quantity:
        """Latent infiltration/ventilation load in W of the space (see
        Principles of Heating, Ventilation, and Air Conditioning in Buildings,
        Chapter 10.5).
        """
        h_o = HumidAir(
            Tdb=self.thz.weather_data.T_db(t_sol_sec),
            Twb=self.thz.weather_data.T_wb(t_sol_sec)
        ).h.to('J / kg')
        h_x = HumidAir(
            Tdb=self.thz.weather_data.T_db(t_sol_sec),
            RH=self.thz.RH_zone
        ).h.to('J / kg')
        rho_o = HumidAir(
            Tdb=self.thz.weather_data.T_db(t_sol_sec),
            Twb=self.thz.weather_data.T_wb(t_sol_sec)
        ).rho.to('kg / m ** 3')
        V_dot_inf = max(
            self.V_dot_env + self.V_dot_open,
            self.V_dot_min - self.V_dot_tech
        ).to('m ** 3 / s')
        Q_dot_lat = rho_o * V_dot_inf * (h_o - h_x)
        return Q_dot_lat
