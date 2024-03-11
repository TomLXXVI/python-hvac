from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from hvac import Quantity

if TYPE_CHECKING:
    from hvac.cooling_load_calc.building.building_entity import BuildingEntity
    from hvac.cooling_load_calc.building.fixed_temperature_zone import FixedTemperatureZone


Q_ = Quantity


class VentilationZone:
    """Represents a single ventilation zone which is part of a building entity.
    A ventilation zone contains all the single conditioned or unconditioned
    thermal zones (spaces) that are served by the same ventilation system.
    It is possible that a building entity has more than one ventilation zone,
    but most often there will be one single ventilation zone per building
    entity.

    The class `VentilationZone` is used for the calculation of the sensible and
    latent heat gain due to space ventilation according to standard EN 12831-1
    (2017).
    """
    def __init__(self):
        self.ID: str = ''
        self.building_entity: BuildingEntity | None = None
        self.q_env_50: float = 0.0
        self.dP_ATD_d: float = 0.0
        self.v_leak: float = 0.0
        self.f_fac: float = 0.0
        self.f_V: float = 0.0
        self.f_dir: float = 0.0
        self.f_iz: float = 0.0
        self.thermal_zones: dict[str, FixedTemperatureZone] = {}
        self._is_solved: bool = False
        self._cooling_load_table: pd.DataFrame | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        q_env_50: Quantity = Q_(3.0, 'm ** 3 / (hr * m ** 2)'),
        dP_ATD_d: Quantity = Q_(4.0, 'Pa'),
        v_leak: float = 0.67,
        f_fac: float = 12,
        f_V: float = 0.05,
        f_dir: float = 2.0,
        f_iz: float = 0.5
    ) -> VentilationZone:
        """Creates a `VentilationZone` object.

        Parameters
        ----------
        ID:
            Identifies the ventilation zone.
        q_env_50:
            Air permeability of the building envelope at a pressure difference
            of 50 Pa between the interior and exterior with any ATDs closed or
            sealed (see EN 12831-1 (2017), B.2.10).
        dP_ATD_d:
            Design pressure difference of the ATDs in the zone (see EN 12831-1
            (2017), B.2.12).
        v_leak:
            Pressure exponent for air leakages (see EN 12831-1 (2017), B.2.13).
        f_fac:
            Adjustment factor for the number of wind exposed facades of the zone
            (see EN 12831-1 (2017), B.2.15). The default value applies to 1 wind
            exposed facade.
        f_V:
            Coefficient for the volume flow ratio of the zone (see EN 12831-1
            (2017), B.2.11, Table B.8). The default value applies to more than
            1 external facade, height of the zone above ground level between 0
            and 50 m, normal shielding, and a zone height between 5 and 10 m.
        f_dir:
            Factor for the orientation of the zone (see EN 12831-1 (2017),
            B.2.14). Default value according to B.2.14.
        f_iz:
            Ratio between the minimum air volume flow rates of single
            conditioned spaces and the air volume flow of the entire zone (see
            EN 12831-1 (2017), B.2.9, Table B.5). The default value applies
            to a zone with 2 or more spaces.

        Returns
        -------
        `VentilationZone` object
        """
        obj = cls()
        obj.ID = ID
        obj.q_env_50 = q_env_50.to('m ** 3 / (m ** 2 * hr)').m
        obj.dP_ATD_d = dP_ATD_d.to('Pa').m
        obj.v_leak = v_leak
        obj.f_fac = f_fac
        obj.f_V = f_V
        obj.f_dir = f_dir
        obj.f_iz = f_iz
        return obj

    def add_thermal_zone(self, *thz: FixedTemperatureZone) -> None:
        """Adds a thermal zone (or multiple zones) to the ventilation zone."""
        for thz_ in thz: thz_.ventilation_zone = self
        self.thermal_zones.update({thz_.ID: thz_ for thz_ in thz})

    def get_cooling_load_table(
        self,
        dt_hr: float = 1.0,
        num_cycles: int = 5,
        unit: str = 'W',
        do_recalculation: bool = False
    ) -> pd.DataFrame:
        """Returns the total cooling load of the ventilation zone.

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
        do_recalculation:
            If True, indicates that a full recalculation must be performed of
            the heat gains and cooling loads of all thermal zones in the
            ventilation zone, e.g. after modifications to thermal zones in the
            ventilation zone or after a modification of the ventilation zone
            itself.
            By default, if the cooling loads of the thermal zones were already
            calculated once, the original table is returned.

        Returns
        -------
        A Pandas DataFrame with the heat gains and the total cooling load
        of the ventilation zone for each time index k of the design day.
        """
        if not self._is_solved or do_recalculation:
            self._cooling_load_table = sum(
                thz.solve(dt_hr, num_cycles, unit)
                for thz in self.thermal_zones.values()
            )
            self._is_solved = True
        return self._cooling_load_table

    @property
    def floor_area(self) -> Quantity:
        """Returns the total floor area of the ventilation zone."""
        A_fl = sum(thz.floor_area for thz in self.thermal_zones.values())
        return A_fl
    
    @property
    def envelope_area(self) -> Quantity:
        """Returns the exterior envelope surface area of the ventilation zone."""
        A_env = sum(thz.envelope_area for thz in self.thermal_zones.values())
        return A_env
    
    @property
    def volume(self) -> Quantity:
        """Returns the volume of the ventilation zone."""
        V = sum(thz.volume for thz in self.thermal_zones.values())
        return V

    @property
    def V_dot_ATD_d(self) -> Quantity:
        """Design air volume flow rate in m³/h of the ATDs in the ventilation
        zone (EN 12831-1, Annex B.2.12). This is optional: only required when
        ATDs are present.
        """
        V_dot_ATD_d = sum(
            space.local_ventilation.V_dot_ATD_d.m
            for space in self.thermal_zones.values()
        )
        return Q_(V_dot_ATD_d, 'm**3 / hr')

    @property
    def V_dot_ATD_50(self) -> Quantity:
        """Airflow rate in m³/h into the ventilation zone through ATDs at a
        pressure difference of 50 Pa (EN 12831-1 eq. 30).
        """
        V_dot_ATD_50 = self.V_dot_ATD_d.m * (50.0 / self.dP_ATD_d) ** self.v_leak
        return Q_(V_dot_ATD_50, 'm**3 / hr')

    @property
    def f_ez(self) -> float:
        """Adjustment factor taking into account the additional pressure 
        difference due to unbalanced ventilation (EN 12831-1 eq. 29).
        """
        n = self.V_dot_exh.m + self.V_dot_comb.m - self.V_dot_sup.m
        d = self.q_env_50 * self.A_env + self.V_dot_ATD_50.m
        try:
            return 1 / (1 + (self.f_fac / self.f_V) * (n / d) ** 2)
        except ZeroDivisionError:
            return 1.0

    @property
    def A_env(self) -> float:
        """Envelope surface area in m² of the ventilation zone."""
        A_env = self.envelope_area.to('m ** 2').m
        return A_env

    @property
    def V_dot_exh(self) -> Quantity:
        """Exhaust airflow rate in m³/h from the ventilation zone."""
        V_dot_exh = sum(
            space.local_ventilation.V_dot_exh.m
            for space in self.thermal_zones.values()
        )
        return Q_(V_dot_exh, 'm**3 / hr')

    @property
    def V_dot_comb(self) -> Quantity:
        """Combustion airflow rate in m³/h into the ventilation zone."""
        V_dot_comb = sum(
            space.local_ventilation.V_dot_comb.m
            for space in self.thermal_zones.values()
        )
        return Q_(V_dot_comb, 'm**3 / hr')

    @property
    def V_dot_sup(self) -> Quantity:
        """Supply airflow rate in m³/h into the ventilation zone."""
        V_dot_sup = sum(
            space.local_ventilation.V_dot_sup.m
            for space in self.thermal_zones.values()
        )
        return Q_(V_dot_sup, 'm**3 / hr')

    @property
    def V_dot_inf_add(self) -> Quantity:
        """Airflow rate in m³/h due to additional infiltration into the
        ventilation zone, determined based on the air permeability (parameter
        `q_env_50`) and the airflow rate through ATDs (property `V_dot_ATD_50`).
        (EN 12831-1 eq. 28)
        """
        V_dot_inf_add = (
            (self.q_env_50 * self.A_env + self.V_dot_ATD_50.m)
            * self.f_V * self.f_ez
        )
        return Q_(V_dot_inf_add, 'm**3 / hr')

    @property
    def V_dot_env(self) -> Quantity:
        """External airflow rate in m³/h into the ventilation zone through the
        building envelope, determined taking into consideration the exhaust and
        supply airflow rates, the demand for combustion air and the airflow rate
        through additional infiltration. (EN 12831-1 eq. 24)
        """
        V_dot_env = max(
            self.V_dot_exh + self.V_dot_comb - self.V_dot_sup,
            Q_(0.0, 'm**3 / hr')
        )
        V_dot_env += self.V_dot_inf_add
        return V_dot_env.to('m**3 / hr')

    @property
    def a_ATD(self) -> float:
        """Authority of the ATDs in the zone (i.e., the ratio of airflow rate
        through ATDs only to the sum of the airflow rates through ATDs and
        through envelope air infiltration, both at a pressure difference of
        50 Pa). (EN 12831-1 eq. 22)
        """
        V_dot_leakage_50 = self.q_env_50 * self.A_env
        try:
            return self.V_dot_ATD_50.m / (self.V_dot_ATD_50.m + V_dot_leakage_50)
        except ZeroDivisionError:
            return 0.0

    @property
    def V_dot_leak(self) -> Quantity:
        """External airflow rate in m³/h into the ventilation zone through
        leakages (EN 12831-1 eq. 20).
        """
        V_dot_leak = (1 - self.a_ATD) * self.V_dot_env
        return V_dot_leak.to('m**3 / hr')

    @property
    def V_dot_ATD(self) -> Quantity:
        """External airflow rate in m³/h into the ventilation zone through ATDs
        (EN 12831-1 eq. 21).
        """
        V_dot_ATD = self.a_ATD * self.V_dot_env
        return V_dot_ATD.to('m**3 / hr')
