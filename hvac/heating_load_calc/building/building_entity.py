from __future__ import annotations

import pandas as pd
from hvac import Quantity

from .ventilation_zone import VentilationZone


Q_ = Quantity


class BuildingEntity:

    def __init__(self):
        self.ID: str = ''
        self.T_ext_d: Quantity | None = None
        self.ventilation_zones: dict[str, VentilationZone] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        T_ext_d: Quantity
    ) -> BuildingEntity:
        """Creates a new building entity.

        Parameters
        ----------
        ID: str
            Name for the building entity.
        T_ext_d: Quantity
            Design value of the outdoor air temperature.
        """
        self = cls()
        self.ID = ID
        self.T_ext_d = T_ext_d
        return self

    def add_ventilation_zone(
        self,
        ID: str,
        V_atd_d: Quantity | None = None,
        dp_atd_d: Quantity = Q_(4.0, 'Pa'),
        n_leak: float = 0.67,
        q_env_50: Quantity = Q_(6.0, 'm ** 3 / (hr * m ** 2)'),
        f_fac: float = 8.0,
        f_qv: float = 0.05,
        f_dir: float = 2,
        f_iz: Quantity = Q_(0.5, 'frac'),
        T_sup: Quantity | None = None,
        eff_heat_recover: Quantity = Q_(90.0, 'pct')
    ) -> VentilationZone:
        """Adds a new ventilation zone to the building entity.

        Parameters
        ----------
        ID: str
            Name that identifies the ventilation zone.
        V_atd_d: Quantity, default None
            The design air volume flow rate of all the ATDs in the zone
            (EN 12831-1, B.2.12). Only required if ATDs are used for
            ventilation.
        dp_atd_d: Quantity, default 4 Pa
            Design pressure difference of the ATDs in the zone
            (EN 12831-1, B.2.12)
        n_leak: float, default 0.67
            Pressure exponent for air leakages (EN 12831-1, B.2.13)
        q_env_50: Quantity, default 6 m³/hr/m²
            Air permeability of building envelope at a pressure difference of
            50 Pa between interior and exterior with any ATDs closed or sealed
            (EN 12831-1, B.2.10). The default value of 6 m³/hr/m² corresponds
            with air tightness class III (i.e. an air tightness test has not been
            and will not be performed, while the requirement regarding air
            tightness is considered "mid-level")
        f_fac: float, default 8.0
            Adjustment factor for the number of wind exposed facades of the
            zone (EN 12831-1, B.2.15). The default value of 8 applies to more
            than 1 exposed facade. In case of 1 exposed facade `f_fac` is 12.
        f_qv: float, default 0.05
            Coefficient for the volume flow ratio of the zone (EN 12831-1,
            B.2.11 - Table B.8). Default value applies to more than 1 exposed
            facade, height of the zone above ground level between 0 and 50
            m, normal shielding, and a zone height between 5 and 10 m.
        f_dir: float, default 2.0
            Factor for the orientation of the zone (EN 12831-1, B.2.14).
            Default value according to B.2.14.
        f_iz: Quantity, default 0.5 frac
            Ratio between the minimum air volume flow rates of single heated
            spaces and the air volume flow of the entire zone (EN 12831-1,
            B.2.9 - Table B.5). Default value applies to a zone with 2 or more
            spaces.
        T_sup: Quantity, default None
            Temperature of the ventilation supply air after passing heat
            recovery or after passive preheating.
            Any temperature rise by active preheating, which requires power from
            a heat generator, shall not be taken into account (ventilation load
            is not a space load in that case).
            If `T_sup` is None, but `V_exh` of rooms in the ventilation zone is
            specified, then `T_sup` will be estimated according to
            EN 12831-1 §6.3.3.7
        eff_heat_recover: Quantity, default 90.0 pct
            Efficiency of the heat recovery of the ventilation system under
            design external conditions.
        """
        vz = VentilationZone.create(
            ID=ID,
            T_ext_d=self.T_ext_d,
            V_atd_d=V_atd_d,
            dp_atd_d=dp_atd_d,
            n_leak=n_leak,
            q_env_50=q_env_50,
            f_fac=f_fac,
            f_qv=f_qv,
            f_dir=f_dir,
            f_iz=f_iz,
            T_sup=T_sup,
            eff_heat_recover=eff_heat_recover
        )
        self.ventilation_zones[vz.ID] = vz
        return vz

    def get_transmission_heat_loss(self) -> Quantity:
        """Returns the heat loss of the building entity due to heat conduction
        through the building elements.
        """
        Q_trm = sum(
            hs.get_transmission_heat_loss()
            for vz in self.ventilation_zones.values()
            for hs in vz.spaces.values()
        )
        return Q_trm

    def get_ventilation_heat_loss(self) -> Quantity:
        """Returns the heat loss of the building entity due to space
        ventilation and outdoor air infiltration.
        """
        Q_ven = sum(
            vz.get_ventilation_heat_loss()
            for vz in self.ventilation_zones.values()
        )
        return Q_ven

    def get_additional_heating_up_power(self) -> Quantity:
        """Returns the total heat rate needed to heat-up the spaces in the
        building entity after a temperature set-back period.
        """
        Q_hu = sum(
            hs.get_additional_heating_up_power()
            for vz in self.ventilation_zones.values()
            for hs in vz.spaces.values()
        )
        return Q_hu

    def get_heat_load(self) -> Quantity:
        """Returns the total heating load of the building entity."""
        Q_trm = self.get_transmission_heat_loss()
        Q_ven = self.get_ventilation_heat_loss()
        Q_hu = self.get_additional_heating_up_power()
        return Q_trm + Q_ven + Q_hu

    def get_summary(self, unit: str = 'kW', n_digits: int = 3) -> pd.DataFrame:
        """Returns a Pandas DataFrame with a summary of the ventilation zones
        in the building entity.
        """
        col_1 = 'ventilation zone'
        col_2 = f'Q transmission [{unit}]'
        col_3 = f'Q ventilation [{unit}]'
        col_4 = f'Q heating-up [{unit}]'
        col_5 = f'Q total [{unit}]'
        d = {
            col_1: [],
            col_2: [],
            col_3: [],
            col_4: [],
            col_5: []
        }
        for vz in self.ventilation_zones.values():
            d[col_1].append(vz.ID)
            d[col_2].append(round(vz.get_transmission_heat_loss().to(unit).m, n_digits))
            d[col_3].append(round(vz.get_ventilation_heat_loss().to(unit).m, n_digits))
            d[col_4].append(round(vz.get_additional_heating_up_power().to(unit).m, n_digits))
            d[col_5].append(round(vz.get_heat_load().to(unit).m, n_digits))
        df = pd.DataFrame(d)
        return df
