from __future__ import annotations

import pandas as pd
from hvac import Quantity

from .space import HeatedSpace, UnheatedSpace


Q_ = Quantity


class VentilationZone:

    def __init__(self):
        self.ID: str = ''
        self.T_ext_d: Quantity | None = None
        self._V_atd_d: Quantity | None = None
        self.dp_atd_d: Quantity = Q_(4.0, 'Pa')
        self.n_leak: float = 0.67
        self.q_env_50: Quantity = Q_(6, 'm ** 3 / (hr * m ** 2)')
        self.f_fac: float = 8
        self.f_qv: float = 0.05
        self.f_dir: float = 2
        self.f_iz: Quantity = Q_(0.5, 'frac')
        self._T_sup: Quantity | None = None
        self.eff_heat_recover: Quantity = Q_(90.0, 'pct')
        self.heated_spaces: dict[str, HeatedSpace] = {}
        self.unheated_spaces: dict[str, UnheatedSpace] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        T_ext_d: Quantity,
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
    ) -> 'VentilationZone':
        """Creates a ventilation zone.

        Parameters
        ----------
        ID: str
            Name that identifies the ventilation zone.
        T_ext_d:
            Design value of outdoor air temperature.
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
        self = cls()
        self.ID = ID
        self.T_ext_d = T_ext_d
        self._V_atd_d = V_atd_d
        self.dp_atd_d = dp_atd_d
        self.n_leak = n_leak
        self.q_env_50 = q_env_50
        self.f_fac = f_fac
        self.f_qv = f_qv
        self.f_dir = f_dir
        self.f_iz = f_iz
        self._T_sup = T_sup
        self.eff_heat_recover = eff_heat_recover
        return self

    def add_heated_space(
        self,
        ID: str,
        height: Quantity,
        area: Quantity,
        volume: Quantity | None = None,
        T_int_d: Quantity = Q_(20.0, 'degC'),
        grad_T_air: Quantity = Q_(1.0, 'K / m'),
        height_occ_zone: Quantity = Q_(1.0, 'm'),
        dT_surf: Quantity = Q_(0.0, 'K'),
        dT_rad: Quantity = Q_(0.0, 'K'),
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_atd_d: Quantity | None = None,
        V_exh: Quantity | None = None,
        V_comb: Quantity | None = None,
        V_sup: Quantity | None = None,
        V_trf: Quantity | None = None,
        V_open: Quantity | None = None,
        T_trf: Quantity | None = None,
        q_hu: Quantity | None = None
    ) -> HeatedSpace:
        """Adds a new heated space to the ventilation zone.

        Parameters
        ----------
        ID: str
            Name of heated space.
        height: Quantity
            Mean height of the heated space.
        area: Quantity
            Floor area of the heated space.
        volume: Quantity, default None
            Volume of the heated space. If not given, will be calculated as the
            product of `height` and `area`.
        T_int_d: Quantity, default 20 degC
            Internal design air temperature of the heated space.
        grad_T_air: Quantity, default 1 K/m
            Air temperature gradient of the heat emission system used in the room.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        height_occ_zone: Quantity, default 1 m
            Height of the occupied zone in the heated space.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        dT_surf: Quantity, default 0.0 K
            Correction term for the influence of the heat emission system on
            surface temperatures.
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        dT_rad: Quantity, default 0.0 K
            Difference between air and operative temperature (which can be
            approximated as the arithmetic average of air and mean radiant
            temperature).
            Only relevant if `height` >= 4 m (see EN 12831-1, B.2.6).
        n_min: Quantity, default 0.5 1/hr
            Minimum air change rate required for the heated space for reasons of
            air quality/hygiene and comfort (EN 12831-1, B.2.10 - Table B.7).
            Default value applies to permanent dwelling areas (living rooms,
            offices) and a ceiling height less than 3 m.
        V_atd_d: Quantity, default None
            Design air volume flow rate of the ATDs in the room
            (EN 12831-1, B.2.12).
            Only required if ATDs are used for ventilation.
        V_exh: Quantity, default None
            Exhaust ventilation air volume flow rate from the heated space.
        V_comb: Quantity, default None
            Air volume flow rate exhausted from the heated space that has not
            been included in the exhaust air volume flow of the ventilation
            system (typically, but not necessarily, combustion air if an open
            flue heater is present in the heated space).
        V_sup: Quantity, default None
            Supply air volume flow rate from the ventilation system into the
            heated space.
        V_trf: Quantity, default None
            Transfer air volume flow rate into the heated space from adjacent
            spaces.
        V_open: Quantity, default None
            External air volume flow rate into the heated space through large
            openings (EN 12831-1, Annex G).
        T_trf: Quantity, default None
            Temperature of the transfer air volume flow into the heated space
            from another space. In case the room height of the other space is
            less than 4 m, it is equal to the internal design temperature of
            the other space; otherwise, it is equal to mean air temperature of
            the other space (see EN 12831-1 §6.3.8.3).
        q_hu: Quantity, default None
            Additional heating-up power per unit floor area (EN 12831-1, Annex F).
        """
        hs = HeatedSpace.create(
            ID=ID,
            height=height,
            area=area,
            volume=volume,
            T_int_d=T_int_d,
            T_ext_d=self.T_ext_d,
            grad_T_air=grad_T_air,
            height_occ_zone=height_occ_zone,
            dT_surf=dT_surf,
            dT_rad=dT_rad,
            n_min=n_min,
            V_atd_d=V_atd_d,
            V_exh=V_exh,
            V_comb=V_comb,
            V_sup=V_sup,
            V_trf=V_trf,
            V_open=V_open,
            T_trf=T_trf,
            q_hu=q_hu,
            vz=self
        )
        self.heated_spaces[hs.ID] = hs
        return hs

    def add_unheated_space(
        self,
        ID: str,
        height: Quantity,
        area: Quantity,
        volume: Quantity | None = None,
        T_int_d: Quantity = Q_(20.0, 'degC'),
        grad_T_air: Quantity = Q_(1.0, 'K / m'),
        height_occ_zone: Quantity = Q_(1.0, 'm'),
        dT_surf: Quantity = Q_(0.0, 'K'),
        dT_rad: Quantity = Q_(0.0, 'K'),
        n_min: Quantity = Q_(0.5, '1 / hr'),
        V_atd_d: Quantity | None = None,
        V_exh: Quantity | None = None,
        V_comb: Quantity | None = None,
        V_sup: Quantity | None = None,
        V_trf: Quantity | None = None,
        V_open: Quantity | None = None,
        T_trf: Quantity | None = None,
    ) -> UnheatedSpace:
        """
        Adds an unheated space to ventilation zone.

        For an explanation of the parameters, see the docstring of method
        `add_heated_space`.
        """
        uhs = UnheatedSpace.create(
            ID=ID,
            height=height,
            area=area,
            volume=volume,
            T_int_d=T_int_d,
            T_ext_d=self.T_ext_d,
            grad_T_air=grad_T_air,
            height_occ_zone=height_occ_zone,
            dT_surf=dT_surf,
            dT_rad=dT_rad,
            n_min=n_min,
            V_atd_d=V_atd_d,
            V_exh=V_exh,
            V_comb=V_comb,
            V_sup=V_sup,
            V_trf=V_trf,
            V_open=V_open,
            T_trf=T_trf,
            vz=self
        )
        self.unheated_spaces[uhs.ID] = uhs
        return uhs

    @property
    def spaces(self) -> dict[str, HeatedSpace | UnheatedSpace]:
        """Returns the heated and unheated spaces in the ventilation zone."""
        spaces = dict(self.heated_spaces)
        spaces.update(self.unheated_spaces)
        return spaces

    @property
    def V_atd_d(self) -> Quantity:
        """Air volume flow rate through ATDs into the ventilation zone at design
        pressure difference.
        """
        if self._V_atd_d:
            V_atd_d = self._V_atd_d
        else:
            V_atd_d = sum(
                sp.V_atd_d
                for sp in self.spaces.values()
                if sp.V_atd_d is not None
            )
            if isinstance(V_atd_d, int):
                V_atd_d = Q_(0.0, 'm ** 3 / hr')
        return V_atd_d

    @property
    def V_atd_50(self) -> Quantity:
        """Air volume flow rate through ATDs into the ventilation zone at a
        pressure difference of 50 Pa (F.30).
        """
        V_atd_50 = self.V_atd_d * (Q_(50, 'Pa') / self.dp_atd_d) ** self.n_leak
        return V_atd_50

    @property
    def A_env(self) -> Quantity:
        """Envelope area of the zone."""
        A_env = sum(
            sp.A_env for sp in self.spaces.values()
        )
        return A_env or Q_(0.0, 'm ** 2')

    @property
    def V_exh(self) -> Quantity:
        """Total volume flow rate of ventilation exhaust air in the ventilation
        zone.
        """
        V_exh = sum(
            sp.V_exh
            for sp in self.spaces.values()
            if sp.V_exh is not None
        )
        return V_exh or Q_(0.0, 'm ** 3 / hr')

    @property
    def V_comb(self) -> Quantity:
        """Total volume flow rate of air exhausted from the ventilation zone
        that has not been included in the exhaust air volume flow of the
        ventilation system (typically, but not necessarily, combustion air if an
        open flue heater is present in the heated space).
        """
        V_comb = sum(
            sp.V_comb
            for sp in self.spaces.values()
            if sp.V_comb is not None
        )
        return V_comb or Q_(0.0, 'm ** 3 / hr')

    @property
    def V_sup(self) -> Quantity:
        """Total volume flow rate of supply air from the ventilation system into
        the heated space.
        """
        V_sup = sum(
            sp.V_sup
            for sp in self.spaces.values()
            if sp.V_sup is not None
        )
        return V_sup or Q_(0.0, 'm ** 3 / hr')

    @property
    def f_ez(self):
        """Adjustment factor taking into account the additional pressure
        difference due to unbalanced ventilation.
        """
        V_unbal = self.V_exh + self.V_comb - self.V_sup
        # external air volume flow rate into the ventilation zone due to
        # unbalanced ventilation
        V_inf = self.q_env_50 * self.A_env + self.V_atd_50
        # external air volume flow rate into the ventilation zone due to envelope
        # leaks and openings (ATDs) at 50 Pa
        f_ez = 1 / (1 + (self.f_fac / self.f_qv) * (V_unbal / V_inf) ** 2)
        return f_ez

    @property
    def V_inf_add(self) -> Quantity:
        """Air volume flow rate through additional infiltration into the
        ventilation zone.
        """
        V_inf_add = self.q_env_50 * self.A_env + self.V_atd_50
        V_inf_add *= self.f_qv * self.f_ez
        return V_inf_add

    @property
    def V_env(self) -> Quantity:
        """External air volume flow rate into the ventilation zone through the
        building envelope, being a combination of infiltration air due to
        unbalanced ventilation and due to envelope leaks and openings (ATDs)."""
        V_env = max(self.V_exh + self.V_comb - self.V_sup, 0) + self.V_inf_add
        return V_env

    @property
    def a_atd(self) -> Quantity:
        """Authority of ATDs (fraction of infiltration air that comes from
        ATDs)."""
        a_atd = self.V_atd_50 / (self.V_atd_50 + self.q_env_50 * self.A_env)
        return a_atd

    @property
    def V_leak(self) -> Quantity:
        """External air volume flow rate into the ventilation zone through
        leakages."""
        V_leak = (1 - self.a_atd) * self.V_env
        return V_leak

    @property
    def V_atd(self) -> Quantity:
        """External air volume flow into the ventilation zone through ATDs."""
        V_atd = self.a_atd * self.V_env
        return V_atd

    @property
    def T_exh(self) -> Quantity | None:
        """Average temperature of all the air exhausted in the spaces of the
        ventilation zone.
        """
        if V_exh := sum(sp.V_exh for sp in self.spaces.values()):
            T_exh = sum(
                sp.V_exh * sp.T_int_air.to('K')
                for sp in self.spaces.values()
            ) / V_exh
            return T_exh
        return None

    @property
    def T_sup(self) -> Quantity | None:
        """Temperature of the ventilation supply air to the spaces in the
        ventilation zone.
        """
        if self.T_exh is not None and self._T_sup is None:
            T_sup = self.T_ext_d.to('K') + self.eff_heat_recover * (self.T_exh - self.T_ext_d)
            return T_sup
        elif isinstance(self._T_sup, Quantity):
            return self._T_sup
        return None

    def get_ventilation_heat_loss(self) -> Quantity:
        """Returns the heat loss of the ventilation zone due to space
        ventilation and outdoor air infiltration.
        """
        rho_cp = Q_(0.34, 'W * hr / (m ** 3 * K)')
        Q_ven = rho_cp * sum(
            max(sp.V_leak_atd + sp.V_open, self.f_iz * sp.V_min - sp.V_tech) *
            (sp.T_int_air - sp.T_ext_d) +
            sp.V_sup * (sp.T_int_air - sp.T_sup) +
            sp.V_trf * (sp.T_int_air - sp.T_trf)
            for sp in self.spaces.values()
        )
        return Q_ven

    def get_transmission_heat_loss(self) -> Quantity:
        """Returns the heat loss of the ventilation zone due to heat conduction
        through the building elements.
        """
        Q_trm = sum(
            hs.get_transmission_heat_loss()
            for hs in self.heated_spaces.values()
        )
        return Q_trm

    def get_additional_heating_up_power(self) -> Quantity:
        """Returns the total heat rate needed to heat-up the spaces in the
        ventilation zone after a temperature set-back period.
        """
        Q_hu = sum(
            hs.get_additional_heating_up_power()
            for hs in self.heated_spaces.values()
        )
        return Q_hu

    def get_heat_load(self) -> Quantity:
        """Returns the total heating load of the ventilation zone."""
        Q_trm = self.get_transmission_heat_loss()
        Q_ven = self.get_ventilation_heat_loss()
        Q_hu = self.get_additional_heating_up_power()
        return Q_trm + Q_ven + Q_hu

    def get_summary(self, unit: str = 'kW', n_digits: int = 3) -> pd.DataFrame:
        """Returns a pandas DataFrame with an overview of the heated spaces in the
        ventilation zone together with their transmission heat loss, ventilation
        heat loss, additional heating up power, and total heating load.
        """
        col_1 = 'heated space'
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
        for hs in self.heated_spaces.values():
            d[col_1].append(hs.ID)
            d[col_2].append(round(hs.get_transmission_heat_loss().to(unit).m, n_digits))
            d[col_3].append(round(hs.get_ventilation_heat_loss().to(unit).m, n_digits))
            d[col_4].append(round(hs.get_additional_heating_up_power().to(unit).m, n_digits))
            d[col_5].append(round(hs.get_heat_load().to(unit).m, n_digits))
        df = pd.DataFrame(d)
        return df
