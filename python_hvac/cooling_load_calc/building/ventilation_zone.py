from typing import Optional, Dict, TYPE_CHECKING
import pandas as pd
from hvac import Quantity

if TYPE_CHECKING:
    from .building_entity import BuildingEntity
    from .space import Space

Q_ = Quantity


class VentilationZone:

    def __init__(self):
        self.ID: str = ''
        self.building_entity: Optional[BuildingEntity] = None
        self.spaces: Dict[str, Space] = {}
        self.heat_gains: Optional[pd.DataFrame] = None

        self.q_env_50: Optional[float] = 0.0
        self.dP_ATD_d: Optional[float] = 0.0
        self.v_leak: Optional[float] = 0.0
        self.f_fac: Optional[float] = 0.0
        self.f_V: Optional[float] = 0.0
        self.f_dir: Optional[float] = 0.0
        self.f_iz: Optional[float] = 0.0

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
    ) -> 'VentilationZone':
        """
        Create ventilation zone.

        Parameters
        ----------
        ID:
            Name of the ventilation zone.
        q_env_50:
            Air permeability of building envelope at a pressure difference of 50
            Pa between interior and exterior with any ATDs closed or sealed
            (NBN EN 12831-1, B.2.10).
        dP_ATD_d:
            Design pressure difference of the ATDs in the zone
            (NBN EN 12831-1, B.2.12).
        v_leak: float
            Pressure exponent for air leakages (NBN EN 12831-1, B.2.13).
        f_fac: float
            Adjustment factor for the number of wind exposed facades of the zone
            (NBN EN 12831-1, B.2.15). Default value applies to more than 1
            exposed facade.
        f_V: float
            Coefficient for the volume flow ratio of the zone (NBN EN 12831-1,
            B.2.11 - Table B.8). Default value applies to more than 1 exposed
            facade, height of the zone above ground level between 0 and 50
            m, normal shielding, and a zone height between 5 and 10 m.
        f_dir: float
            Factor for the orientation of the zone (NBN EN 12831-1, B.2.14).
            Default value according to B.2.14.
        f_iz: Quantity
            Ratio between the minimum air volume flow rates of single heated
            spaces and the air volume flow of the entire zone (NBN EN 12831-1,
            B.2.9 - Table B.5). Default value applies to a zone with 2 or more
            spaces.

        Returns
        -------
        Instance of `VentilationZone`.
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

    def add_space(self, space: 'Space'):
        self.spaces[space.ID] = space
        space.ventilation_zone = self

    def get_heat_gains(self, unit: str = 'W') -> pd.DataFrame:
        self.heat_gains = sum(
            space.get_heat_gains(unit)
            for space in self.spaces.values()
        )
        return self.heat_gains

    @property
    def floor_area(self) -> Quantity:
        A_floor = sum(
            space.floor_area
            for space in self.spaces.values()
        )
        return A_floor

    @property
    def envelope_area(self) -> Quantity:
        A_env = sum(
            space.envelope_area
            for space in self.spaces.values()
        )
        return A_env

    @property
    def volume(self) -> Quantity:
        V = sum(
            space.volume
            for space in self.spaces.values()
        )
        return V

    @property
    def V_ATD_d(self) -> float:
        """
        Design air volume flow of the ATDs in the ventilation zone.
        Optional: only required when ATDs are present.
        (EN 12831-1, Annex B.2.12).
        """
        V_ATD_d = sum(
            space.ventilation.V_ATD_d
            for space in self.spaces.values()
        )
        return V_ATD_d

    @property
    def V_ATD_50(self) -> float:
        """
        Airflow rate into the ventilation zone through ATDs at a pressure
        difference of 50 Pa.
        (EN 12831-1 eq. 30)
        """
        return self.V_ATD_d * (50.0 / self.dP_ATD_d) ** self.v_leak

    @property
    def f_ez(self) -> float:
        """
        Adjustment factor taking into account the additional pressure difference
        due to unbalanced ventilation.
        (EN 12831-1 eq. 29)
        """
        n = self.V_exh + self.V_comb - self.V_sup
        d = self.q_env_50 * self.A_env + self.V_ATD_50
        try:
            return 1 / (1 + (self.f_fac / self.f_V) * (n / d) ** 2)
        except ZeroDivisionError:
            return 1.0

    @property
    def A_env(self) -> float:
        """
        Envelope surface area of the ventilation zone in square meters (float).
        """
        A_env = self.envelope_area.to('m ** 2').m
        return A_env

    @property
    def V_exh(self) -> float:
        """
        Exhaust airflow rate from the ventilation zone.
        """
        V_exh = sum(
            space.ventilation.V_exh
            for space in self.spaces.values()
        )
        return V_exh

    @property
    def V_comb(self) -> float:
        """
        Combustion airflow rate into the ventilation zone.
        """
        V_comb = sum(
            space.ventilation.V_comb
            for space in self.spaces.values()
        )
        return V_comb

    @property
    def V_sup(self) -> float:
        """
        Supply airflow rate into the ventilation zone.
        """
        V_sup = sum(
            space.ventilation.V_sup
            for space in self.spaces.values()
        )
        return V_sup

    @property
    def V_inf_add(self) -> float:
        """
        Airflow rate through additional infiltration into the ventilation zone,
        determined based on the air permeability (parameter `q_env_50`) and the
        airflow rate through ATDs (property `V_ATD_50`).
        (EN 12831-1 eq. 28)
        """
        return (self.q_env_50 * self.A_env + self.V_ATD_50) * self.f_V * self.f_ez

    @property
    def V_env(self) -> float:
        """
        External airflow rate into the ventilation zone through the building
        envelope, determined taking into consideration the exhaust and supply
        airflow rates, the demand for combustion air and the airflow rate
        through additional infiltration.
        (EN 12831-1 eq. 24)
        """
        return max(self.V_exh + self.V_comb - self.V_sup, 0.0) + self.V_inf_add

    @property
    def a_ATD(self) -> float:
        """
        Authority of the ATDs in the zone (i.e. the ratio of airflow rate
        through ATDs at 50 Pa on the sum of airflow rate through ATDs at 50 Pa
        and leakage airflow rate at 50 Pa, which is equal to `q_env_50 * A_env`).
        (EN 12831-1 eq. 22)
        """
        try:
            return self.V_ATD_50 / (self.V_ATD_50 + self.q_env_50 * self.A_env)
        except ZeroDivisionError:
            return 0.0

    @property
    def V_leak(self) -> float:
        """
        External airflow rate into the ventilation zone through leakages.
        (EN 12831-1 eq. 20)
        """
        return (1 - self.a_ATD) * self.V_env

    @property
    def V_ATD(self) -> float:
        """
        External airflow rate into the ventilation zone through ATDs.
        (EN 12831-1 eq. 21)
        """
        return self.a_ATD * self.V_env
