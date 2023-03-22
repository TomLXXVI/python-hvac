from typing import Optional
import math
from hvac import Quantity
from hvac.fluids import Fluid, FluidState
from hvac.fluid_flow.schedule import pipe_schedule_40
from hvac.fluid_flow.conduit import Pipe
from ..abstract_fitting import AbstractFitting


Q_ = Quantity
Water = Fluid('Water')


class FlowCoefficient:
    water_15degC = Water(T=Q_(15.0, 'degC'), P=Q_(101.325, 'kPa'))
    rho_15degC = water_15degC.rho.to('kg / m ** 3').magnitude

    @classmethod
    def to_Kv(cls, Av: float) -> float:
        Kv = Av * 3.6e5 * math.sqrt(10) / math.sqrt(cls.rho_15degC)
        return Kv

    @classmethod
    def to_Av(cls, Kv: float) -> float:
        Av = Kv * math.sqrt(cls.rho_15degC) / (3.6e5 * math.sqrt(10))
        return Av

    @classmethod
    def get_Kv(cls, volume_flow_rate: Quantity, pressure_drop: Quantity, fluid: Optional[FluidState] = None) -> float:
        V = volume_flow_rate.to('m ** 3 / hr').magnitude
        dp = pressure_drop.to('bar').magnitude
        specific_gravity = 1.0
        if isinstance(fluid, FluidState):
            rho_fluid = fluid.rho.to('kg / m ** 3').magnitude
            specific_gravity = rho_fluid / cls.rho_15degC
        Kv = V / math.sqrt(dp / specific_gravity)
        return Kv


class ResistanceCoefficient:

    @staticmethod
    def from_Av(Av: float, diameter: Quantity) -> float:
        D = diameter.to('m').magnitude
        zeta = math.pi ** 2 * D ** 4 / (8.0 * Av ** 2)
        return zeta

    @staticmethod
    def from_Kv(Kv: float, diameter: Quantity) -> float:
        Av = FlowCoefficient.to_Av(Kv)
        return ResistanceCoefficient.from_Av(Av, diameter)

    @staticmethod
    def from_ELR(ELR: float, diameter: Quantity) -> float:
        D = diameter.to('mm').magnitude
        D_nom = pipe_schedule_40.get_closest_nominal_diameter(diameter)
        D_sch40 = pipe_schedule_40.get_internal_diameter(D_nom).to('mm').magnitude
        e = 0.046  # wall roughness of steel pipe schedule 40
        f = 0.25 / (math.log10(e / (3.7 * D_sch40))) ** 4
        return f * ELR * (D / D_sch40) ** 4


class PipeFitting(AbstractFitting):

    def __init__(
            self,
            pipe: Pipe,
            ID: str = '',
            Kv: float = float('nan'),
            zeta: float = float('nan'),
            zeta_inf: float = float('nan'),
            zeta_d: float = float('nan'),
            ELR: float = float('nan')
    ):
        super().__init__(ID)
        self.pipe = pipe
        self.Kv = Kv
        self._zeta = zeta
        self.zeta_inf = zeta_inf
        self.zeta_d = zeta_d
        self.ELR = ELR

    def _deltaP_from_Kv(self) -> Quantity:
        zeta = ResistanceCoefficient.from_Kv(self.Kv, self.pipe.cross_section.hydraulic_diameter)
        rho = self.pipe.fluid.rho
        v = self.pipe.velocity
        dp = zeta * rho * v ** 2 / 2
        return dp.to('Pa')

    def _deltaP_from_1K(self) -> Quantity:
        rho = self.pipe.fluid.rho
        v = self.pipe.velocity
        dp = self._zeta * rho * v ** 2 / 2
        return dp.to('Pa')

    def _deltaP_from_3K(self) -> Quantity:
        rho = self.pipe.fluid.rho
        mu = self.pipe.fluid.absolute_viscosity
        v = self.pipe.velocity
        pv = rho * v ** 2 / 2
        D = self.pipe.cross_section.hydraulic_diameter
        Re = v * D * rho / mu
        dp = (self._zeta / Re + self.zeta_inf * (1 + self.zeta_d / D ** 0.3)) * pv
        return dp.to('Pa')

    def _deltaP_from_ELR(self) -> Quantity:
        rho = self.pipe.fluid.rho
        v = self.pipe.velocity
        zeta = ResistanceCoefficient.from_ELR(self.ELR, self.pipe.cross_section.hydraulic_diameter)
        dp = zeta * rho * v ** 2 / 2
        return dp.to('Pa')

    @property
    def pressure_drop(self) -> Quantity:
        if not math.isnan(self.Kv):
            return self._deltaP_from_Kv()
        if not math.isnan(self.zeta_inf):
            return self._deltaP_from_3K()
        if not math.isnan(self.ELR):
            return self._deltaP_from_ELR()
        if not math.isnan(self._zeta):
            return self._deltaP_from_1K()

    @property
    def zeta(self) -> float:
        if not math.isnan(self.Kv):
            return ResistanceCoefficient.from_Kv(self.Kv, self.pipe.cross_section.hydraulic_diameter)
        if not math.isnan(self.zeta_inf):
            dp = self._deltaP_from_3K()
            rho = self.pipe.fluid.rho
            v = self.pipe.velocity
            pv = rho * v ** 2 / 2
            return (dp / pv).magnitude
        if not math.isnan(self.ELR):
            return ResistanceCoefficient.from_ELR(self.ELR, self.pipe.cross_section.hydraulic_diameter)
        if not math.isnan(self._zeta):
            return self._zeta


class BalancingValve(AbstractFitting):

    def __init__(
            self,
            pipe: Pipe,
            ID: str = '',
            pressure_drop_full_open: Quantity = Q_(float('nan'), 'Pa')
    ):
        super().__init__(ID)
        self.pipe = pipe
        self.fluid = self.pipe.fluid
        self.diameter = self.pipe.cross_section.hydraulic_diameter
        self.volume_flow_rate = self.pipe.volume_flow_rate
        self.pressure_drop_full_open = pressure_drop_full_open
        self._pressure_drop: Quantity = Q_(float('nan'), 'Pa')
        self.Kvs: float = float('nan')
        self.Kvr: float = float('nan')

    def calculate_preliminary_Kvs(self) -> float:
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        dp_100 = self.pressure_drop_full_open.to('Pa').magnitude
        Avs = V / math.sqrt(dp_100 / rho)
        Kvs = FlowCoefficient.to_Kv(Avs)
        return Kvs

    def set_Kvs(self, Kvs: float):
        self.Kvs = Kvs
        self._calculate_pressure_drop_full_open()

    def _calculate_pressure_drop_full_open(self):
        zeta = ResistanceCoefficient.from_Kv(self.Kvs, self.diameter)
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        D = self.diameter.to('m').magnitude
        A = math.pi * D ** 2 / 4
        v = V / A
        dp = zeta * rho * v ** 2 / 2
        self._pressure_drop = Q_(dp, 'Pa')

    def calculate_Kv_setting(self, extra_pressure_drop: Quantity) -> float:
        dp_extra = extra_pressure_drop.to('Pa').magnitude
        dp_open = self._pressure_drop.to('Pa').magnitude
        dp_final = dp_open + dp_extra
        self._pressure_drop = Q_(dp_final, 'Pa')
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        Avr = V / math.sqrt(dp_final / rho)
        self.Kvr = FlowCoefficient.to_Kv(Avr)
        return self.Kvr

    @property
    def zeta(self) -> float:
        if not math.isnan(self.Kvr):
            return ResistanceCoefficient.from_Kv(self.Kvr, self.pipe.cross_section.hydraulic_diameter)
        if not math.isnan(self.Kvs):
            return ResistanceCoefficient.from_Kv(self.Kvs, self.pipe.cross_section.hydraulic_diameter)
        return float('nan')

    @property
    def pressure_drop(self) -> Quantity:
        return self._pressure_drop


class ControlValve(AbstractFitting):

    def __init__(
            self,
            pipe: Pipe,
            ID: str = '',
            target_authority: float = 0.5,
            pressure_drop_crit_path: Quantity = Q_(float('nan'), 'Pa')
    ):
        super().__init__(ID)
        self.pipe = pipe
        self.fluid = self.pipe.fluid
        self.diameter = self.pipe.cross_section.hydraulic_diameter
        self.volume_flow_rate = self.pipe.volume_flow_rate
        self.target_authority = target_authority
        self.pressure_drop_crit_path = pressure_drop_crit_path
        self._pressure_drop: Quantity = Q_(float('nan'), 'Pa')
        self.Kvs: float = float('nan')

    def calculate_preliminary_Kvs(self) -> float:
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        dp_crit_path = self.pressure_drop_crit_path.to('Pa').magnitude  # still wo control valve
        a = self.target_authority
        dp_valve_full_open = a * dp_crit_path / (1 - a)
        Avs = V / math.sqrt(dp_valve_full_open / rho)
        Kvs = FlowCoefficient.to_Kv(Avs)
        return Kvs

    def set_Kvs(self, Kvs: float):
        self.Kvs = Kvs
        self._calculate_pressure_drop_full_open()

    def _calculate_pressure_drop_full_open(self):
        zeta = ResistanceCoefficient.from_Kv(self.Kvs, self.diameter)
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        D = self.diameter.to('m').magnitude
        A = math.pi * D ** 2 / 4
        v = V / A
        dp = zeta * rho * v ** 2 / 2
        self._pressure_drop = Q_(dp, 'Pa')

    def get_valve_authority(self, pressure_drop_crit_path: Quantity) -> float:
        # to be called when control valve has been added to the cross-over
        dp_valve_full_open = self._pressure_drop.to('Pa').magnitude
        dp_crit_path = pressure_drop_crit_path.to('Pa').magnitude
        return dp_valve_full_open / dp_crit_path

    @property
    def zeta(self) -> float:
        return ResistanceCoefficient.from_Kv(self.Kvs, self.pipe.cross_section.hydraulic_diameter)

    @property
    def pressure_drop(self) -> Quantity:
        return self._pressure_drop

    def set_valve_opening(self, percent_open: int) -> float:
        Kv = (percent_open / 100) * self.Kvs
        zeta = ResistanceCoefficient.from_Kv(Kv, self.pipe.cross_section.hydraulic_diameter)
        return zeta
