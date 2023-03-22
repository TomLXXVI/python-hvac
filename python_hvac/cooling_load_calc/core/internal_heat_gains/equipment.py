from abc import ABC, abstractmethod
from enum import Enum
from hvac import Quantity
from ..schedule import OnOffSchedule

Q_ = Quantity


class Equipment(ABC):

    class Category(Enum):
        MACHINE = 'machine'
        HOODED_COOKING_APPLIANCE = 'hooded cooking appliance'
        OFFICE_APPLIANCE = 'office appliance'
        OFFICE_EQUIPMENT = 'office equipment'
        GENERIC_APPLIANCE = 'generic appliance'

    def __init__(self):
        self.ID: str = ''
        self.schedule: OnOffSchedule | None = None
        self.F_rad: Quantity = Q_(0, 'frac')
        self.Q_sen_rad: Quantity = Q_(0, 'W')
        self.Q_sen_conv: Quantity = Q_(0, 'W')
        self.Q_lat: Quantity = Q_(0, 'W')

    @classmethod
    @abstractmethod
    def create(cls, *args, **kwargs) -> 'Equipment':
        ...

    @abstractmethod
    def calculate_heat_gain(self) -> None:
        ...


class Machine(Equipment):

    class Configuration(Enum):
        MOTOR_AND_MACHINE = 'motor and machine'
        ONLY_MACHINE = 'only machine'
        ONLY_MOTOR = 'only motor'

    def __init__(self):
        super().__init__()
        self.P_motor: Quantity = Q_(0, 'W')
        self.eta_motor: Quantity = Q_(0, 'frac')
        self.configuration: Machine.Configuration = self.Configuration.ONLY_MACHINE

    @classmethod
    def create(
        cls,
        ID: str,
        P_motor: Quantity,
        eta_motor: Quantity,
        configuration: Configuration,
        schedule: OnOffSchedule,
        F_rad: Quantity = Q_(0.5, 'frac')
    ) -> 'Machine':
        machine = cls()
        machine.ID = ID
        machine.P_motor = P_motor
        machine.eta_motor = eta_motor.to('frac')
        machine.configuration = configuration
        machine.schedule = schedule
        machine.F_rad = F_rad.to('frac')
        return machine

    def calculate_heat_gain(self) -> None:
        if self.configuration == self.Configuration.ONLY_MACHINE:
            Q_sen = self.P_motor
        elif self.configuration == self.Configuration.ONLY_MOTOR:
            Q_sen = (1.0 - self.eta_motor) * (self.P_motor / self.eta_motor)
        else:
            Q_sen = self.P_motor / self.eta_motor
        self.Q_sen_rad = self.F_rad * Q_sen
        self.Q_sen_conv = Q_sen - self.Q_sen_rad


class HoodedCookingAppliance(Equipment):

    def __init__(self):
        super().__init__()
        self.P_rated: Quantity = Q_(0.0, 'W')
        self.F_U: Quantity = Q_(0.0, 'frac')

    @classmethod
    def create(
        cls,
        ID: str,
        P_rated: Quantity,
        F_U: Quantity,
        schedule: OnOffSchedule,
        F_rad: Quantity
    ) -> 'HoodedCookingAppliance':
        appliance = HoodedCookingAppliance()
        appliance.ID = ID
        appliance.P_rated = P_rated
        appliance.F_U = F_U
        appliance.F_rad = F_rad
        appliance.schedule = schedule
        return appliance

    def calculate_heat_gain(self):
        self.Q_sen_rad = self.F_rad * self.F_U * self.P_rated


class OfficeAppliance(Equipment):

    def __init__(self):
        super().__init__()
        self.P_cons: Quantity = Q_(0.0, 'W')

    @classmethod
    def create(
        cls,
        ID: str,
        P_peak: Quantity,
        schedule: OnOffSchedule,
        F_rad: Quantity = Q_(30.0, 'pct')
    ) -> 'OfficeAppliance':
        appliance = OfficeAppliance()
        appliance.ID = ID
        appliance.P_cons = P_peak
        appliance.schedule = schedule
        appliance.F_rad = F_rad
        return appliance

    def calculate_heat_gain(self):
        self.Q_sen_rad = self.F_rad * self.P_cons
        self.Q_sen_conv = self.P_cons - self.Q_sen_rad


class OfficeEquipment(Equipment):

    def __init__(self):
        super().__init__()
        self.heat_density: Quantity = Q_(0.0, 'W / m ** 2')
        self.A_floor: Quantity = Q_(0.0, 'm ** 2')

    @classmethod
    def create(
        cls,
        name: str,
        heat_density: Quantity,
        schedule: OnOffSchedule,
        F_rad: Quantity = Q_(30.0, 'pct'),
        floor_area: Quantity = Q_(1.0, 'm ** 2')
    ) -> 'OfficeEquipment':
        eqp = cls()
        eqp.name = name
        eqp.heat_density = heat_density
        eqp.A_floor = floor_area
        eqp.schedule = schedule
        eqp.F_rad = F_rad
        return eqp

    def calculate_heat_gain(self):
        Q_sen = self.heat_density * self.A_floor
        self.Q_sen_rad = self.F_rad * Q_sen
        self.Q_sen_conv = Q_sen - self.Q_sen_rad


class GenericAppliance(Equipment):

    def __init__(self):
        super().__init__()

    @classmethod
    def create(
        cls,
        name: str,
        Q_sen_rad: Quantity,
        Q_sen_con: Quantity,
        Q_lat: Quantity,
        schedule: OnOffSchedule
    ) -> 'GenericAppliance':
        eqp = cls()
        eqp.name = name
        eqp.Q_sen_rad = Q_sen_rad
        eqp.Q_sen_con = Q_sen_con
        eqp.Q_lat = Q_lat
        eqp.schedule = schedule
        return eqp

    def calculate_heat_gain(self):
        pass
