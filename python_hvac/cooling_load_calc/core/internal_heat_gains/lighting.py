from enum import Enum
from abc import ABC, abstractmethod
from hvac import Quantity
from ..schedule import OnOffSchedule

Q_ = Quantity


class Lighting(ABC):
    class Category(Enum):
        FIXTURE = 'lighting fixture'
        SPACE_LIGHTING = 'space lighting'

    def __init__(self):
        self.ID: str = ''
        self.schedule: OnOffSchedule | None = None
        self.F_space: Quantity = Q_(0.0, 'frac')
        self.F_rad: Quantity = Q_(0.0, 'frac')
        self.Q_light: Quantity = Q_(0.0, 'W')
    
    @classmethod
    @abstractmethod
    def create(cls, *args, **kwargs) -> 'Lighting':
        ...

    @abstractmethod
    def calculate_heat_gain(self):
        ...


class LightingFixture(Lighting):

    def __init__(self):
        super().__init__()
        self.P_lamp: Quantity = Q_(0.0, 'W')
        self.F_allowance: Quantity = Q_(0.0, 'frac')

    @classmethod
    def create(
        cls,
        ID: str,
        P_lamp: Quantity,
        F_allowance: Quantity,
        F_space: Quantity,
        schedule: OnOffSchedule,
        F_rad: Quantity
    ) -> 'LightingFixture':
        fixture = cls()
        fixture.ID = ID
        fixture.P_lamp = P_lamp
        fixture.F_allowance = F_allowance.to('frac')
        fixture.F_space = F_space.to('frac')
        fixture.schedule = schedule
        fixture.F_rad = F_rad.to('frac')
        return fixture

    def calculate_heat_gain(self):
        self.Q_light = self.F_space * self.P_lamp * self.F_allowance


class SpaceLighting(Lighting):

    def __init__(self):
        super().__init__()
        self.power_density: Quantity = Q_(0.0, 'W / m ** 2')
        self.A_floor: Quantity = Q_(0.0, 'm ** 2')
        
    @classmethod
    def create(
        cls,
        ID: str,
        power_density: Quantity,
        F_space: Quantity,
        schedule: OnOffSchedule,
        F_rad: Quantity,
        floor_area: Quantity = Q_(1.0, 'm ** 2')
    ) -> 'SpaceLighting':
        lighting = cls()
        lighting.ID = ID
        lighting.power_density = power_density
        lighting.A_floor = floor_area
        lighting.F_space = F_space.to('frac')
        lighting.schedule = schedule
        lighting.F_rad = F_rad.to('frac')
        return lighting

    def calculate_heat_gain(self):
        self.Q_light = self.F_space * self.power_density * self.A_floor
