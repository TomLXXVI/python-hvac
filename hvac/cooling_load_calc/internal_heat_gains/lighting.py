from typing import Callable
from enum import Enum
from abc import ABC, abstractmethod
from hvac import Quantity


Q_ = Quantity


class Lighting(ABC):

    class Category(Enum):
        FIXTURE = 'lighting fixture'
        SPACE_LIGHTING = 'space lighting'

    def __init__(self):
        self.name: str = ''
        self.schedule: Callable[[float], float] | None = None
        self.F_space: float = 0.0
        self.F_rad: float = 0.0
        self.Q_dot_light: Quantity = Q_(0.0, 'W')
    
    @classmethod
    @abstractmethod
    def create(cls, *args, **kwargs) -> 'Lighting':
        ...

    @abstractmethod
    def calculate_heat_gain(self, t_sol_sec: float) -> None:
        ...


class LightingFixture(Lighting):
    """Represents a single lighting fixture or a group of lighting fixtures in
    a room. The light heat gain is calculated based on the lighting fixture's
    technical specifications.
    """
    def __init__(self):
        super().__init__()
        self.P_lamp: Quantity = Q_(0.0, 'W')
        self.F_allowance: float = 0.0

    @classmethod
    def create(
        cls,
        name: str,
        P_lamp: Quantity,
        F_allowance: float,
        F_use: float,
        F_rad: float,
        schedule: Callable[[float], float],
    ) -> 'LightingFixture':
        """Creates a `LightingFixture` object.
        (see ASHRAE Fundamentals 2017, Chapter 18, 2.2 Lighting).

        Parameters
        ----------
        name:
            Identifies the lighting fixture.
        P_lamp:
            The nominal wattage of the lamp(s) in the fixture without the power
            consumption of ballasts.
        F_allowance:
            The ratio of the lighting fixture's power consumption including
            lamps and ballast to the nominal power consumption of the lamps.
            For incandescent lamps, this factor is 1. For fluorescent and other
            lights, it accounts for power consumed by the ballast as well as the
            ballast's effect on lamp power consumption.
        F_use:
            The ratio of the lamp wattage in use to the total installed wattage
            for the conditions under which the cooling load estimate is being
            made. For commercial applications such as stores, the use factor is
            usually 1.0.
        F_rad:
            The radiative fraction is the radiative part of the lighting heat
            gain that goes to the room (see ASHRAE Fundamentals, Chapter 18,
            Table 3).
        schedule:
            Function with signature `f(t_sol_sec: float) -> float` that takes the
            solar time `t_sol_sec` in seconds from midnight (0 s) and returns
            a float between 0 and 1, where 0 stands for completely turned off
            and 1 for maximum lighting power.

        Notes
        -----
        A single `LightingFixture` object could also be used to represent a
        group of multiple lighting fixtures instead of only a single fixture.
        """
        fixture = cls()
        fixture.name = name
        fixture.P_lamp = P_lamp
        fixture.F_allowance = F_allowance
        fixture.F_space = F_use
        fixture.schedule = schedule
        fixture.F_rad = F_rad
        return fixture

    def calculate_heat_gain(self, t_sol_sec: float):
        div_fac = self.schedule(t_sol_sec)
        self.Q_dot_light = div_fac * self.F_space * self.P_lamp * self.F_allowance


class SpaceLighting(Lighting):
    """This class can be used to estimate the lighting gain on a per-unit area 
    basis (e.g., when final lighting plans are not available).
    """
    def __init__(self):
        super().__init__()
        self.power_density: Quantity = Q_(0.0, 'W / m ** 2')
        self.A_floor: Quantity = Q_(0.0, 'm ** 2')
        
    @classmethod
    def create(
        cls,
        name: str,
        power_density: Quantity,
        F_space: float,
        F_rad: float,
        schedule: Callable[[float], float],
        floor_area: Quantity = Q_(1.0, 'm ** 2')
    ) -> 'SpaceLighting':
        """Creates a `SpaceLighting` object.
        (see ASHRAE Fundamentals 2017, Chapter 18, Tables 2 and 3).

        Parameters
        ----------
        name:
            Identifies the space lighting.
        power_density:
            The maximum lighting power density, i.e. the lighting heat gain per
            square metre (see ASHRAE Fundamentals 2017, Chapter 18, Table 2).
        F_space:
            The space fraction, i.e., the fraction of lighting heat gain that
            goes to the room (the other fraction goes to the plenum) (see ASHRAE
            Fundamentals 2017, Chapter 18, Table 3 and Figure 3).
        F_rad:
            The radiative fraction is the radiative part of the lighting heat
            gain that goes to the room (see ASHRAE Fundamentals 2017, Chapter
            18, Table 3 and Figure 3).
        schedule:
            Function with signature `f(t_sol_sec: float) -> float` that takes the
            solar time `t_sol_sec` in seconds from midnight (0 s) and returns
            a float between 0 and 1 which indicates the diversity factor.
        floor_area
            The floor area of the zone.
        """
        lighting = cls()
        lighting.name = name
        lighting.power_density = power_density
        lighting.A_floor = floor_area
        lighting.F_space = F_space
        lighting.schedule = schedule
        lighting.F_rad = F_rad
        return lighting

    def calculate_heat_gain(self, t_sol_sec: float):
        div_fac = self.schedule(t_sol_sec)
        self.Q_dot_light = div_fac * self.F_space * self.power_density * self.A_floor
