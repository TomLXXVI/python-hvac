from __future__ import annotations
from typing import Callable
from abc import ABC, abstractmethod
from enum import Enum
from hvac import Quantity


Q_ = Quantity


class Equipment(ABC):

    class Category(Enum):
        MACHINE = 'machine'
        HOODED_COOKING_APPLIANCE = 'hooded_cooking_appliance'
        OFFICE_APPLIANCE = 'office_appliance'
        OFFICE_EQUIPMENT = 'office_equipment'
        GENERIC_APPLIANCE = 'generic_appliance'

    def __init__(self):
        self.name: str = ''
        self.schedule: Callable[[float], float] | None = None
        self.F_rad: float = 0.0
        self.Q_dot_sen_rd: Quantity = Q_(0, 'W')
        self.Q_dot_sen_cv: Quantity = Q_(0, 'W')
        self.Q_dot_lat: Quantity = Q_(0, 'W')

    @classmethod
    @abstractmethod
    def create(cls, *args, **kwargs) -> 'Equipment':
        ...

    @abstractmethod
    def calculate_heat_gain(self, t_sol_sec: float) -> None:
        """Calculates the convective sensible, the radiative sensible and the
        latent heat gain of the equipment or appliance, but doesn't return any
        results. The results are stored in attributes `Q_dot_sen_cv` (convective
        sensible heat gain), `Q_dot_sen_rd` (radiative sensible heat gain), and
        `Q_dot_lat` (latent heat gain).
        """
        ...


class Machine(Equipment):
    """Represents any machine driven by an electric motor (e.g., a fan, a pump,
    etc.).
    """
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
        name: str,
        P_motor: Quantity,
        eta_motor: Quantity,
        configuration: Configuration,
        schedule: Callable[[float], float],
        F_rad: float = 0.5
    ) -> 'Machine':
        """Creates a `Machine` object.

        Parameters
        ----------
        name:
            Identifies the machine.
        P_motor:
            Motor power nameplate rating.
        eta_motor:
            Motor efficiency. See e.g. ASHRAE Fundamentals 2017, Chapter 18,
            Tables 4A and 4B.
        configuration:
            See `Enum` subclass `Configuration` in class `Machine`:
            - `Configuration.MOTOR_AND_MACHINE`: the motor and the driven
               machine are both in the conditioned space.
            - `Configuration.ONLY_MACHINE`: only the driven machine is in the
               conditioned space.
            - `Configuration.ONLY_MOTOR`: only the motor is in the conditioned
               space (e.g. a fan in the conditioned space that exhausts air
               outside that space).
        schedule:
            Function with signature `f(t_sol_sec: float) -> float` that takes the
            solar time `t_sol_sec` in seconds since midnight (0 s) and returns
            a float between 0 and 1, where 0 stands for completely off and 1
            for running at full power.
        F_rad:
            The radiative fraction is the radiative part of the machine heat
            gain that goes to the room.
        """
        machine = cls()
        machine.name = name
        machine.P_motor = P_motor
        machine.eta_motor = eta_motor.to('frac')
        machine.configuration = configuration
        machine.schedule = schedule
        machine.F_rad = F_rad
        return machine

    def calculate_heat_gain(self, t_sol_sec: float) -> None:
        if self.configuration == self.Configuration.ONLY_MACHINE:
            Q_dot_sen = self.P_motor
        elif self.configuration == self.Configuration.ONLY_MOTOR:
            Q_dot_sen = (1.0 - self.eta_motor) * (self.P_motor / self.eta_motor)
        else:
            Q_dot_sen = self.P_motor / self.eta_motor
        div_fac = self.schedule(t_sol_sec)
        Q_dot_sen *= div_fac
        self.Q_dot_sen_rd = self.F_rad * Q_dot_sen
        self.Q_dot_sen_cv = Q_dot_sen - self.Q_dot_sen_rd


class HoodedCookingAppliance(Equipment):
    """Represents a cooking appliance installed under an effective hood; only
    radiant gain adds to the cooling load of the space.
    """

    def __init__(self):
        super().__init__()
        self.P_rated: Quantity = Q_(0.0, 'W')
        self.F_U: float = 0.0

    @classmethod
    def create(
        cls,
        name: str,
        P_rated: Quantity,
        F_U: float,
        F_rad: float,
        schedule: Callable[[float], float]
    ) -> 'HoodedCookingAppliance':
        """Creates a `HoodedCookingAppliance` object.
        (See ASHRAE Fundamentals 2017, Chapter 18, Tables 5A to 5E).

        Parameters
        ----------
        name:
            Identifies the appliance.
        P_rated:
            Nameplate or rated energy input rating.
        F_U:
            Usage factor applied to the nameplate rating that determines the
            average rate of appliance energy consumption.
        F_rad:
            The radiative fraction is the radiative part of the cooking
            appliance heat gain that goes to the space.
        schedule:
            Function with signature `f(t_sol_sec: float) -> float` that takes the
            solar time `t_sol_sec` in seconds from midnight (0 s) and returns
            a float between 0 and 1, where 0 stands for completely off and 1
            for running at full power.
        """
        appliance = HoodedCookingAppliance()
        appliance.name = name
        appliance.P_rated = P_rated
        appliance.F_U = F_U
        appliance.F_rad = F_rad
        appliance.schedule = schedule
        return appliance

    def calculate_heat_gain(self, t_sol_sec: float):
        div_fac = self.schedule(t_sol_sec)
        self.Q_dot_sen_cv = Q_(0.0, 'W')
        self.Q_dot_sen_rd = div_fac * self.F_rad * self.F_U * self.P_rated


class OfficeAppliance(Equipment):
    """Represents a single office appliance or a group of appliances in
    a room. The heat gain is calculated based on the technical specifications
    of the appliance.
    """

    def __init__(self):
        super().__init__()
        self.P_peak: Quantity = Q_(0.0, 'W')

    @classmethod
    def create(
        cls,
        name: str,
        P_peak: Quantity,
        F_rad: float,
        schedule: Callable[[float], float],
    ) -> 'OfficeAppliance':
        """Creates an `OfficeAppliance` object.

        Parameters
        ----------
        name:
            Identifies the office appliance.
        P_peak:
            Peak heat gain.
            - Computers: see ASHRAE Fundamentals 2017, Chapter 18, Tables 8A.
              Approximately 90% convective heat gain and 10% radiative heat
              gain.
            - Laptops and laptop docking station: see ASHRAE Fundamentals 2017,
              Chapter 18, Tables 8B. Approximately 75% convective heat gain and
              25% radiative heat gain.
            - Tablet PC: see ASHRAE Fundamentals 2017, Chapter 18, Tables 8C.
            - Monitors: see ASHRAE Fundamentals 2017, Chapter 18, Table 8D.
              Approximately 60% convective heat gain and 40% radiative heat
              gain.
            - Printers and copiers: see ASHRAE Fundamentals 2017, Chapter 18,
              Table 9. Approximately 70% convective heat gain and 30% radiative
              heat gain.
            - Miscellaneous office equipment: see ASHRAE Fundamentals 2017,
              Chapter 18, Table 10.
        F_rad:
            The radiative fraction is the radiative part of the office appliance
            heat gain that goes to the space.
        schedule:
            Function with signature `f(t_sol_sec: float) -> float` that takes the
            solar time `t_sol_sec` in seconds from midnight (0 s) and returns
            a float between 0 and 1, where 0 stands for completely off and 1
            for running at full power.
        """
        appliance = OfficeAppliance()
        appliance.name = name
        appliance.P_peak = P_peak
        appliance.schedule = schedule
        appliance.F_rad = F_rad
        return appliance

    def calculate_heat_gain(self, t_sol_sec: float):
        div_fac = self.schedule(t_sol_sec)
        P_peak = div_fac * self.P_peak
        self.Q_dot_sen_rd = self.F_rad * P_peak
        self.Q_dot_sen_cv = P_peak - self.Q_dot_sen_rd


class OfficeEquipment(Equipment):
    """This class can be used to estimate the heat gain of office equipment on a
    per-square-metre basis.
    """

    def __init__(self):
        super().__init__()
        self.heat_density: Quantity = Q_(0.0, 'W / m ** 2')
        self.A_floor: Quantity = Q_(0.0, 'm ** 2')

    @classmethod
    def create(
        cls,
        name: str,
        heat_gain_density: Quantity,
        floor_area: Quantity,
        schedule: Callable[[float], float],
        F_rad: float = 0.3,
    ) -> 'OfficeEquipment':
        """Creates an `OfficeEquipment` object.

        Parameters
        ----------
        name:
            Identifies the object.
        heat_gain_density:
            Heat gain per unit area, aka load factor. See ASHRAE Fundamentals
            2017, Chapter 18, Table 11. The medium heat gain density is likely
            to be appropriate for most standard offices.
        floor_area:
            The floor area of the space.
        schedule:
            Function with signature `f(t_sol_sec: float) -> float` that takes the
            solar time `t_sol_sec` in seconds from midnight (0 s) and returns
            a float between 0 and 1 which indicates the diversity factor.
        F_rad:
            The radiative fraction is the radiative part of the office equipment
            heat gain that goes to the space.
        """
        eqp = cls()
        eqp.name = name
        eqp.heat_density = heat_gain_density.to('W / m**2')
        eqp.A_floor = floor_area.to('m**2')
        eqp.schedule = schedule
        eqp.F_rad = F_rad
        return eqp

    def calculate_heat_gain(self, t_sol_sec: float):
        div_fac = self.schedule(t_sol_sec)
        Q_dot_sen = div_fac * self.heat_density * self.A_floor
        self.Q_dot_sen_rd = self.F_rad * Q_dot_sen
        self.Q_dot_sen_cv = Q_dot_sen - self.Q_dot_sen_rd


class GenericAppliance(Equipment):
    """Represents a generic appliance that doesn't belong to any of the
    other categories implemented above.
    """
    def __init__(self):
        super().__init__()

    @classmethod
    def create(
        cls,
        name: str,
        Q_dot_sen_rd: Quantity,
        Q_dot_sen_cv: Quantity,
        Q_dot_lat: Quantity,
        schedule: Callable[[float], float]
    ) -> 'GenericAppliance':
        """Creates a `GenericAppliance` object of which the heat gain components
        are already known.

        Parameters
        ----------
        name:
            Identifies the appliance.
        Q_dot_sen_rd:
            The radiative component of the sensible heat gain.
        Q_dot_sen_cv:
            The convective component of the sensible heat gain.
        Q_dot_lat:
            The latent heat gain from the appliance.
        schedule:
            Function with signature `f(t_sol_sec: float) -> float` that takes the
            solar time `t_sol_sec` in seconds from midnight (0 s) and returns
            a float between 0 and 1, where 0 stands for completely off and 1
            for running at full power.
        """
        eqp = cls()
        eqp.name = name
        eqp.Q_dot_sen_rd = Q_dot_sen_rd
        eqp.Q_dot_sen_cv = Q_dot_sen_cv
        eqp.Q_dot_lat = Q_dot_lat
        eqp.schedule = schedule
        return eqp

    def calculate_heat_gain(self, t_sol_sec: float):
        div_fac = self.schedule(t_sol_sec)
        self.Q_dot_sen_cv *= div_fac
        self.Q_dot_sen_rd *= div_fac
