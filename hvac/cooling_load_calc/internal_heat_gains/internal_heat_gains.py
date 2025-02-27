from typing import Callable
from abc import ABC, abstractmethod
from hvac import Quantity
from .equipment import Equipment
from .lighting import Lighting


Q_ = Quantity


class InternalHeatGain(ABC):

    def __init__(self, name: str):
        self.name = name

    @abstractmethod
    def Q_dot(self, t_sol_sec: float) -> tuple[Quantity, ...]:
        """Returns a 3-tuple with the convective sensible, the radiative
        sensible, and the latent internal heat gain in Watts at solar time
        `t_sol_sec` in seconds from 00:00:00.
        """
        ...


class EquipmentHeatGain(InternalHeatGain):
    """Represents the internal heat gain from equipment in the space.

    An `EquipmentHeatGain` object contains one or more `Equipment` objects
    (see module equipment.py). The user must create these objects with the
    necessary input data so that the heat gain of this equipment can be
    calculated.
    """

    def __init__(self, name: str):
        super().__init__(name)
        self.equipment: dict[str, Equipment] = {}

    def add_equipment(self, eqp: Equipment):
        """Adds an `Equipment` object (any object of subclass `Machine`,
        `HoodedCookingAppliance`, `OfficeAppliance`, `OfficeEquipment` or
        `GenericAppliance) to this `EquipmentHeatGain` object.
        """
        self.equipment[eqp.name] = eqp

    def remove_equipment(self, eqp_ID: str):
        """Removes the equipment with name `eqp_ID` from this `EquipmentHeatGain`
        object.
        """
        self.equipment.pop(eqp_ID)

    def get_equipment(self, eqp_ID: str) -> Equipment:
        """Returns the equipment with name `eqp_ID`."""
        return self.equipment.get(eqp_ID)

    def Q_dot(self, t_sol_sec: float) -> tuple[Quantity, ...]:
        Q_dot_sen_rd = Q_(0.0, 'W')
        Q_dot_sen_cv = Q_(0.0, 'W')
        Q_dot_lat = Q_(0.0, 'W')
        for eqp in self.equipment.values():
            eqp.calculate_heat_gain(t_sol_sec)
            Q_dot_sen_rd += eqp.Q_dot_sen_rd.to('W')
            Q_dot_sen_cv += eqp.Q_dot_sen_cv.to('W')
            Q_dot_lat += eqp.Q_dot_lat.to('W')
        return Q_dot_sen_cv, Q_dot_sen_rd, Q_dot_lat


class LightingHeatGain(InternalHeatGain):
    """Represents the internal heat gain from space lighting.

    A `LightingHeatGain` object contains one or more `Lighting` objects
    (see module lighting.py). The user must create these objects with the
    necessary input data so that the heat gain of the lighting can be
    calculated.
    """

    def __init__(self, name: str):
        super().__init__(name)
        self.lighting: dict[str, Lighting] = {}

    def add_lighting(self, light: Lighting):
        """Adds a `Lighting` object (object of class `LightingFixture` or class
        `SpaceLighting`) to this `LightingHeatGain` object.
        """
        self.lighting[light.name] = light

    def remove_lighting(self, name: str):
        """Removes the `Lighting` object with given name from this
        `LightingHeatGain` object.
        """
        self.lighting.pop(name)

    def get_lighting(self, name: str) -> Lighting:
        """Returns the `Lighting` object with given name."""
        return self.lighting[name]

    def Q_dot(self, t_sol_sec: float) -> tuple[Quantity, ...]:
        Q_dot_sen_rd = Q_(0.0, 'W')
        Q_dot_sen_cv = Q_(0.0, 'W')
        for light in self.lighting.values():
            light.calculate_heat_gain(t_sol_sec)
            Q_dot_light = light.Q_dot_light.to('W')
            Q_dot_sen_rd += light.F_rad * Q_dot_light
            Q_dot_sen_cv += Q_dot_light - Q_dot_sen_rd
        return Q_dot_sen_cv, Q_dot_sen_rd, Q_(0.0, 'W')


class PeopleHeatGain(InternalHeatGain):
    """Represents the heat gain from people in the space."""

    def __init__(self, name: str):
        super().__init__(name)
        self.Q_dot_sen_person: Quantity = Q_(0.0, 'W')
        self.Q_dot_lat_person: Quantity = Q_(0.0, 'W')
        self.F_rad: float = 0.0
        self.schedule: Callable[[float], int] | None = None

    @classmethod
    def create(
        cls,
        name: str,
        Q_dot_sen_person: Quantity,
        Q_dot_lat_person: Quantity,
        F_rad: float,
        schedule: Callable[[float], int]
    ) -> 'PeopleHeatGain':
        """
        Creates a `PeopleHeatGain` object.
        (see ASHRAE Fundamentals 2017, Chapter 18, table 1).

        Parameters
        ----------
        name:
            Identifier for the internal heat gain
        Q_dot_sen_person :
            Sensible heat release per person.
        Q_dot_lat_person :
            Latent heat release per person.
        F_rad :
            Fraction of sensible heat release that is radiant.
        schedule :
            Function with signature `f(t_sol_sec: float) -> int` that takes the
            solar time in seconds from midnight (0 s) and returns the number
            of people in the thermal zone.
        """
        phg = cls(name)
        phg.Q_dot_sen_person = Q_dot_sen_person.to('W')
        phg.Q_dot_lat_person = Q_dot_lat_person.to('W')
        phg.F_rad = F_rad
        phg.schedule = schedule
        return phg

    def Q_dot(self, t_sol_sec: float) -> tuple[Quantity, ...]:
        Q_dot_sen = self.schedule(t_sol_sec) * self.Q_dot_sen_person
        Q_dot_sen_rd = self.F_rad * Q_dot_sen
        Q_dot_sen_cv = Q_dot_sen - Q_dot_sen_rd
        Q_dot_lat = self.schedule(t_sol_sec) * self.Q_dot_lat_person
        return Q_dot_sen_cv, Q_dot_sen_rd, Q_dot_lat
