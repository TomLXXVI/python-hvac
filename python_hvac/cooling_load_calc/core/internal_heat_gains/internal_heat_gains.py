from abc import ABC, abstractmethod
from hvac import Quantity
from ..schedule import OccupancySchedule
from .equipment import Equipment
from .lighting import Lighting


Q_ = Quantity


class InternalHeatGain(ABC):

    def __init__(self, ID: str):
        self.ID = ID

    @abstractmethod
    def Q_sen(self, t: float) -> dict[str, float]:
        """
        Get sensible fraction of internal heat gain at time moment t in seconds 
        from 00:00:00. 
        
        Returns a dictionary with the radiative component (key 'rad') and 
        convective component (key 'conv') of the sensible heat, expressed
        in Watts.
        """
        ...

    @abstractmethod
    def Q_lat(self, t: float) -> float:
        """
        Get latent internal heat gain at time moment t in seconds from 00:00:00 
        expressed in Watts.
        """
        ...


class EquipmentHeatGain(InternalHeatGain):

    def __init__(self, ID: str):
        super().__init__(ID)
        self.equipment: dict[str, Equipment] = {}

    def add_equipment(self, eqp: Equipment):
        """
        Add object of type `Equipment` (`Machine`, `HoodedCookingAppliance`, 
        `OfficeAppliance`, `OfficeEquipment` or `GenericAppliance).
        """
        self.equipment[eqp.ID] = eqp

    def remove_equipment(self, eqp_ID: str):
        """
        Remove equipment with ID `eqp_ID`.
        """
        self.equipment.pop(eqp_ID)

    def get_equipment(self, eqp_ID: str) -> Equipment:
        """
        Get equipment with ID `eqp_ID`.
        """
        return self.equipment.get(eqp_ID)

    def Q_sen(self, t: float) -> dict[str, float]:
        Q_sen_rad = 0.0
        Q_sen_conv = 0.0
        for eqp in self.equipment.values():
            if eqp.schedule(t):
                eqp.calculate_heat_gain()
                Q_sen_rad += eqp.Q_sen_rad.to('W').m
                Q_sen_conv += eqp.Q_sen_conv.to('W').m
        return {'rad': Q_sen_rad, 'conv': Q_sen_conv}

    def Q_lat(self, t: float) -> float:
        Q_lat = 0.0
        for eqp in self.equipment.values():
            if eqp.schedule(t):
                eqp.calculate_heat_gain()
                Q_lat += eqp.Q_lat.to('W').m
        return Q_lat


class LightingHeatGain(InternalHeatGain):

    def __init__(self, ID: str):
        super().__init__(ID)
        self.lighting: dict[str, Lighting] = {}

    def add_lighting(self, light: Lighting):
        """
        Add object of type `Lighting` (`LightingFixture` object or 
        `SpaceLighting` object).
        """
        self.lighting[light.ID] = light

    def remove_lighting(self, ID: str):
        """
        Remove `Lighting`-object with given ID.
        """
        self.lighting.pop(ID)

    def get_lighting(self, ID: str) -> Lighting:
        """
        Get `Lighting`-type object with given ID.
        """
        return self.lighting[ID]

    def Q_sen(self, t: float) -> dict[str, float]:
        Q_sen_rad = 0.0
        Q_sen_conv = 0.0
        for light in self.lighting.values():
            if light.schedule(t):
                light.calculate_heat_gain()
                Q_light = light.Q_light.to('W').m
                Q_sen_rad += light.F_rad * Q_light
                Q_sen_conv += Q_light - Q_sen_rad
        return {'rad': Q_sen_rad, 'conv': Q_sen_conv}

    def Q_lat(self, t: float) -> float:
        return 0.0


class PeopleHeatGain(InternalHeatGain):

    def __init__(self, ID: str):
        super().__init__(ID)
        self.Q_sen_person: Quantity = Q_(0.0, 'W')
        self.Q_lat_person: Quantity = Q_(0.0, 'W')
        self.F_rad: Quantity = Q_(0.0, 'frac')
        self.schedule: OccupancySchedule | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        Q_sen_person: Quantity,
        Q_lat_person: Quantity,
        F_rad: Quantity,
        schedule: OccupancySchedule
    ) -> 'PeopleHeatGain':
        """
        Create `PeopleHeatGain` object (see ASHRAE 2017, Ch. 18, table 1).

        Parameters
        ----------
        ID:
            Identifier for the internal heat gain
        Q_sen_person :
            Sensible heat release per person.
        Q_lat_person :
            Latent heat release per person.
        F_rad :
            Fraction of sensible heat release that is radiant.
        schedule :
            `OccupancySchedule` object that gives the number of people in the
            space depending on the time of the day.
        """
        phg = cls(ID)
        phg.Q_sen_person = Q_sen_person.to('W')
        phg.Q_lat_person = Q_lat_person.to('W')
        phg.F_rad = F_rad.to('frac')
        phg.schedule = schedule
        return phg

    def Q_sen(self, t: float) -> dict[str, float]:
        Q_sen = self.schedule(t) * self.Q_sen_person.m
        Q_sen_rad = self.F_rad * Q_sen
        Q_sen_conv = Q_sen - Q_sen_rad
        return {'rad': Q_sen_rad, 'conv': Q_sen_conv}

    def Q_lat(self, t: float) -> float:
        Q_lat = self.schedule(t) * self.Q_lat_person.m
        return Q_lat
