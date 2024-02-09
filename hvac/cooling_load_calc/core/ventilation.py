"""
Calculation of the minimum required ventilation air volume flow rates in
non-residential spaces.
"""
from abc import ABC, abstractmethod
import math
from hvac import Quantity


Q_ = Quantity


class VentilatedSpace(ABC):

    def __init__(
        self,
        ID: str,
        floor_area: Quantity | None = None
    ) -> None:
        self.ID = ID
        self.floor_area = floor_area
        pass

    @property
    @abstractmethod
    def V_dot_min(self) -> Quantity:
        """Returns the minimum required ventilation volume flow rate."""
        ...


class OccupiedVentilatedSpace(VentilatedSpace):
    """
    Ventilated space intended for human occupancy.
    """
    V_dot_pp: Quantity = Q_(22, 'm**3 / hr')
    # minimum required ventilation volume flow rate per person

    def __init__(
        self,
        ID: str,
        floor_area: Quantity,
        floor_area_pp_des: Quantity,
        num_persons_des: int | None = None
    ) -> None:
        """Instantiates an object of class `OccupiedVentilatedSpace`.

        Parameters
        ----------
        ID:
            Name of the space:
        floor_area:
            The total floor area of the space.
        floor_area_pp_des:
            Design value for the floor area per person in the space (according
            to an applicable standard).
        num_persons_des: optional
            Design value for the number of persons in the space (i.e., the
            design occupation of the space).

        Notes
        -----
        If the design value for the number of people in the space is smaller
        than the calculated minimum occupancy of this space, the actual number
        of people in the space is set to the minimum occupancy.
        """
        super().__init__(ID, floor_area)
        self.num_persons_des = num_persons_des
        self.floor_area_pp = floor_area_pp_des
        self.num_persons_min = math.ceil(self.floor_area / self.floor_area_pp)
        if self.num_persons_des is not None:
            self.num_persons = max(self.num_persons_des, self.num_persons_min)
        else:
            self.num_persons = self.num_persons_min

    @property
    def V_dot_min(self) -> Quantity:
        return self.num_persons * self.V_dot_pp


class NonOccupiedVentilatedSpace(VentilatedSpace):
    """
    Ventilated space not intended for human occupancy.
    """
    V_dot_pa: Quantity = Q_(1.3, 'm**3 / (hr * m**2)')

    @property
    def V_dot_min(self) -> Quantity:
        return self.floor_area * self.V_dot_pa


class Toilets(NonOccupiedVentilatedSpace):
    """
    Ventilated toilet space.
    """
    V_dot_pt: Quantity = Q_(25, 'm**3 / hr')   # per toilet
    V_dot_pa: Quantity = Q_(15, 'm**3 / (hr * m**2)')  # per unit floor area

    def __init__(
        self,
        ID: str,
        floor_area: Quantity | None,
        num_toilets: int | None = None
    ) -> None:
        """Creates a `Toilets` object.

        Parameters
        ----------
        ID:
            Name of the space:
        floor_area: optional
            The total floor area of the space.
        num_toilets: optional
            The number of toilets.

        Notes
        -----
        Either the floor area, or the number of toilets needs to be specified.
        If neither is specified, property `V_dot_min` cannot be determined and
        will return `None`.
        """
        super().__init__(ID, floor_area)
        self.num_toilets = num_toilets

    @property
    def V_dot_min(self) -> Quantity:
        if isinstance(self.floor_area, Quantity):
            return self.floor_area * self.V_dot_pa
        elif isinstance(self.num_toilets, int):
            return self.num_toilets * self.V_dot_pt
        return None


class Bathroom(NonOccupiedVentilatedSpace):
    """
    Ventilated bathroom or shower room.
    """
    V_dot_pa: Quantity = Q_(5, 'm**3 / (hr * m**2)')  # per unit floor area
    V_dot_min_br: Quantity = Q_(50, 'm**3 / hr')

    @property
    def V_dot_min(self) -> Quantity:
        return max(self.floor_area * self.V_dot_pa, self.V_dot_min_br)
