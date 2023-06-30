from typing import Union, Optional, TypeVar
import math
from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import fsolve
from hvac import Quantity, UNITS
from ..schedule import (
    PipeSchedule,
    DuctSchedule
)

u = UNITS


class CrossSection(ABC):

    def __init__(self):
        self.schedule: Optional[Union[PipeSchedule, DuctSchedule]] = None

    @classmethod
    @abstractmethod
    def create(cls, *args, **kwargs):
        """Create `CrossSection` instance."""
        ...

    @property
    @abstractmethod
    def area(self) -> Quantity:
        """Get cross-sectional area."""
        ...

    @property
    @abstractmethod
    def equivalent_diameter(self) -> Quantity:
        """Get equivalent diameter of cross-section."""
        ...

    @equivalent_diameter.setter
    @abstractmethod
    def equivalent_diameter(self, value: Quantity):
        """Set equivalent diameter of cross-section."""
        ...

    @property
    @abstractmethod
    def hydraulic_diameter(self) -> Quantity:
        """Get hydraulic diameter of cross-section."""
        ...


TCrossSection = TypeVar('TCrossSection', bound=CrossSection)


class Circular(CrossSection):

    def __init__(self):
        super().__init__()
        self.internal_diameter: Quantity = float('nan') * u.mm
        self.nominal_diameter: Quantity = float('nan') * u.mm

    @classmethod
    def create(
        cls,
        internal_diameter: Quantity = float('nan') * u.mm,
        nominal_diameter: Quantity = float('nan') * u.mm,
        schedule: Optional[Union[PipeSchedule, DuctSchedule]] = None
    ) -> 'Circular':
        """
        Create a circular cross-section.

        Either `internal_diameter` or `nominal diameter` must be specified;
        `internal_diameter` has precedence over `nominal_diameter` if both are
        given.

        For ducts: If `internal_diameter` is specified and `schedule` is also
        specified, the closest commercially available internal diameter is
        looked up for that duct schedule.

        For pipes: If `nominal_diameter` is specified and `schedule` is of type
        `PipeSchedule`, the corresponding `internal_diameter` is looked up for
        that pipe schedule. If `schedule` is not of type `PipeSchedule`, the
        value of argument `nominal_diameter` will be assigned to instance
        attribute `internal_diameter`.
        """
        obj = cls()
        obj.schedule = schedule
        if not math.isnan(internal_diameter.magnitude):
            if isinstance(obj.schedule, DuctSchedule):
                obj.internal_diameter = obj.schedule.get_closest_internal_size(internal_diameter)
            else:
                obj.internal_diameter = internal_diameter
        elif not math.isnan(nominal_diameter.magnitude):
            if isinstance(obj.schedule, PipeSchedule):
                obj.nominal_diameter = nominal_diameter
                obj.internal_diameter = obj.schedule.get_internal_diameter(nominal_diameter)
            else:
                obj.internal_diameter = nominal_diameter
        return obj

    @property
    def area(self) -> Quantity:
        return math.pi * self.internal_diameter ** 2 / 4.0

    @property
    def equivalent_diameter(self) -> Quantity:
        return self.internal_diameter

    @equivalent_diameter.setter
    def equivalent_diameter(self, value: Quantity):
        """Set diameter of circular cross-section.

        For pipes: If `self.schedule` is of type `PipeSchedule`, the closest
        commercially available nominal diameter and its corresponding internal
        diameter is looked up for that pipe schedule.

        For ducts: If `self.schedule` is of type `DuctSchedule`, the
        closest commercially available internal diameter is looked for that duct
        schedule.
        """
        if isinstance(self.schedule, PipeSchedule):
            self.nominal_diameter = self.schedule.get_closest_nominal_diameter(value)
            self.internal_diameter = self.schedule.get_internal_diameter(self.nominal_diameter)
        elif isinstance(self.schedule, DuctSchedule):
            self.internal_diameter = self.schedule.get_closest_internal_size(value)
        else:
            self.internal_diameter = value

    @property
    def hydraulic_diameter(self) -> Quantity:
        return self.internal_diameter


class Rectangular(CrossSection):

    def __init__(self):
        super().__init__()
        self.height: Quantity = float('nan') * u.mm
        self.width: Quantity = float('nan') * u.mm

    @classmethod
    def create(
        cls,
        height: Quantity,
        width: Quantity = float('nan') * u.mm,
        equivalent_diameter: Quantity = float('nan') * u.mm,
        schedule: Optional[DuctSchedule] = None
    ) -> 'Rectangular':
        """Creates a rectangular cross-section.

        To create a rectangular cross-section there are two options:
        1. Just specify height (and width) of the cross-section and leave
        `equivalent_diameter` and `schedule` to be `NaN` and `None`. Height
        must always be specified.
        2. Specify height and equivalent diameter of the cross-section.
        The width will be calculated from these specifications. If `schedule` is
        not None, the commercially available width closest to the calculated
        width will be used.
        """
        obj = cls()
        obj.schedule = schedule
        obj.height = height
        if not math.isnan(equivalent_diameter.magnitude):
            obj.width = obj._calculate_width(height, equivalent_diameter)
        else:
            obj.width = width
        return obj

    def _calculate_width(
        self,
        height: Quantity,
        equivalent_diameter: Quantity
    ) -> Quantity:
        h = height.to('m').magnitude
        d_eq = equivalent_diameter.to('m').magnitude

        def eq(unknowns: np.ndarray) -> np.ndarray:
            w = unknowns[0]
            out = 1.3 * (w * h) ** 0.625 / (w + h) ** 0.25 - d_eq
            return np.array([out])

        roots = fsolve(eq, np.array([0.01]))
        width = roots[0]
        width = width * u.m
        if isinstance(self.schedule, DuctSchedule):
            width = self.schedule.get_closest_internal_size(width)
        return width

    @property
    def area(self) -> Quantity:
        return self.height * self.width

    @property
    def equivalent_diameter(self) -> Quantity:
        return (
            1.3 * (self.width * self.height) ** 0.625
            / (self.width + self.height) ** 0.25
        )

    @equivalent_diameter.setter
    def equivalent_diameter(self, value: Quantity):
        """Setting the equivalent diameter of the rectangular cross-section will
        calculate the corresponding width of the cross-section. The height will
        remain unaffected. If attribute `schedule` is specified, the commercially
        available width closest to the calculated width will be used.
        """
        self.width = self._calculate_width(self.height, value)

    @property
    def hydraulic_diameter(self) -> Quantity:
        return 2.0 * self.width * self.height / (self.width + self.height)


class FlatOval(CrossSection):

    def __init__(self):
        super().__init__()
        self.height: Quantity = float('nan') * u.mm
        self.width: Quantity = float('nan') * u.mm

    @classmethod
    def create(
        cls,
        height: Quantity,
        width: Quantity = float('nan') * u.mm,
        equivalent_diameter: Quantity = float('nan') * u.mm,
        schedule: Optional[DuctSchedule] = None
    ) -> 'FlatOval':
        """Creates a flat-oval cross-section.

        To create a flat-oval cross-section there are two options:
        1. Just specify height (and width) of the cross-section and leave
        `equivalent_diameter` and `schedule` to be `NaN` and `None`. Height must
        always be specified.
        2. Specify height and equivalent diameter of the cross-section. The width
        will be calculated from these specifications. If `schedule` is specified,
        the commercially available width closest to the calculated width will
        be used.
        """
        obj = cls()
        obj.height = height
        obj.schedule = schedule
        if not math.isnan(equivalent_diameter.magnitude):
            obj.width = obj._calculate_width(height, equivalent_diameter)
        else:
            obj.width = width
        return obj

    def _calculate_width(
        self,
        height: Quantity,
        equivalent_diameter: Quantity
    ) -> Quantity:
        h = height.to('m').magnitude
        d_eq = equivalent_diameter.to('m').magnitude

        def eq(unknowns: np.ndarray) -> np.ndarray:
            w = unknowns[0]
            a = math.pi * h ** 2 / 4 + h * (w - h)
            p = math.pi * h + 2 * (w - h)
            out = 1.55 * a ** 0.625 / p ** 0.25 - d_eq
            return np.ndarray([out])

        roots = fsolve(eq, np.array([0.01]))
        width = roots[0]
        width = width * u.m
        if isinstance(self.schedule, DuctSchedule):
            width = self.schedule.get_closest_internal_size(width)
        return width

    @property
    def area(self) -> Quantity:
        return (
            math.pi * self.height ** 2 / 4.0
            + self.height * (self.width - self.height)
        )

    @property
    def equivalent_diameter(self) -> Quantity:
        perimeter = math.pi * self.height + 2 * (self.width - self.height)
        return 1.55 * self.area ** 0.625 / perimeter ** 0.25

    @equivalent_diameter.setter
    def equivalent_diameter(self, value: Quantity):
        """Setting the equivalent diameter of the flat-oval cross-section will
        calculate the width of the cross-section. The height will remain
        unaffected. If attribute `schedule` is also specified, the commercially
        available width closest to the calculated width will be used.
        """
        self.width = self._calculate_width(self.height, value)

    @property
    def hydraulic_diameter(self) -> Quantity:
        perimeter = math.pi * self.height + 2 * (self.width - self.height)
        return 4.0 * self.area / perimeter
