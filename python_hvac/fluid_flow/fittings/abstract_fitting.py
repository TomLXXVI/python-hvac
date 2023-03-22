from abc import ABC
from hvac import Quantity


class AbstractFitting(ABC):

    def __init__(self, ID: str = ''):
        self.ID = ID

    @property
    def zeta(self) -> float:
        return 0.0

    @property
    def zeta_c(self) -> float:
        return 0.0

    @property
    def zeta_b(self) -> float:
        return 0.0

    @property
    def zeta_s(self) -> float:
        return 0.0

    @property
    def zeta_large(self) -> float:
        return 0.0

    @property
    def zeta_small(self) -> float:
        return 0.0

    @property
    def pressure_drop(self) -> Quantity:
        return Quantity(0.0, 'Pa')

    @property
    def pressure_drop_main(self) -> Quantity:
        return Quantity(0.0, 'Pa')

    @property
    def pressure_drop_branch(self) -> Quantity:
        return Quantity(0.0, 'Pa')
