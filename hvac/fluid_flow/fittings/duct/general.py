import math
from hvac import Quantity
from hvac.fluids import Fluid
from hvac.fluid_flow.fittings.abstract_fitting import AbstractFitting
from hvac.fluid_flow.conduit import Duct


Q_ = Quantity

Air = Fluid('Air')
STANDARD_AIR = Air(T=Q_(20, 'degC'), P=Q_(101_325, 'Pa'))


class FlowCoefficient:

    @staticmethod
    def get_Av(volume_flow_rate: Quantity, pressure_drop: Quantity, rho: Quantity = STANDARD_AIR.rho) -> float:
        V = volume_flow_rate.to('m ** 3 / s').m
        dP = pressure_drop.to('Pa').m
        rho = rho.to('kg / m ** 3').m
        Av = V / math.sqrt(dP / rho)
        return Av


class ResistanceCoefficient:

    @staticmethod
    def from_Av(Av: float, diameter: Quantity) -> float:
        D = diameter.to('m').m
        zeta = math.pi ** 2 * D ** 4 / (8.0 * Av ** 2)
        return zeta


class DuctFitting(AbstractFitting):

    def __init__(
            self,
            duct: Duct,
            ID: str = '',
            zeta: float = float('nan'),  # resistance coefficient (linked to velocity pressure)
            Av: float = float('nan')     # flow coefficient (linked to volume flow rate)
    ):
        super().__init__(ID)
        self.duct = duct
        if math.isnan(zeta) and not math.isnan(Av):
            self._zeta = ResistanceCoefficient.from_Av(Av, self.duct.cross_section.equivalent_diameter)
        else:
            self._zeta = zeta

    @property
    def pressure_drop(self) -> Quantity:
        Pv = self.duct.velocity_pressure
        dP = self._zeta * Pv
        return dP

    @property
    def zeta(self) -> float:
        return self._zeta
