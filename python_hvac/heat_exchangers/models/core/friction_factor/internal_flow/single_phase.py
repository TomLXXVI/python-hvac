import numpy as np
from heat_exchanger import Quantity
from heat_exchanger.fluids import FluidState
from heat_exchanger.geometry import CoilGeometry
from heat_exchanger.heat_transfer.internal_flow import single_phase

Q_ = Quantity


class SinglePhaseFriction:

    def __init__(self, coolant: FluidState, geometry: CoilGeometry, mc: Quantity, Ao: Quantity):
        self.coolant = coolant
        self.Dh = geometry.Di
        Ai = (1 / geometry.Ao_to_Ai) * Ao
        self.L = Ai / (np.pi * geometry.Di)
        e = Q_(0, 'mm')  # smooth duct
        self.v = mc / (coolant.rho * np.pi * (geometry.Di ** 2) / 4)
        self._single_phase_flow = single_phase.SinglePhaseFlow(
            fluid=coolant,
            v=self.v,
            Dh=self.Dh,
            L=self.L,
            e=e,
            shape_of_duct=None
        )

    @property
    def f_dar(self) -> float:
        """Average Darcy (aka Moody) friction factor."""
        return self._single_phase_flow.f_avg

    @property
    def f_fan(self) -> float:
        """Average Fanning friction factor."""
        return self.f_dar / 4

    @property
    def delta_P(self) -> Quantity:
        """Pressure drop."""
        delta_P = self.f_dar * self.L / self.Dh * self.coolant.rho * (self.v ** 2) / 2
        return delta_P
