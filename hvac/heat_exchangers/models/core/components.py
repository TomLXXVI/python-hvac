from typing import Optional
from dataclasses import dataclass
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState, Fluid

Q_ = Quantity


@dataclass
class AirStream:
    inlet: Optional[HumidAir] = None
    outlet: Optional[HumidAir] = None
    m: Optional[Quantity] = None
    vfa: Optional[Quantity] = None


@dataclass
class CoolantStream:
    inlet: Optional[FluidState] = None
    outlet: Optional[FluidState] = None
    m: Optional[Quantity] = None
    m_tube: Optional[Quantity] = None
    coolant_type: Optional[Fluid] = None
    e: Optional[Quantity] = None


@dataclass
class HeatTransferParameters:
    he: Optional[Quantity] = None
    hi: Optional[Quantity] = None
    eta_surf: Optional[Quantity] = None
    xi: Optional[Quantity] = None
    Uo: Optional[Quantity] = None
    R: Optional[Quantity] = None


@dataclass
class Surface:
    Ti: Optional[Quantity] = None
    To: Optional[Quantity] = None

    @property
    def inlet(self) -> HumidAir:
        inlet = HumidAir(Tdb=self.Ti, RH=Q_(100, 'pct'))
        return inlet

    @property
    def outlet(self) -> HumidAir:
        outlet = HumidAir(Tdb=self.To, RH=Q_(100, 'pct'))
        return outlet


@dataclass
class HeatTransfer:
    delta_Qs: Optional[Quantity] = None
    delta_Ql: Optional[Quantity] = None
    delta_Q: Optional[Quantity] = None
    delta_mw: Optional[Quantity] = None
