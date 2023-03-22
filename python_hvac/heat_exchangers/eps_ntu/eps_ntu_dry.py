"""
ε-NTU method for dry heat exchangers.

User input data:
- global heat transfer coefficient `Uo` of heat exchanger
- external heat transfer surface `Ao` of heat exchanger
- entering hot stream `hot_in`
- entering cold stream `cold_in`
- mass flow rate of hot stream `mh`
- mass flow rate of cold stream `mc`

Output:
- effectiveness of heat exchanger `eps`
- actual heating/cooling capacity of heat exchanger `Q`
- temperature of leaving hot stream `Tho`
- temperature of leaving cold stream `Tco`

Use the appropriate class depending on the flow arrangement of the heat exchanger.
"""
from typing import Union
from abc import ABC, abstractmethod
import numpy as np
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState


class DryHeatExchanger(ABC):

    def __init__(
        self,
        Uo: Quantity,
        Ao: Quantity,
        hot_in: Union[HumidAir, FluidState],
        mh: Quantity,
        cold_in: Union[HumidAir, FluidState],
        mc: Quantity
    ) -> None:
        self.Uo = Uo
        self.Ao = Ao
        self.hot_in = hot_in
        self.mh = mh
        self.cold_in = cold_in
        self.mc = mc

        self._eps: Quantity = None
        self._Q: Quantity = None
        self._Tho: Quantity = None
        self._Tco: Quantity = None

    def _get_Ch(self) -> Quantity:
        return (self.hot_in.cp * self.mh).to('W / K')

    def _get_Cc(self) -> Quantity:
        return (self.cold_in.cp * self.mc).to('W / K')

    @staticmethod
    def _get_C_min(Ch: Quantity, Cc: Quantity) -> Quantity:
        return min(Cc, Ch)

    @staticmethod
    def _get_C_max(Ch: Quantity, Cc: Quantity) -> Quantity:
        return max(Cc, Ch)

    def _get_C_star(self, Ch: Quantity, Cc: Quantity) -> Quantity:
        C_min = self._get_C_min(Ch, Cc)
        C_max = self._get_C_max(Ch, Cc)
        C_star = C_min / C_max
        return C_star

    def _get_Q_max(self, C_min: Quantity) -> Quantity:
        if isinstance(self.hot_in, HumidAir):
            Thi = self.hot_in.Tdb
        else:
            Thi = self.hot_in.T
        if isinstance(self.cold_in, HumidAir):
            Tci = self.cold_in.Tdb
        else:
            Tci = self.cold_in.T
        Q_max = C_min * (Thi - Tci)
        return Q_max

    def _get_ntu(self, C_min: Quantity):
        ntu = self.Uo * self.Ao / C_min
        return ntu

    @abstractmethod
    def _get_eps(self, ntu: Quantity, C_star: Quantity) -> Quantity:
        ...

    @property
    def eps(self) -> Quantity:
        if self._eps is None:
            Ch = self._get_Ch()
            Cc = self._get_Cc()
            C_min = self._get_C_min(Ch, Cc)
            ntu = self._get_ntu(C_min)
            C_star = self._get_C_star(Ch, Cc)
            self._eps = self._get_eps(ntu, C_star)
        return self._eps

    @property
    def Q(self) -> Quantity:
        if self._Q is None:
            Ch = self._get_Ch()
            Cc = self._get_Cc()
            C_min = self._get_C_min(Ch, Cc)
            Q_max = self._get_Q_max(C_min)
            self._Q = self.eps * Q_max
        return self._Q

    @property
    def Tho(self) -> Quantity:
        if self._Tho is None:
            Ch = self._get_Ch()
            if isinstance(self.hot_in, HumidAir):
                Thi = self.hot_in.Tdb
            else:
                Thi = self.hot_in.T
            self._Tho = Thi - self.Q / Ch
        return self._Tho

    @property
    def Tco(self) -> Quantity:
        if self._Tco is None:
            Cc = self._get_Cc()
            if isinstance(self.cold_in, HumidAir):
                Tci = self.cold_in.Tdb
            else:
                Tci = self.cold_in.T
            self._Tco = Tci + self.Q / Cc
        return self._Tco


class DryCounterFlowHeatExchanger(DryHeatExchanger):

    def _get_eps(self, ntu: Quantity, C_star: Quantity) -> Quantity:
        k = np.exp(-ntu * (1 - C_star))
        eps = (1 - k) / (1 - C_star * k)
        return eps


class DryDXCoil(DryHeatExchanger):

    def _get_eps(self, ntu: Quantity, C_star: Quantity) -> Quantity:
        eps = 1 - np.exp(-ntu)
        return eps

    def _get_C_min(self, Ch: Quantity, Cc: Quantity) -> Quantity:
        if isinstance(self.hot_in, FluidState):
            return Cc
        elif isinstance(self.cold_in, FluidState):
            return Ch
        else:
            return min(Cc, Ch)
