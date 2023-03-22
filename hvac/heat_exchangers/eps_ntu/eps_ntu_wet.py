"""
ε-NTU method for wet air cooling heat exchangers.

User input data:
- average external energy transfer coefficient `he` (for both sensible and latent heat)
- average surface efficiency of external heat transfer surface `eta_surf`
- ratio of external to internal heat transfer surface `Ao_to_Ai`
- average internal heat transfer coefficient `hi` (convection)
- external heat transfer surface area `Ao`
- entering air condition `air_in`
- entering coolant condition `coolant_in`
- mass flow rate of air stream `ma`
- mass flow rate of coolant stream `mc`

Output:
- effectiveness of heat exchanger `eps`
- actual cooling capacity `Q`
- coolant leaving temperature `Tco`
- leaving air condition `air_out`
- moisture removal rate `mw`
- apparatus dew point `ADP`
- by-pass factor `bp`
- contact factor `cf`

Use the appropriate class depending on the flow arrangement of the heat exchanger.
"""
from typing import Optional
from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import fsolve
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState

Q_ = Quantity


class WetHeatExchanger(ABC):

    def __init__(
        self,
        he: Quantity,
        eta_surf: Quantity,
        Ao_to_Ai: float,
        hi: Quantity,
        Ao: Quantity,
        air_in: HumidAir,
        ma: Quantity,
        coolant_in: FluidState,
        mc: Quantity,
    ) -> None:
        self.he = he
        self.eta_surf = eta_surf
        self.Ao_to_Ai = Ao_to_Ai
        self.hi = hi
        self.Ao = Ao
        self.air_in = air_in
        self.ma = ma
        self.coolant_in = coolant_in
        self.mc = mc

        self._eps: Quantity = None
        self._Tco: Quantity = None
        self._hao: Quantity = None
        self._bf: Quantity = None
        self._cf: Quantity = None
        self._ADP: Quantity = None
        self._mw: Quantity = None
        self._air_out: Optional[HumidAir] = None
        self._Q: Quantity = None

    def _get_xi(self, Tco: Quantity) -> Quantity:
        Tci = self.coolant_in.T
        hci_sat = HumidAir(Tdb=Tci, RH=Q_(100, 'pct')).h
        hco_sat = HumidAir(Tdb=Tco, RH=Q_(100, 'pct')).h
        xi = (hci_sat - hco_sat) / (Tci - Tco)
        return xi

    def _get_Uow(self, xi: Quantity) -> Quantity:
        Uow = 1 / (1 / (self.he * self.eta_surf) + xi * self.Ao_to_Ai / self.hi)
        return Uow

    def _get_ntu(self, Uow: Quantity) -> Quantity:
        ntu = Uow * self.Ao / self.ma
        return ntu

    def _get_m_star(self, xi: Quantity) -> Quantity:
        cpc = self.coolant_in.cp
        m_star = self.ma * xi / (self.mc * cpc)
        return m_star

    def _get_Q_max(self) -> Quantity:
        hai = self.air_in.h
        Tci = self.coolant_in.T
        hci_sat = HumidAir(Tdb=Tci, RH=Q_(100, 'pct')).h
        Q_max = self.ma * (hai - hci_sat)
        return Q_max

    @abstractmethod
    def _get_eps(self, ntu: Quantity, m_star: Quantity) -> Quantity:
        ...

    @property
    @abstractmethod
    def Tco(self) -> Quantity:
        ...

    @property
    def eps(self) -> Quantity:
        if self._eps is None:
            xi = self._get_xi(self.Tco)
            Uow = self._get_Uow(xi)
            ntu = self._get_ntu(Uow)
            m_star = self._get_m_star(xi)
            eps = self._get_eps(ntu, m_star)
            self._eps = eps
        return self._eps

    def _get_hao(self) -> Quantity:
        if self._hao is None:
            hai = self.air_in.h
            Tci = self.coolant_in.T
            hci_sat = HumidAir(Tdb=Tci, RH=Q_(100, 'pct')).h
            hao = hai - self.eps * (hai - hci_sat)
            self._hao = hao
        return self._hao

    @property
    def bf(self) -> Quantity:
        if self._bf is None:
            bf = np.exp(-(self.he * self.Ao) / self.ma)
            self._bf = bf
        return self._bf

    @property
    def cf(self) -> Quantity:
        if self._cf is None:
            self._cf = 1 - self.bf
        return self._cf

    @property
    def ADP(self) -> HumidAir:
        if self._ADP is None:
            hai = self.air_in.h
            hao = self._get_hao()
            h_adp = hai - (hai - hao) / self.cf
            self._ADP = HumidAir(h=h_adp, RH=Q_(100, 'pct'))
        return self._ADP

    @property
    def air_out(self) -> HumidAir:
        if self._air_out is None:
            Tadp = self.ADP.Tdb
            Tai = self.air_in.Tdb
            Tao = Tadp + self.bf * (Tai - Tadp)
            hao = self._get_hao()
            self._air_out = HumidAir(Tdb=Tao, h=hao)
        return self._air_out

    @property
    def mw(self) -> Quantity:
        if self._mw is None:
            Wai = self.air_in.W
            Wao = self.air_out.W
            self._mw = self.ma * (Wai - Wao)
        return self._mw

    @property
    def Q(self) -> Quantity:
        if self._Q is None:
            Q_max = self._get_Q_max()
            self._Q = self.eps * Q_max
        return self._Q


class WetCounterFlowHeatExchanger(WetHeatExchanger):

    def _get_eps(self, ntu: Quantity, m_star: Quantity) -> Quantity:
        eps = (1 - np.exp(-ntu * (1 - m_star))) / (1 - m_star * np.exp(-ntu * (1 - m_star)))
        return eps

    @property
    def Tco(self) -> Quantity:
        if self._Tco is None:
            def eq(unknown: np.ndarray) -> np.ndarray:
                Tco = Q_(unknown[0], 'K')
                xi = self._get_xi(Tco)
                Uow = self._get_Uow(xi)
                ntu = self._get_ntu(Uow)
                m_star = self._get_m_star(xi)
                eps = self._get_eps(ntu, m_star)
                Tci = self.coolant_in.T
                Q_max = self._get_Q_max()
                cpc = self.coolant_in.cp
                out = Tco - (Tci + eps * Q_max / (self.mc * cpc))
                return np.array(out.to('K').m)

            T_co_ini = self.coolant_in.T.to('K') + Q_(5, 'K')
            roots = fsolve(eq, np.array([T_co_ini.m]))
            T_co = Q_(roots[0], 'K')
            self._Tco = T_co

        return self._Tco


class WetDXCoil(WetHeatExchanger):

    def _get_eps(self, ntu: Quantity, m_star: Quantity) -> Quantity:
        eps = 1 - np.exp(-ntu)
        return eps

    @property
    def Tco(self) -> Quantity:
        if self._Tco is None:
            def eq(unknown: np.ndarray) -> np.ndarray:
                Tco = Q_(unknown[0], 'K')
                xi = self._get_xi(Tco)
                Uow = self._get_Uow(xi)
                ntu = self._get_ntu(Uow)
                m_star = self._get_m_star(xi)
                eps = self._get_eps(ntu, m_star)
                Q_max = self._get_Q_max()
                hci = self.coolant_in.h
                hco = hci + eps * Q_max / self.mc
                Tco_new = self.coolant_in.fluid(h=hco, P=self.coolant_in.P).T
                out = Tco_new - Tco
                return np.array(out.to('K').m)

            T_co_ini = self.coolant_in.T.to('K') + Q_(5, 'K')
            roots = fsolve(eq, np.array([T_co_ini.m]))
            T_co = Q_(roots[0], 'K')
            self._Tco = T_co

        return self._Tco
