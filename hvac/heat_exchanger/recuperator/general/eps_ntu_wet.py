"""Module for designing and analyzing different types of heat exchangers based
on the effectiveness-NTU method for wet external surfaces (air coolers). The
following general types of heat exchangers are implemented:
1. counter-flow heat exchanger
2. parallel-flow heat exchanger
3. cross-flow heat exchanger
4. shell-and-tube heat exchanger.

The designing problem involves determining the required overall conductance UA
of a heat exchanger knowing the inlet conditions of the hot and the cold stream
and the required rate of heat transfer from the hot to the cold stream.

The analysis problem involves finding the resulting heat transfer from the hot
to the cold stream if the conductance UA of the heat exchanger and the inlet
conditions of the hot and cold stream are known.
"""
from abc import ABC, abstractmethod
from enum import IntEnum
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid, HumidAir, CP_HUMID_AIR, CoolPropError


Q_ = Quantity


class AbstractHeatExchanger(ABC):
    """Abstract base class from which concrete heat exchanger models are
    derived, implementing the attributes and methods common to all heat
    exchanger models.
    """
    def __init__(
        self,
        m_dot_r: Quantity,
        m_dot_a: Quantity,
        T_r_in: Quantity,
        T_r_out: Quantity | None,
        T_r_out_ini: Quantity | None,
        P_r: Quantity,
        refrigerant: Fluid,
        air_in: HumidAir,
        h_ext: Quantity | None = None,
        h_int: Quantity | None = None,
        eta_surf_wet: float | None = None,
        A_ext_to_A_int: float | None = None,
        A_ext: Quantity | None = None,
        Q: Quantity | None = None,
        **kwargs
    ) -> None:
        """Creates the common part of an object instantiated from any of the
        `AbstractHeatExchanger` derived classes.

        Parameters
        ----------
        m_dot_r:
            Mass flow rate of refrigerant.
        m_dot_a:
            Mass flow rate of humid air.
        T_r_in:
            Inlet temperature of refrigerant.
        T_r_out:
            Outlet temperature of refrigerant (if known, else set to None). In
            the case of a boiling refrigerant, set T_r_out equal to T_r_in.
        T_r_out_ini:
            Initial guess of the refrigerant outlet temperature (if it is not
            known, else set to None). In this case, also parameters h_ext,
            eta_surf_wet, A_ext_to_A_int, and h_int must be specified, so that
            the leaving refrigerant temperature can be calculated internally.
        P_r:
            Pressure of refrigerant.
        refrigerant:
            Type of refrigerant.
        air_in:
            State of entering humid air.
        h_ext: optional
            Convection heat transfer coefficient at wetted air-side surface.
        h_int: optional
            Convection heat transfer coefficient at internal, refrigerant-side
            surface.
        eta_surf_wet: optional
            Efficiency of wetted air-side finned surface.
        A_ext_to_A_int: optional
            Ratio of external to internal heat transfer surface area.
        A_ext: optional
            Area of external heat transfer surface.
        Q: optional
            Heat transfer rate between humid air and refrigerant.
        kwargs:
            Additional parameters.
            In the case of a cross-flow heat exchanger: parameter `type_`
            indicates the type of cross-flow arrangement, which must be one of
            the types defined by `CrossFlowHeatExchanger.Type`.
            In the case of a shell-and-tube heat exchanger: parameter `N`
            indicates the number of tube passes.
        """
        # properties of fluids
        self.m_dot_r = m_dot_r.to('kg / s')
        self.m_dot_a = m_dot_a.to('kg / s')
        self.T_r_in = T_r_in.to('K')
        self.T_r_out = T_r_out.to('K') if T_r_out is not None else None
        self.air_in = air_in
        self._P_r = P_r.to('Pa')
        self._Rfg = refrigerant
        # heat transfer characteristics
        if h_ext is not None:
            self._h_ext = h_ext.to('W / (m ** 2 * K)')
            self._h_a_wet = self._h_ext / CP_HUMID_AIR  # kg / (m ** 2 * s)
        else:
            self._h_ext = None
            self._h_a_wet = None
        self._h_int = h_int.to('W / (m ** 2 * K)') if h_int is not None else None
        self._eta_surf_wet = eta_surf_wet
        self._A_ext_to_A_int = A_ext_to_A_int
        self._A_ext = A_ext.to('m ** 2') if A_ext is not None else None
        self._Q = Q.to('W') if Q is not None else None
        # other parameters
        self._type = kwargs.get('type_')  # cross-flow heat exchanger
        self._N = kwargs.get('N')  # shell-tube heat exchanger
        # if the refrigerant outlet temperature is not known:
        if self.T_r_out is None and T_r_out_ini is not None:
            self.T_r_out = self._find_T_r_out(T_r_out_ini)

    def _xi(self, T_r_out: Quantity) -> Quantity:
        i_a_Tr_in = HumidAir(Tdb=self.T_r_in, RH=Q_(100, 'pct')).h
        i_a_Tr_out = HumidAir(Tdb=T_r_out, RH=Q_(100, 'pct')).h
        xi = (i_a_Tr_in - i_a_Tr_out) / (self.T_r_in - T_r_out)
        return xi

    def _C_r(self, xi: Quantity, T_r_out: Quantity) -> float:
        T_r_avg = (self.T_r_in + T_r_out) / 2
        try:
            cp_r = self._Rfg(T=T_r_avg, P=self._P_r).cp
        except CoolPropError:
            # refrigerant is a 2-phase mixture (boiling or condensing)
            return 0.0
        else:
            return (self.m_dot_a * xi) / (self.m_dot_r * cp_r)

    def _U_ext_wet(self, xi: Quantity) -> Quantity | None:
        if all([
            self._h_a_wet, self._eta_surf_wet,
            self._A_ext_to_A_int, self._h_int
        ]):
            R_ext = 1 / (self._h_a_wet * self._eta_surf_wet)  # (m ** 2 * s) / kg
            R_int = xi * self._A_ext_to_A_int / self._h_int  # (m ** 2 * s) / kg
            R_tot = R_ext + R_int  # (m ** 2 * s) / kg
            U_ext_wet = 1 / R_tot  # kg / (m ** 2 * s)
            return U_ext_wet.to('kg / (m ** 2 * s)')
        else:
            return None

    @property
    def Q_max(self) -> Quantity:
        """Maximum heat transfer rate between humid air and refrigerant
        that is theoretically possible.
        """
        i_a_Tr_in = HumidAir(Tdb=self.T_r_in, RH=Q_(100, 'pct')).h
        Q_max = self.m_dot_a * (self.air_in.h - i_a_Tr_in)
        return Q_max.to('W')

    @abstractmethod
    def __NTU__(self, C_r: float, eps: float) -> float:
        ...

    @abstractmethod
    def __eps__(self, C_r: float, NTU: float) -> float:
        ...

    def _find_T_r_out(self, T_r_out_ini: Quantity) -> Quantity:
        T_r_out = T_r_out_ini
        i_max = 10
        i = 0
        tol_T_r_out = Q_(0.1, 'K')
        while i <= i_max:
            xi = self._xi(T_r_out)
            UA_ext_wet = self._U_ext_wet(xi) * self._A_ext
            NTU = (UA_ext_wet / self.m_dot_a).m
            C_r = self._C_r(xi, T_r_out)
            eps = self.__eps__(C_r, NTU)
            T_r_avg = (self.T_r_in + T_r_out) / 2
            try:
                cp_r = self._Rfg(T=T_r_avg, P=self._P_r).cp
            except CoolPropError:
                # refrigerant is a 2-phase mixture (boiling or condensing)
                T_r_out_new = self.T_r_in
            else:
                T_r_out_new = self.T_r_in + (eps * self.Q_max) / (self.m_dot_r * cp_r)
            dev_T_r_out = abs(T_r_out_new.to('K') - T_r_out)
            if dev_T_r_out <= tol_T_r_out:
                self._T_r_out = T_r_out_new.to('K')
                return self._T_r_out
            T_r_out = T_r_out_new.to('K')
            i += 1
        else:
            raise ValueError(f'no acceptable solution found after {i_max} iterations')

    @property
    def xi(self) -> Quantity:
        dT_r = self.T_r_in - self.T_r_out
        if abs(dT_r.m) > 1.e-6:
            i_a_T_r_in = HumidAir(Tdb=self.T_r_in, RH=Q_(100, 'pct')).h
            i_a_T_r_out = HumidAir(Tdb=self.T_r_out, RH=Q_(100, 'pct')).h
            xi = (i_a_T_r_in - i_a_T_r_out) / dT_r
            return xi
        else:
            # If `self.T_r_in` and `self.T_r_out` are equal or very close to
            # each other (according to CoolProp <= 1.e-6), the refrigerant is a
            # 2-phase mixture (boiling or condensing).
            return Q_(0.0, 'J / (kg * K)')

    @property
    def C_r(self) -> float:
        T_r_avg = (self.T_r_in + self.T_r_out) / 2
        try:
            cp_r = self._Rfg(T=T_r_avg, P=self._P_r).cp
        except CoolPropError:
            # Refrigerant is a 2-phase mixture (boiling or condensing):
            return 0.0
        else:
            C_r = (self.m_dot_a * self.xi) / (self.m_dot_r * cp_r)
            return C_r.m

    @property
    def NTU(self) -> float:
        U_ext_wet = self._U_ext_wet(self.xi)
        if U_ext_wet is not None:
            # solving rating problem
            UA_ext_wet = U_ext_wet * self._A_ext  # kg / s
            NTU = UA_ext_wet / self.m_dot_a
            return NTU.m
        elif self.Q is not None:
            # solving sizing problem
            eps = (self.Q / self.Q_max).to('W / W').m
            return self.__NTU__(self.C_r, eps)
        return float('nan')

    @property
    def eps(self) -> float:
        eps = self.__eps__(self.C_r, self.NTU)
        return eps

    @property
    def Q(self) -> Quantity:
        if self._Q is None:
            Q = self.eps * self.Q_max
            return Q
        return self._Q

    @property
    def BF(self) -> float:
        n = -(self._h_ext * self._A_ext).to('W / K').m
        d = (self.m_dot_a * CP_HUMID_AIR).to('W / K').m
        bf = np.exp(n / d)  # bypass factor
        return bf

    @property
    def CF(self) -> float:
        return 1 - self.BF

    @property
    def air_out(self) -> HumidAir:
        i_a_out = self.air_in.h - self.Q / self.m_dot_a
        i_a_adp = self.air_in.h - (self.air_in.h - i_a_out) / self.CF
        air_adp = HumidAir(h=i_a_adp, RH=Q_(100, 'pct'))
        T_a_out = (
            air_adp.Tdb.to('K')
            + self.BF * (self.air_in.Tdb.to('K') - air_adp.Tdb.to('K'))
        )
        air_out = HumidAir(Tdb=T_a_out, h=i_a_out)
        return air_out

    @property
    def UA_ext_wet(self) -> Quantity:
        UA_ext_wet = self.NTU * self.m_dot_a
        return UA_ext_wet.to('kg / s')


class CounterFlowHeatExchanger(AbstractHeatExchanger):

    def __eps__(self, C_r: float, NTU: float) -> float:
        eps = float('nan')
        if C_r != 1.0:
            n = 1 - np.exp(-NTU * (1 - C_r))
            d = 1 - C_r * np.exp(-NTU * (1 - C_r))
            eps = n / d
        if C_r == 1.0:
            eps = NTU / (1 + NTU)
        return eps

    def __NTU__(self, C_r: float, eps: float) -> float:
        NTU = float('nan')
        if C_r != 1.0:
            NTU = (
                np.log((1 - eps * C_r) / (1 - eps))
                / (1 - C_r)
            )
        if C_r == 1.0:
            NTU = eps / (1 - eps)
        return NTU


class ParallelFlowHeatExchanger(AbstractHeatExchanger):

    def __eps__(self, C_r: float, NTU: float) -> float:
        n = 1 - np.exp(-NTU * (1 + C_r))
        d = 1 + C_r
        eps = n / d
        return eps

    def __NTU__(self, C_r: float, eps: float) -> float:
        n = np.log(1 - eps * (1 + C_r))
        d = 1 + C_r
        NTU = -n / d
        return NTU


class CrossFlowHeatExchanger(AbstractHeatExchanger):

    class Type(IntEnum):
        BOTH_FLUIDS_UNMIXED = 1  # both fluids cannot mix well in the direction perpendicular to the main flow
        BOTH_FLUIDS_MIXED = 2
        C_max_MIXED_C_min_UNMIXED = 3
        C_min_MIXED_C_max_UNMIXED = 4

    def __eps__(self, C_r: float, NTU: float) -> float:
        match self._type:
            case self.Type.BOTH_FLUIDS_UNMIXED:
                return self.__eps_01__(C_r, NTU)
            case self.Type.BOTH_FLUIDS_MIXED:
                return self.__eps_02__(C_r, NTU)
            case self.Type.C_max_MIXED_C_min_UNMIXED:
                return self.__eps_03__(C_r, NTU)
            case self.Type.C_min_MIXED_C_max_UNMIXED:
                return self.__eps_04__(C_r, NTU)
            case _:
                return float('nan')
    
    @staticmethod
    def __eps_01__(C_r: float, NTU: float) -> float:
        k = (
            (NTU ** 0.22) / C_r
            * (np.exp(-C_r * NTU ** 0.78) - 1)
        )
        eps = 1 - np.exp(k)
        return eps
    
    @staticmethod
    def __eps_02__(C_r: float, NTU: float) -> float:
        a = 1 / (1 - np.exp(-NTU))
        b = C_r / (1 - np.exp(-C_r * NTU))
        c = 1 / NTU
        eps = 1 / (a + b - c)
        return eps
    
    @staticmethod
    def __eps_03__(C_r: float, NTU: float) -> float:
        k = C_r * (np.exp(-NTU) - 1)
        eps = (1 - np.exp(k)) / C_r
        return eps
    
    @staticmethod
    def __eps_04__(C_r: float, NTU: float) -> float:
        k = (1 - np.exp(-C_r * NTU)) / C_r
        eps = 1 - np.exp(-k)
        return eps

    def __NTU__(self, C_r: float, eps: float) -> float:
        match self._type:
            case self.Type.BOTH_FLUIDS_UNMIXED:
                return self.__NTU_01__()
            case self.Type.BOTH_FLUIDS_MIXED:
                return self.__NTU_02__()
            case self.Type.C_max_MIXED_C_min_UNMIXED:
                return self.__NTU_03__(C_r, eps)
            case self.Type.C_min_MIXED_C_max_UNMIXED:
                return self.__NTU_04__(C_r, eps)
            case _:
                return float('nan')

    @staticmethod
    def __NTU_01__() -> float:
        return float('nan')

    @staticmethod
    def __NTU_02__() -> float:
        return float('nan')
    
    @staticmethod
    def __NTU_03__(C_r: float, eps: float) -> float:
        a = np.log(1 - eps * C_r) / C_r
        NTU = -np.log(1 + a)
        return NTU
    
    @staticmethod
    def __NTU_04__(C_r: float, eps: float) -> float:
        k = C_r * np.log(1 - eps) + 1
        NTU = -np.log(k) / C_r
        return NTU


class ShellAndTubeHeatExchanger(AbstractHeatExchanger):

    def __eps__(self, C_r: float, NTU: float) -> float:
        NTU_1 = NTU / self._N
        k = np.exp(-NTU_1 * np.sqrt(1 + C_r ** 2))
        n = 1 + k
        d = 1 - k
        z = 1 + C_r + np.sqrt(1 + C_r ** 2) * (n / d)
        eps_1 = 2 / z
        k = (1 - eps_1 * C_r) / (1 - eps_1)
        eps = ((k ** self._N) - 1) / ((k ** self._N) - C_r)
        return eps

    def __NTU__(self, C_r: float, eps: float) -> float:
        F = ((eps * C_r - 1) / (eps - 1)) ** (1 / self._N)
        eps_1 = (F - 1) / (F - C_r)
        E = (
            (2 - eps_1 * (1 + C_r))
            / (eps_1 * np.sqrt(1 + (C_r ** 2)))
        )
        n = np.log((E + 1) / (E - 1))
        d = np.sqrt(1 + (C_r ** 2))
        NTU_1 = n / d
        NTU = self._N * NTU_1
        return NTU
