"""Module for designing and analyzing different types of heat exchangers based
on the effectiveness-NTU method.

The following general types of heat exchangers are implemented:
1. counter-flow heat exchanger
2. parallel-flow heat exchanger
3. cross-flow heat exchanger
4. shell-and-tube heat exchanger.

The designing problem involves determining the required overall conductance UA
of a heat exchanger knowing the inlet conditions of the hot and the cold stream
and the required rate of heat transfer from the hot to the cold stream.

The analysis problem involves finding the resulting heat transfer rate from the
hot to the cold stream if the conductance UA of the heat exchanger and the inlet
conditions of the hot and cold stream are known.
"""
import math
from abc import ABC, abstractmethod
from enum import IntEnum
from hvac import Quantity

Q_ = Quantity


class AbstractHeatExchanger(ABC):
    """Abstract base class from which concrete heat exchanger models are
    derived, implementing the attributes and methods common to all heat
    exchanger models.
    """

    def __init__(
        self,
        C_cold: Quantity,
        C_hot: Quantity,
        T_cold_in: Quantity,
        T_hot_in: Quantity,
        UA: Quantity | None = None,
        Q: Quantity | None = None
    ) -> None:
        """Creates the common part of an object instantiated from any of the
        `AbstractHeatExchanger` derived classes.

        Parameters
        ----------
        C_cold:
            Capacitance rate of the cold flow stream (product of mass flow rate
            and specific heat at constant pressure). The cold flow stream
            receives heat from the hot stream and gets hotter while flowing
            through the heat exchanger.
        C_hot:
            Capacitance rate of the hot flow stream. The hot flow stream loses
            heat to the cold stream and gets cooler while flowing through the
            heat exchanger.
        T_cold_in:
            Inlet temperature of cold stream.
        T_hot_in:
            Inlet temperature of hot flow stream.
        UA: optional
            Overall conductance of the heat exchanger.
        Q: optional
            Heat transfer from hot to cold stream.
        """
        self._C_cold = C_cold.to('W / K').m
        self._C_hot = C_hot.to('W / K').m
        self._T_cold_in = T_cold_in.to('degC').m
        self._T_hot_in = T_hot_in.to('degC').m
        self._UA = UA.to('W / K').m if UA is not None else None
        self._Q = Q.to('W').m if Q is not None else None

    @property
    def C_min(self) -> Quantity:
        C_min = min(self._C_cold, self._C_hot)
        return Q_(C_min, 'W / K')

    @property
    def C_max(self) -> Quantity:
        C_max = max(self._C_cold, self._C_hot)
        return Q_(C_max, 'W / K')

    @property
    def C_r(self) -> float:
        return (self.C_min / self.C_max).m

    @property
    def NTU(self) -> float:
        """Number of transfer units."""
        if self._UA is not None:
            return self._UA / self.C_min.m
        else:
            return self.__NTU__()

    @abstractmethod
    def __NTU__(self) -> float:
        ...

    @property
    def Q_max(self) -> Quantity:
        """Maximum heat transfer from hot to cold stream
        that is theoretically possible.
        """
        Q_max = self.C_min.m * (self._T_hot_in - self._T_cold_in)
        return Q_(Q_max, 'W')

    @property
    def Q(self) -> Quantity:
        """Actual heat transfer from hot to cold stream."""
        return self.eps * self.Q_max

    @property
    def UA(self) -> Quantity:
        """Overall conductance of the heat exchanger."""
        if self._UA is not None:
            return Q_(self._UA, 'W / K')
        else:
            return Q_(self.NTU * self.C_min.m, 'W / K')

    @property
    def eps(self) -> float:
        """Effectiveness of the heat exchanger."""
        if self._Q is not None:
            return self._Q / self.Q_max.m
        else:
            return self.__eps__()

    @abstractmethod
    def __eps__(self) -> float:
        ...

    @property
    def T_hot_out(self) -> Quantity:
        """Outlet temperature of hot flow stream."""
        T_hot_out = self._T_hot_in - self.Q.m / self._C_hot
        return Q_(T_hot_out, 'degC')

    @property
    def T_cold_out(self) -> Quantity:
        """Outlet temperature of cold flow stream."""
        T_cold_out = self._T_cold_in + self.Q.m / self._C_cold
        return Q_(T_cold_out, 'degC')


class CounterFlowHeatExchanger(AbstractHeatExchanger):

    def __eps__(self) -> float:
        eps = float('nan')
        if self.C_r < 1.0:
            n = 1 - math.exp(-self.NTU * (1 - self.C_r))
            d = 1 - self.C_r * math.exp(-self.NTU * (1 - self.C_r))
            eps = n / d
        if self.C_r == 1.0:
            eps = self.NTU / (1 + self.NTU)
        return eps

    def __NTU__(self) -> float:
        NTU = float('nan')
        if self.C_r < 1.0:
            NTU = (
                math.log((1 - self.eps * self.C_r) / (1 - self.eps))
                / (1 - self.C_r)
            )
        if self.C_r == 1.0:
            NTU = self.eps / (1 - self.eps)
        return NTU


class ParallelFlowHeatExchanger(AbstractHeatExchanger):

    def __eps__(self) -> float:
        n = 1 - math.exp(-self.NTU * (1 + self.C_r))
        d = 1 + self.C_r
        eps = n / d
        return eps

    def __NTU__(self) -> float:
        n = math.log(1 - self.eps * (1 + self.C_r))
        d = 1 + self.C_r
        NTU = -n / d
        return NTU


class CrossFlowHeatExchanger(AbstractHeatExchanger):

    class Type(IntEnum):
        BOTH_FLUIDS_UNMIXED = 1  # both fluids cannot mix well in the direction perpendicular to the main flow
        BOTH_FLUIDS_MIXED = 2
        C_max_MIXED_C_min_UNMIXED = 3
        C_min_MIXED_C_max_UNMIXED = 4

    def __init__(
        self,
        type_: Type,
        C_cold: Quantity,
        C_hot: Quantity,
        T_cold_in: Quantity,
        T_hot_in: Quantity,
        UA: Quantity | None = None,
        Q: Quantity | None = None
    ) -> None:
        """Creates a `CrossFlowHeatExchanger` instance.

        Parameters
        ----------
        type_:
            The type of cross-flow heat exchanger. See IntEnum-class
            `CrossFlowHeatExchanger.Type`.

        Other Parameters
        ----------------
        See `__init__` method of class `AbstractHeatExchanger`
        """
        super().__init__(C_cold, C_hot, T_cold_in, T_hot_in, UA, Q)
        self.type = type_

    def __eps__(self) -> float:
        match self.type:
            case self.Type.BOTH_FLUIDS_UNMIXED:
                return self.__eps_01__()
            case self.Type.BOTH_FLUIDS_MIXED:
                return self.__eps_02__()
            case self.Type.C_max_MIXED_C_min_UNMIXED:
                return self.__eps_03__()
            case self.Type.C_min_MIXED_C_max_UNMIXED:
                return self.__eps_04__()

    def __eps_01__(self) -> float:
        k = (
            (self.NTU ** 0.22) / self.C_r
            * (math.exp(-self.C_r * self.NTU ** 0.78) - 1)
        )
        eps = 1 - math.exp(k)
        return eps

    def __eps_02__(self) -> float:
        a = 1 / (1 - math.exp(-self.NTU))
        b = self.C_r / (1 - math.exp(-self.C_r * self.NTU))
        c = 1 / self.NTU
        eps = 1 / (a + b - c)
        return eps

    def __eps_03__(self) -> float:
        k = self.C_r * (math.exp(-self.NTU) - 1)
        eps = (1 - math.exp(k)) / self.C_r
        return eps

    def __eps_04__(self) -> float:
        k = (1 - math.exp(-self.C_r * self.NTU)) / self.C_r
        eps = 1 - math.exp(-k)
        return eps

    def __NTU__(self) -> float:
        match self.type:
            case self.Type.BOTH_FLUIDS_UNMIXED:
                return self.__NTU_01__()
            case self.Type.BOTH_FLUIDS_MIXED:
                return self.__NTU_02__()
            case self.Type.C_max_MIXED_C_min_UNMIXED:
                return self.__NTU_03__()
            case self.Type.C_min_MIXED_C_max_UNMIXED:
                return self.__NTU_04__()

    # noinspection PyMethodMayBeStatic
    def __NTU_01__(self) -> float:
        return float('nan')

    # noinspection PyMethodMayBeStatic
    def __NTU_02__(self) -> float:
        return float('nan')

    def __NTU_03__(self) -> float:
        a = math.log(1 - self.eps * self.C_r) / self.C_r
        NTU = -math.log(1 + a)
        return NTU

    def __NTU_04__(self) -> float:
        k = self.C_r * math.log(1 - self.eps) + 1
        NTU = -math.log(k) / self.C_r
        return NTU


class ShellAndTubeHeatExchanger(AbstractHeatExchanger):

    def __init__(
        self,
        C_cold: Quantity,
        C_hot: Quantity,
        T_hot_in: Quantity,
        T_cold_in: Quantity,
        UA: Quantity | None = None,
        Q: Quantity | None = None,
        N: int = 1
    ) -> None:
        """Creates a `ShellAndTubeHeatExchanger` instance.

        Parameters
        ----------
        N: default 1
            The number of shell passes, i.e., the number of times the fluid in
            the shell passes along the length of the heat exchanger.

        Other Parameters
        ----------------
        See `__init__` method of class `AbstractHeatExchanger`.

        Notes
        -----
        The number of tube passes, i.e., the number of times the fluid in the
        tubes passes along the length of the heat exchanger, is assumed to be
        always an even number.
        """
        super().__init__(C_cold, C_hot, T_cold_in, T_hot_in, UA, Q)
        self.N = N

    def __eps__(self) -> float:
        NTU_1 = self.NTU / self.N
        k = math.exp(-NTU_1 * math.sqrt(1 + self.C_r ** 2))
        n = 1 + k
        d = 1 - k
        z = 1 + self.C_r + math.sqrt(1 + self.C_r ** 2) * (n / d)
        eps_1 = 2 / z
        k = (1 - eps_1 * self.C_r) / (1 - eps_1)
        eps = ((k ** self.N) - 1) / ((k ** self.N) - self.C_r)
        return eps

    def __NTU__(self) -> float:
        F = ((self.eps * self.C_r - 1) / (self.eps - 1)) ** (1 / self.N)
        eps_1 = (F - 1) / (F - self.C_r)
        E = (
            (2 - eps_1 * (1 + self.C_r))
            / (eps_1 * math.sqrt(1 + (self.C_r ** 2)))
        )
        n = math.log((E + 1) / (E - 1))
        d = math.sqrt(1 + (self.C_r ** 2))
        NTU_1 = n / d
        NTU = self.N * NTU_1
        return NTU
