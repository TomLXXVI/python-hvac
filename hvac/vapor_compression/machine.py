from typing import Union, Optional
import numpy as np
from scipy.optimize import fsolve
from .. import Quantity
from ..fluids import HumidAir, CP_HUMID_AIR, FluidState
from .real_compressor import FixedSpeedCompressor, VariableSpeedCompressor

Q_ = Quantity

Compressor = Union[FixedSpeedCompressor, VariableSpeedCompressor]


class CondensingUnit:

    def __init__(
        self,
        compressor: Compressor,
        eps_con: Quantity,
        air_in: HumidAir,
        ma: Quantity
    ) -> None:
        """
        Parameters
        ----------
        compressor: Compressor
            Model of the compressor.
        eps_con: Quantity
            Heat exchanger effectiveness of condenser.
        air_in: HumidAir
            State of air at condenser entrance.
        ma: Quantity
            Mass flow rate of air through condenser.
        """
        self.compressor = compressor
        self.eps_con = eps_con
        self.air_in = air_in
        self.ma = ma

    def find_Tc(self, Te: Quantity, Tc_ini: Quantity) -> Quantity:

        self.compressor.Te = Te

        def eq(unknowns: np.ndarray) -> np.ndarray:
            Tc = Q_(unknowns[0], 'K')
            self.compressor.Tc = Tc
            Qh_dot = self.compressor.Qh_dot
            Tc_new = self.air_in.Tdb + Qh_dot / (self.ma * CP_HUMID_AIR * self.eps_con)
            out = Tc_new.to('K') - Tc
            return np.array([out.m])

        roots = fsolve(eq, np.array([Tc_ini.to('K').m]))
        Tc = Q_(roots[0], 'K')
        return Tc


class SingleStageVCMachine:

    def __init__(
        self,
        compressor: Compressor,
        eps_con: Quantity,
        eps_eva: Quantity,
        lt_ma: Quantity,
        ht_ma: Quantity,
        lt_air_in: Optional[HumidAir] = None,
        ht_air_in: Optional[HumidAir] = None,
    ) -> None:
        """
        Parameters
        ----------
        compressor: Compressor
            Model of the compressor. Can be an instance of `FixSpeedCompressor`-
            or `VariableSpeedCompressor`-class.
        eps_con: Quantity
            Heat exchanger effectiveness of condenser.
        eps_eva: Quantity
            Heat exchanger effectiveness of evaporator.
        lt_ma: Quantity
            Mass flow rate of low temperature air through evaporator.
        ht_ma: Quantity
            Mass flow rate of high temperature air through condenser.
        lt_air_in: HumidAir, optional
            State of low temperature inlet air at evaporator entrance.
        ht_air_in: HumidAir, optional
            State of high temperature inlet air at condenser entrance.
        """
        self.compressor = compressor
        self.eps_con = eps_con
        self.eps_eva = eps_eva
        self.lt_ma = lt_ma
        self.ht_ma = ht_ma
        self._lt_air_in = lt_air_in
        self._ht_air_in = ht_air_in
        if lt_air_in is not None and ht_air_in is not None:
            Te_ini = lt_air_in.Tdb - Q_(10, 'K')
            Tc_ini = ht_air_in.Tdb + Q_(10, 'K')
            self._find_steady_equilibrium(Te_ini, Tc_ini)

    @property
    def lt_air_in(self) -> HumidAir:
        return self._lt_air_in

    @lt_air_in.setter
    def lt_air_in(self, lt_air_in: HumidAir) -> None:
        """Set state of air at evaporator entrance."""
        self._lt_air_in = lt_air_in

    @property
    def ht_air_in(self) -> HumidAir:
        return self._ht_air_in

    @ht_air_in.setter
    def ht_air_in(self, ht_air_in: HumidAir) -> None:
        """Set state of air at condenser entrance."""
        self._ht_air_in = ht_air_in

    @property
    def speed(self) -> Optional[Quantity]:
        if isinstance(self.compressor, VariableSpeedCompressor):
            return self.compressor.speed
        else:
            return None

    @speed.setter
    def speed(self, rpm: Quantity) -> None:
        """Set compressor speed (in case `compressor` attribute is an instance
        of `VariableSpeedCompressor`-class)."""
        if isinstance(self.compressor, VariableSpeedCompressor):
            self.compressor.speed = rpm

    def simulate(self) -> None:
        """Run the calculations to find evaporator temperature and condenser
         temperature at the current set working conditions. Call this method after
         first setting `lt_air_in`, `ht_air_in` and/or `rpm`."""
        if self._lt_air_in is not None:
            Te_ini = self._lt_air_in.Tdb - Q_(10, 'K')
        else:
            Te_ini = None
        if self._ht_air_in is not None:
            Tc_ini = self._ht_air_in.Tdb + Q_(10, 'K')
        else:
            Tc_ini = None
        if Te_ini is not None and Tc_ini is not None:
            self._find_steady_equilibrium(Te_ini, Tc_ini)

    def _find_steady_equilibrium(self, Te_ini: Quantity, Tc_ini: Quantity) -> None:
        """Find evaporator and condenser temperature under the given working
        conditions."""

        def eq(unknowns: np.ndarray) -> np.ndarray:
            Te = Q_(unknowns[0], 'K')
            Tc = Q_(unknowns[1], 'K')
            self.compressor.Te = Te
            self.compressor.Tc = Tc
            Qc_dot = self.compressor.Qc_dot
            Qh_dot = self.compressor.Qh_dot
            Te_new = self._get_Te(Qc_dot)
            Tc_new = self._get_Tc(Qh_dot)
            out1 = Te_new.to('K') - Te
            out2 = Tc_new.to('K') - Tc
            return np.array([out1.m, out2.m])

        roots = fsolve(eq, np.array([Te_ini.to('K').m, Tc_ini.to('K').m]))
        Te = Q_(roots[0], 'K')
        Tc = Q_(roots[1], 'K')
        self.compressor.Te = Te
        self.compressor.Tc = Tc

    def _get_Te(self, Qc_dot: Quantity) -> Quantity:
        """Find evaporator temperature for given cooling capacity `Qc_dot` and
        for given conditions at evaporator set at instantiation of the
        `VaporCompressionSystem`-object."""
        ha_sat = self._lt_air_in.h - Qc_dot / (self.eps_eva * self.lt_ma)
        Te = HumidAir(h=ha_sat, RH=Q_(100, 'pct')).Tdb
        return Te

    def _get_Tc(self, Qh_dot: Quantity) -> Quantity:
        """Find condenser temperature for given condenser capacity `Qh_dot` and
        for given conditions at condenser set at instantiation of the
        `VaporCompressionSystem`-object."""
        Tc = self._ht_air_in.Tdb + Qh_dot / (self.ht_ma * CP_HUMID_AIR * self.eps_con)
        return Tc

    @property
    def Qc_dot(self) -> Quantity:
        """Get cooling capacity of vapor compression system under the conditions
        set at instantiation of the `VaporCompressionSystem`-object."""
        return self.compressor.Qc_dot

    @property
    def Wc_dot(self) -> Quantity:
        """Get compressor input power under the conditions set at instantiation
        of the `VaporCompressionSystem`-object."""
        return self.compressor.Wc_dot

    @property
    def m_dot(self) -> Quantity:
        """Get refrigerant mass flow rate under the conditions set at
        instantiation of the `VaporCompressionSystem`-object."""
        return self.compressor.m_dot

    @property
    def Qh_dot(self) -> Quantity:
        """Get condenser heat rejection under the conditions
        set at instantiation of the `VaporCompressionSystem`-object."""
        return self.compressor.Qh_dot

    @property
    def COP(self) -> Quantity:
        """Get COP of vapor compression system under the conditions
        set at instantiation of the `VaporCompressionSystem`-object."""
        return self.compressor.COP

    @property
    def Te(self) -> Quantity:
        """Get evaporator temperature of vapor compression system under the
        conditions set at instantiation of the `VaporCompressionSystem`-object."""
        return self.compressor.Te

    @property
    def Tc(self) -> Quantity:
        """Get condenser temperature of vapor compression system under the
        conditions set at instantiation of the `VaporCompressionSystem`-object."""
        return self.compressor.Tc

    @property
    def suction_gas(self) -> FluidState:
        """Get state of suction gas under the conditions set at instantiation
        of the `VaporCompressionSystem`-object."""
        return self.compressor.suction_gas

    @property
    def discharge_gas(self) -> FluidState:
        """Get state of discharge gas under the conditions set at instantiation
        of the `VaporCompressionSystem`-object."""
        return self.compressor.discharge_gas

    @property
    def liquid(self) -> FluidState:
        """Get state of liquid under the conditions set at instantiation
        of the `VaporCompressionSystem`-object."""
        return self.compressor.liquid

    @property
    def mixture(self) -> FluidState:
        """Get state of mixture under the conditions set at instantiation
        of the `VaporCompressionSystem`-object."""
        return self.compressor.mixture

    @property
    def isentropic_efficiency(self) -> Quantity:
        """Get the isentropic efficiency of the compressor."""
        return self.compressor.eta_is

    @property
    def super_heat(self) -> Quantity:
        """Get the degree of superheat of the suction gas."""
        return self.compressor.dT_sh

    @property
    def sub_cooling(self) -> Quantity:
        """Get the degree of subcooling of the liquid at the inlet of
        the expansion device."""
        return self.compressor.dT_sc
