"""
Module for analyzing the steady-state operation of a single-stage vapor
compression machine.
"""
from dataclasses import dataclass
import functools
import time

import numpy as np
from scipy import optimize

from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import HumidAir, FluidState
from hvac.vapor_compression import VariableSpeedCompressor, FixedSpeedCompressor
from hvac.heat_exchanger.recuperator.fintube.continuous_fin import (
    PlainFinTubeCounterFlowAirEvaporator,
    PlainFinTubeCounterFlowAirCondenser,
    EvaporatorError, 
    CondenserError
)
from hvac.vapor_compression.refrigerant_lines import (
    SuctionLine,
    DischargeLine,
    LiquidLine
)

Q_ = Quantity
logger = ModuleLogger.get_logger(__name__)

Evaporator = PlainFinTubeCounterFlowAirEvaporator
Condenser = PlainFinTubeCounterFlowAirCondenser
Compressor = VariableSpeedCompressor | FixedSpeedCompressor


@dataclass
class Output:
    """
    Dataclass for grouping the simulation results to be returned when
    a simulation is finished.

    Use attribute `units` to set the desired measuring units of the quantities
    and the number of decimal places for rounding. Attribute `units` is a
    dictionary. The keys refer to the physical quantity and map to a tuple
    containing the unit and the number of decimal places after the decimal point.
    The default dictionary is:
    ```
        self.units: dict[str, tuple[str, int]] = {
            'm_dot': ('kg / hr', 3),
            'T': ('degC', 3),
            'W': ('g / kg', 3),
            'n': ('rpm', 3),
            'dT': ('K', 3),
            'Q_dot': ('kW', 3),
            'W_dot': ('kW', 3),
            'P_rfg': ('bar', 3),
            'P_air': ('Pa', 3),
            'h': ('kJ / kg', 3)
        }
    ```
    """
    evp_air_m_dot: Quantity | None
    cnd_air_m_dot: Quantity
    evp_air_in: HumidAir | None
    cnd_air_in: HumidAir
    n_cmp: Quantity
    dT_sh: Quantity
    evp_air_out: HumidAir | None = None
    cnd_air_out: HumidAir | None = None
    evp_Q_dot: Quantity | None = None
    cnd_Q_dot: Quantity | None = None
    cmp_W_dot: Quantity | None = None
    rfg_m_dot: Quantity | None = None
    T_evp: Quantity | None = None
    P_evp: Quantity | None = None
    T_cnd: Quantity | None = None
    P_cnd: Quantity | None = None
    dT_sc: Quantity | None = None
    suction_gas: FluidState | None = None
    discharge_gas: FluidState | None = None
    liquid: FluidState | None = None
    mixture: FluidState | None = None
    evp_eps: Quantity | None = None
    cnd_eps: Quantity | None = None
    COP: Quantity | None = None
    EER: Quantity | None = None
    evp_air_dP: Quantity | None = None
    cnd_air_dP: Quantity | None = None

    def __post_init__(self):
        self.units: dict[str, tuple[str, int]] = {
            'm_dot': ('kg / hr', 3),
            'T': ('degC', 3),
            'W': ('g / kg', 3),
            'n': ('rpm', 3),
            'dT': ('K', 3),
            'Q_dot': ('kW', 3),
            'W_dot': ('kW', 3),
            'P_rfg': ('bar', 3),
            'P_air': ('Pa', 3),
            'h': ('kJ / kg', 3)
        }
        self.success: bool = False

    @classmethod
    def create_from_dict(cls, d: dict) -> 'Output':
        obj = cls(**d)
        return obj

    def to_dict(self) -> dict[str, float]:
        u_m_dot, d_m_dot = self.units['m_dot']
        u_T, d_T = self.units['T']
        u_W, d_W = self.units['W']
        u_n, d_n = self.units['n']
        u_dT, d_dT = self.units['dT']
        u_Q_dot, d_Q_dot = self.units['Q_dot']
        u_W_dot, d_W_dot = self.units['W_dot']
        u_P_rfg, d_P_rfg = self.units['P_rfg']
        u_P_air, d_P_air = self.units['P_air']
        u_h, d_h = self.units['h']

        if self.success:
            d = {
                'evp_air_m_dot': self.evp_air_m_dot.to(u_m_dot).m,
                'cnd_air_m_dot': self.cnd_air_m_dot.to(u_m_dot).m,
                'T_evp_air_in': self.evp_air_in.Tdb.to(u_T).m,
                'W_evp_air_in': self.evp_air_in.W.to(u_W).m,
                'RH_evp_air_in': self.evp_air_in.RH.to('pct').m,
                'T_cnd_air_in': self.cnd_air_in.Tdb.to(u_T).m,
                'W_cnd_air_in': self.cnd_air_in.W.to(u_W).m,
                'RH_cnd_air_in': self.cnd_air_in.RH.to('pct').m,
                'n_cmp': self.n_cmp.to(u_n).m if self.n_cmp is not None else Q_(float('nan'), u_n),
                'dT_sh': self.dT_sh.to(u_dT).m,
                'T_evp_air_out': self.evp_air_out.Tdb.to(u_T).m,
                'W_evp_air_out': self.evp_air_out.W.to(u_W).m,
                'RH_evp_air_out': self.evp_air_out.RH.to('pct').m,
                'T_cnd_air_out': self.cnd_air_out.Tdb.to(u_T).m,
                'W_cnd_air_out': self.cnd_air_out.W.to(u_W).m,
                'RH_cnd_air_out': self.cnd_air_out.RH.to('pct').m,
                'evp_Q_dot': self.evp_Q_dot.to(u_Q_dot).m,
                'cnd_Q_dot': self.cnd_Q_dot.to(u_Q_dot).m,
                'cmp_W_dot': self.cmp_W_dot.to(u_W_dot).m,
                'rfg_m_dot': self.rfg_m_dot.to(u_m_dot).m,
                'COP': self.COP.to('frac').m,
                'EER': self.EER.to('frac').m,
                'evp_eps': self.evp_eps.to('frac').m,
                'cnd_eps': self.cnd_eps.to('frac').m,
                'T_evp': self.T_evp.to(u_T).m,
                'P_evp': self.P_evp.to(u_P_rfg).m,
                'T_cnd': self.T_cnd.to(u_T).m,
                'P_cnd': self.P_cnd.to(u_P_rfg).m,
                'dT_sc': self.dT_sc.to(u_dT).m,
                'T_suction_gas': self.suction_gas.T.to(u_T).m,
                'P_suction_gas': self.suction_gas.P.to(u_P_rfg).m,
                'h_suction_gas': self.suction_gas.h.to(u_h).m,
                'T_discharge_gas': self.discharge_gas.T.to(u_T).m,
                'P_discharge_gas': self.discharge_gas.P.to(u_P_rfg).m,
                'h_discharge_gas': self.discharge_gas.h.to(u_h).m,
                'T_liquid': self.liquid.T.to(u_T).m,
                'P_liquid': self.liquid.P.to(u_P_rfg).m,
                'h_liquid': self.liquid.h.to(u_h).m,
                'T_mixture': self.mixture.T.to(u_T).m,
                'P_mixture': self.mixture.P.to(u_P_rfg).m,
                'h_mixture': self.mixture.h.to(u_h).m,
                'x_mixture': self.mixture.x.to('pct').m,
                'evp_air_dP': self.evp_air_dP.to(u_P_air).m,
                'cnd_air_dP': self.cnd_air_dP.to(u_P_air).m
            }
        else:
            d = {
                'evp_air_m_dot': self.evp_air_m_dot.to(u_m_dot).m,
                'cnd_air_m_dot': self.cnd_air_m_dot.to(u_m_dot).m,
                'T_evp_air_in': self.evp_air_in.Tdb.to(u_T).m,
                'W_evp_air_in': self.evp_air_in.W.to(u_W).m,
                'RH_evp_air_in': self.evp_air_in.RH.to('pct').m,
                'T_cnd_air_in': self.cnd_air_in.Tdb.to(u_T).m,
                'W_cnd_air_in': self.cnd_air_in.W.to(u_W).m,
                'RH_cnd_air_in': self.cnd_air_in.RH.to('pct').m,
                'n_cmp': self.n_cmp.to(u_n).m,
                'dT_sh': self.dT_sh.to(u_dT).m,
                'T_evp_air_out': float('nan'),
                'W_evp_air_out': float('nan'),
                'RH_evp_air_out': float('nan'),
                'T_cnd_air_out': float('nan'),
                'W_cnd_air_out': float('nan'),
                'RH_cnd_air_out': float('nan'),
                'evp_Q_dot': float('nan'),
                'cnd_Q_dot': float('nan'),
                'cmp_W_dot': float('nan'),
                'rfg_m_dot': float('nan'),
                'COP': float('nan'),
                'EER': float('nan'),
                'evp_eps': float('nan'),
                'cnd_eps': float('nan'),
                'T_evp': float('nan'),
                'P_evp': float('nan'),
                'T_cnd': float('nan'),
                'P_cnd': float('nan'),
                'dT_sc': float('nan'),
                'T_suction_gas': float('nan'),
                'P_suction_gas': float('nan'),
                'h_suction_gas': float('nan'),
                'T_discharge_gas': float('nan'),
                'P_discharge_gas': float('nan'),
                'h_discharge_gas': float('nan'),
                'T_liquid': float('nan'),
                'P_liquid': float('nan'),
                'h_liquid': float('nan'),
                'T_mixture': float('nan'),
                'P_mixture': float('nan'),
                'h_mixture': float('nan'),
                'x_mixture': float('nan'),
                'evp_air_dP': float('nan'),
                'cnd_air_dP': float('nan')
            }
        return d

    def to_text(self) -> str:
        if self.success:
            u_m_dot, d_m_dot = self.units['m_dot']
            u_T, d_T = self.units['T']
            u_W, d_W = self.units['W']
            u_n, d_n = self.units['n']
            u_dT, d_dT = self.units['dT']
            u_Q_dot, d_Q_dot = self.units['Q_dot']
            u_W_dot, d_W_dot = self.units['W_dot']
            u_P_rfg, d_P_rfg = self.units['P_rfg']
            u_P_air, d_P_air = self.units['P_air']
            u_h, d_h = self.units['h']

            output = ''
            if self.evp_air_m_dot is not None:
                output += f"evp_air_m_dot = {self.evp_air_m_dot.to(u_m_dot):~P.{d_m_dot}f}\n"
            if self.cnd_air_m_dot is not None:
                output += f"cnd_air_m_dot = {self.cnd_air_m_dot.to(u_m_dot):~P.{d_m_dot}f}\n"
            if self.evp_air_in is not None:
                output += (
                    "evp_air_in = "
                    f"{self.evp_air_in.Tdb.to(u_T):~P.{d_T}f} DB, "
                    f"{self.evp_air_in.W.to(u_W):~P.{d_W}f} AH "
                    f"({self.evp_air_in.RH.to('pct'):~P.0f} RH)\n"
                )
            if self.cnd_air_in is not None:
                output += (
                    "cnd_air_in = "
                    f"{self.cnd_air_in.Tdb.to(u_T):~P.{d_T}f} DB, "
                    f"{self.cnd_air_in.W.to(u_W):~P.{d_W}f} AH "
                    f"({self.cnd_air_in.RH.to('pct'):~P.0f} RH)\n"
                )
            if self.n_cmp is not None:
                output += f"n_cmp = {self.n_cmp.to(u_n):~P.{d_n}f}\n"
            if self.dT_sh is not None:
                output += f"dT_sh = {self.dT_sh.to(u_dT):~P.{d_dT}f}\n"
            if self.evp_air_out is not None:
                output += (
                    "evp_air_out = "
                    f"{self.evp_air_out.Tdb.to(u_T):~P.{d_T}f} DB, "
                    f"{self.evp_air_out.W.to(u_W):~P.{d_W}f} AH "
                    f"({self.evp_air_out.RH.to('pct'):~P.0f} RH)\n"
                )
            if self.cnd_air_out is not None:
                output += (
                    "cnd_air_out = "
                    f"{self.cnd_air_out.Tdb.to(u_T):~P.{d_T}f} DB, "
                    f"{self.cnd_air_out.W.to(u_W):~P.{d_W}f} AH "
                    f"({self.cnd_air_out.RH.to('pct'):~P.0f} RH)\n"
                )
            if self.evp_Q_dot is not None:
                output += f"evp_Q_dot = {self.evp_Q_dot.to(u_Q_dot):~P.{d_Q_dot}f}\n"
            if self.cnd_Q_dot is not None:
                output += f"cnd_Q_dot = {self.cnd_Q_dot.to(u_Q_dot):~P.{d_Q_dot}f}\n"
            if self.cmp_W_dot is not None:
                output += f"cmp_W_dot = {self.cmp_W_dot.to(u_W_dot):~P.{d_W_dot}f}\n"
            if self.rfg_m_dot is not None:
                output += f"rfg_m_dot = {self.rfg_m_dot.to(u_m_dot):~P.{d_m_dot}f}\n"
            if self.COP is not None:
                output += f"COP = {self.COP.to('frac'):~P.3f}\n"
            if self.EER is not None:
                output += f"EER = {self.EER.to('frac'):~P.3f}\n"
            if self.evp_eps is not None:
                output += f"evp_eps = {self.evp_eps.to('frac'):~P.4f}\n"
            if self.cnd_eps is not None:
                output += f"cnd_eps = {self.cnd_eps.to('frac'):~P.4f}\n"
            if self.T_evp is not None:
                output += f"T_evp = {self.T_evp.to(u_T):~P.{d_T}f}\n"
            if self.P_evp is not None:
                output += f"P_evp = {self.P_evp.to(u_P_rfg):~P.{d_P_rfg}f}\n"
            if self.T_cnd is not None:
                output += f"T_cnd = {self.T_cnd.to(u_T):~P.{d_T}f}\n"
            if self.P_cnd is not None:
                output += f"P_cnd = {self.P_cnd.to(u_P_rfg):~P.{d_P_rfg}f}\n"
            if self.dT_sc is not None:
                output += f"dT_sc = {self.dT_sc.to(u_dT):~P.{d_dT}f}\n"
            if self.suction_gas is not None:
                output += (
                    "suction_gas = "
                    f"{self.suction_gas.T.to(u_T):~P.{d_T}f}, "
                    f"{self.suction_gas.h.to(u_h):~P.{d_h}f}, "
                    f"{self.suction_gas.P.to(u_P_rfg):~P.{d_P_rfg}f}\n"
                )
            if self.discharge_gas is not None:
                output += (
                    "discharge_gas = "
                    f"{self.discharge_gas.T.to(u_T):~P.{d_T}f}, "
                    f"{self.discharge_gas.h.to(u_h):~P.{d_h}f}, "
                    f"{self.discharge_gas.P.to(u_P_rfg):~P.{d_P_rfg}f}\n"
                )
            if self.liquid is not None:
                output += (
                    "liquid = "
                    f"{self.liquid.T.to(u_T):~P.{d_T}f}, "
                    f"{self.liquid.h.to(u_h):~P.{d_h}f}, "
                    f"{self.liquid.P.to(u_P_rfg):~P.{d_P_rfg}f}\n"
                )
            if self.mixture is not None:
                output += (
                    "mixture = "
                    f"{self.mixture.T.to(u_T):~P.{d_T}f}, "
                    f"{self.mixture.h.to(u_h):~P.{d_h}f}, "
                    f"{self.mixture.P.to(u_P_rfg):~P.{d_P_rfg}f}, "
                    f"{self.mixture.x.to('pct'):~P.0f}\n"
                )
            if self.evp_air_dP is not None:
                output += f"evp_air_dP = {self.evp_air_dP.to(u_P_air):~P.{d_P_air}f}\n"
            if self.cnd_air_dP is not None:
                output += f"cnd_air_dP = {self.cnd_air_dP.to(u_P_air):~P.{d_P_air}f}\n"
            return output
        else:
            return 'None'


class SingleStageVaporCompressionMachine:
    """
    Model of a single-stage vapor compression machine.

    Method `rate` can be used to retrieve the steady-state machine performance
    at a given compressor speed and for given operating conditions on the
    air-side of the evaporator and condenser. The method searches for the
    evaporation and condensation temperature for which the mass flow rate let
    through by the expansion device (to maintain the set degree of refrigerant
    superheating at the evaporator outlet) balances with the mass flow rate of
    refrigerant displaced by the compressor.

    Method `balance_by_speed` can be used to retrieve the steady-state machine
    performance at a given evaporation temperature and condensation temperature
    and at given operating conditions on the air-side of the evaporator and
    condenser. The method searches for the compressor speed for which the mass
    flow rate let through by the expansion device (to maintain the set degree
    of refrigerant superheating at the evaporator outlet) balances with the mass
    flow rate of refrigerant displaced by the compressor. This method can only
    be used if the compressor has a variable speed drive.
    """
    def __init__(
        self,
        evaporator: Evaporator,
        condenser: Condenser,
        compressor: Compressor,
        dT_sh: Quantity | None = None,
        n_cmp_min: Quantity | None = None,
        n_cmp_max: Quantity | None = None,
        suction_line: SuctionLine | None = None,
        discharge_line: DischargeLine | None = None,
        liquid_line: LiquidLine | None = None,
        solver: str = "least_squares"
    ) -> None:
        """
        Creates a `SingleStageVaporCompressionMachine` object.

        Parameters
        ----------
        evaporator:
            A model instance of the evaporator.
        condenser:
            A model instance of the condenser.
        compressor:
            A model instance of the compressor.
        dT_sh: optional
            The setting on the expansion device of the degree of refrigerant
            superheating. If `None`, the same setting as on the compressor is
            taken.
        n_cmp_min: optional
            The minimum speed that can be set on the compressor if it is a
            variable speed compressor.
        n_cmp_max: optional
            The maximum speed that can be set on the compressor if it is a
            variable speed compressor.
        suction_line: optional
            Suction line between evaporator and compressor.
        discharge_line: optional
            Discharge line between compressor and condenser.
        liquid_line: optional
            Liquid line between condenser and expansion device.
        solver: str, {"minimize", "least_squares" (default)}
            Indicates the solving technique used in method `rate()`: either
            Scipy `minimize()` or Scipy `least_squares()`, which is the default
            and fastest solving technique.
        """
        self.evaporator = evaporator
        self.condenser = condenser
        self.compressor = compressor
        self.refrigerant = compressor.refrigerant
        self.dT_sh = dT_sh if dT_sh is not None else compressor.dT_sh
        self.n_cmp_min = n_cmp_min
        self.n_cmp_max = n_cmp_max

        self.suction_line = suction_line
        self.discharge_line = discharge_line
        self.liquid_line = liquid_line
        self._refrigeration_lines: bool = False
        if all([self.suction_line, self.discharge_line, self.liquid_line]):
            self._refrigeration_lines = True

        self.evp_air_in: HumidAir | None = None
        self.evp_air_m_dot: Quantity | None = None
        self.cnd_air_in: HumidAir | None = None
        self.cnd_air_m_dot: Quantity | None = None
        self.n_cmp: Quantity | None = None

        self.output: Output | None = None

        if solver not in {"least_squares", "minimize"}:
            raise ValueError("solver must be 'least_squares' or 'minimize'")
        self._solver = solver

    @staticmethod
    def time_it(fun):
        @functools.wraps(fun)
        def wrapper(self, *args, **kwargs):
            t_start = time.perf_counter()

            output = fun(self, *args, **kwargs)

            t_finish = time.perf_counter()
            dt_exec = t_finish - t_start

            logger.info(
                f"Execution time: {dt_exec} seconds"
            )

            if output.success:
                mb_err = self._check_mass_balance()
                eb_err = self._check_energy_balance()
                logger.info(
                    "Error mass balance: "
                    f"absolute error = {mb_err[0].to('kg / hr'):~P.3f}, "
                    f"relative error = {mb_err[1].to('pct'):~P.2f}"
                )
                logger.info(
                    "Error energy balance: "
                    f"absolute error = {eb_err[0].to('kW'):~P.3f}, "
                    f"relative error = {eb_err[1].to('pct'):~P.2f}"
                )
            return output
        return wrapper

    @time_it
    def rate(
        self,
        evp_air_in: HumidAir,
        evp_air_m_dot: Quantity,
        cnd_air_in: HumidAir,
        cnd_air_m_dot: Quantity,
        n_cmp: Quantity | None = None,
        T_evp_ini: Quantity | None = None,
        T_cnd_ini: Quantity | None = None,
        solver_options: dict | None = None
    ) -> Output:
        """
        Determines the steady-state performance of the single-stage vapor
        compression machine for the given operating conditions.

        Parameters
        ----------
        evp_air_in: HumidAir
            State of air entering evaporator.
        evp_air_m_dot: Quantity
            Mass flow rate of air through evaporator.
        cnd_air_in: HumidAir
            State of air entering condenser.
        cnd_air_m_dot: Quantity
            Mass flow rate of air through condenser.
        n_cmp: Quantity, optional
            Current compressor speed (only to be used with a variable speed
            compressor).
        T_evp_ini: Quantity, optional
            Initial guess for the evaporation temperature.
        T_cnd_ini: Quantity, optional
            Initial guess for the condensation temperature.
        solver_options: dict[str, int | float], default dict(maxiter=50, maxfev=50, xatol=0.1, fatol=0.05)
            xa_tol:
                Absolute tolerance for the change in evaporation and condensation
                temperature used in the minimization algorithm (see parameter
                `xatol` in `scipy.optimize.minimize` with Nelder-Mead method).
            fa_tol:
                Absolute tolerance for the change in the function value used in the
                minimization algorithm (see parameter `fatol` in
                `scipy.optimize.minimize` with Nelder-Mead method).
            maxiter, maxfev:
                Maximum number of iterations and function evaluations used in the
                minimization algorithm.

        Returns
        -------
        Object of class `Output`.

        Notes
        -----
        This method uses the Nelder-Mead minimization algorithm to determine the
        evaporation and condensation temperature in steady-state machine
        operation under the given operating conditions.
        The algorithm tries to find an evaporation and condensation temperature
        for which the difference is minimal between the mass flow rate of
        refrigerant according to the compressor model and the mass flow rate
        of refrigerant according to the evaporator model (where the expansion
        device regulates the mass flow rate of refrigerant in order to maintain
        the set degree of refrigerant superheating at the evaporator outlet).
        In the end, the mass flow rate of refrigerant displaced by the 
        compressor should balance with the mass flow rate let through by the
        expansion device.
        
        1. The algorithm starts with an initial guess for the evaporation
           temperature and the condensing temperature.
        2. With these two values, the refrigerant mass flow rate displaced by
           the compressor is determined with the compressor model, and also the
           state of the discharge gas, which enters the condenser.
        3. At the condenser, the state of entering air and the air mass flow
           rate are fixed. With the state of the entering refrigerant and the
           refrigerant mass flow rate, a solution is determined from the
           condenser model for the condenser's performance. In particular, the 
           state of refrigerant leaving the condenser and the state of air 
           leaving the condenser are determined.
        4. Considering that the expansion process is an isenthalpic process,
           the state of refrigerant entering the evaporator is now also 
           determined.
        5. At the evaporator, the state of entering air, the air mass flow
           rate, and the degree of refrigerant superheating (being a setting on
           the expansion device) are fixed. The refrigerant mass flow rate is
           determined such that the refrigerant leaves the evaporator with the
           degree of superheating set on the expansion device.
        6. Ultimately, the refrigerant mass flow rate let through by the
           expansion device and the refrigerant mass flow rate displaced by the
           compressor must balance. The deviation is determined between the
           refrigerant mass flow rate let through by the expansion device and
           the mass flow rate displaced by the compressor.
        7. Based on the absolute value of the deviation, the minimization
           algorithm determines a new value for the evaporating temperature and
           a new value for the condensing temperature.
        8. Steps 2 to 7 are repeated until the deviation between the two mass
           flow rates has become sufficiently small, or until the maximum number
           of iterations has been reached.
        """
        self._init(
            evp_air_in, evp_air_m_dot,
            cnd_air_in, cnd_air_m_dot,
            n_cmp  # <-- compressor speed is fixed
        )

        T_evp_ini = self._guess_initial_T_evp(T_evp_ini)
        T_cnd_ini = self._guess_initial_T_cnd(T_cnd_ini)
        bounds = self._get_physical_bounds()
        counter = [0]

        logger.info(
            "Starting machine operation analysis with..."
        )
        if n_cmp is not None:
            logger.info(
                f"n_cmp = {n_cmp:~P.0f}"
            )
        logger.info(
            f"evp_air_in: {evp_air_in.Tdb.to('degC'):~P.2f}, "
            f"{evp_air_in.RH.to('pct'):~P.0f}"
        )
        logger.info(
            f"evp_air_m_dot: {evp_air_m_dot.to('kg / hr'):~P.2f}"
        )
        logger.info(
            f"cnd_air_in: {cnd_air_in.Tdb.to('degC'):~P.2f}, "
            f"{cnd_air_in.RH.to('pct'):~P.0f}"
        )
        logger.info(
            f"cnd_air_m_dot: {cnd_air_m_dot.to('kg / hr'):~P.2f}"
        )
        try:
            if self._solver == "minimize":
                if solver_options is None:
                    dict(
                        maxiter=50,
                        maxfev=50,
                        xatol=0.1,
                        fatol=0.05
                    )
                res = optimize.minimize(
                    self.__fun_rate__,
                    args=(counter,),
                    x0=np.array([T_evp_ini, T_cnd_ini]),
                    method='Nelder-Mead',
                    bounds=bounds,
                    options=solver_options
                )
            else:
                res = optimize.least_squares(
                    self.__residuals__,
                    x0=np.array([T_evp_ini, T_cnd_ini]),
                    args=(counter,),
                    bounds=bounds,
                    method="trf",
                    xtol=1e-2,
                    diff_step=np.array([5e-2, 5e-2]),  # 0.05 °C
                    ftol=1e-8,
                    gtol=1e-8
                )
        except Exception as err:
            logger.error(
                f'Analysis failed due to "{type(err).__name__}: {err}".'
            )
            self.output = Output(
                evp_air_m_dot=self.evp_air_m_dot,
                cnd_air_m_dot=self.cnd_air_m_dot,
                evp_air_in=self.evp_air_in,
                cnd_air_in=self.cnd_air_in,
                n_cmp=self.n_cmp,
                dT_sh=self.dT_sh
            )
            self.output.success = False
            return self.output
        else:
            logger.info(
                f'Analysis finished after {counter[0]} iterations with '
                f'message: "{res.message}"'
            )
            # Update the machine operating state using the evaporation and
            # condensation temperature returned from `optimize.minimize` (these
            # are not necessarily the values that were used in the last
            # iteration).
            self.compressor.T_evp = Q_(res.x[0], 'degC')
            self.compressor.T_cnd = Q_(res.x[1], 'degC')
            cmp_rfg_m_dot = self.compressor.m_dot.to('kg / hr')
            dev = self._get_deviation(cmp_rfg_m_dot, 0, logger_on=False)
            if np.isinf(dev):
                logger.error("Machine performance could not be determined.")
                self.output = Output(
                    evp_air_m_dot=self.evp_air_m_dot,
                    cnd_air_m_dot=self.cnd_air_m_dot,
                    evp_air_in=self.evp_air_in,
                    cnd_air_in=self.cnd_air_in,
                    n_cmp=self.n_cmp,
                    dT_sh=self.dT_sh
                )
                self.output.success = False
                return self.output
            else:
                self.output = Output(
                    evp_air_m_dot=self.evp_air_m_dot,
                    cnd_air_m_dot=self.cnd_air_m_dot,
                    evp_air_in=self.evp_air_in,
                    cnd_air_in=self.cnd_air_in,
                    n_cmp=self.n_cmp,
                    dT_sh=self.dT_sh,
                    evp_air_out=self.evaporator.air_out,
                    cnd_air_out=self.condenser.air_out,
                    evp_Q_dot=self.evaporator.Q_dot,
                    cnd_Q_dot=self.condenser.Q_dot,
                    cmp_W_dot=self.compressor.W_dot,
                    rfg_m_dot=self.compressor.m_dot,
                    T_evp=self.evaporator.T_evp,
                    P_evp=self.evaporator.P_evp,
                    T_cnd=self.condenser.T_cnd,
                    P_cnd=self.condenser.P_cnd,
                    dT_sc=self.condenser.dT_sc,
                    suction_gas=self.evaporator.rfg_out,
                    discharge_gas=self.condenser.rfg_in,
                    liquid=self.condenser.rfg_out,
                    mixture=self.evaporator.rfg_in,
                    COP=self.condenser.Q_dot / self.compressor.W_dot,
                    EER=self.evaporator.Q_dot / self.compressor.W_dot,
                    evp_eps=self.evaporator.eps,
                    cnd_eps=self.condenser.eps,
                    evp_air_dP=self.evaporator.air_dP,
                    cnd_air_dP=self.condenser.air_dP
                )
                self.output.success = True
                return self.output

    def __fun_rate__(self, unknowns: np.ndarray, counter: list[int]) -> float:
        """
        Calculates the deviation between the refrigerant mass flow rate let
        through by the expansion device to maintain the set degree of
        refrigerant superheating, and the refrigerant mass flow rate generated
        by the compressor at the given evaporation temperature and condensing
        temperature.

        Parameters
        ----------
        unknowns:
            Numpy array with two floats: the evaporation temperature, and the
            condensing temperature.
        counter:
            List of a single int to count the number of calls to `__fun__`.

        Notes
        -----
        This function should only be used internally in method `rate`.
        """
        i = counter[0]
        T_evp = Q_(unknowns[0], 'degC')
        T_cnd = Q_(unknowns[1], 'degC')

        logger.info(
            f"Iteration {i + 1}: "
            f"Try with: T_evp = {T_evp:~P.3f}, T_cnd = {T_cnd:~P.3f}"
        )

        self.compressor.T_evp = T_evp  # saturation temperature at compressor inlet
        self.compressor.T_cnd = T_cnd  # saturation temperature at compressor outlet
        cmp_rfg_m_dot = self.compressor.m_dot.to('kg / hr')
        dev = self._get_deviation(cmp_rfg_m_dot, i)
        counter[0] += 1
        return abs(dev)

    def __residuals__(self, unknowns: np.ndarray, counter: list[int]) -> np.ndarray:
        i = counter[0]
        T_evp = Q_(unknowns[0], 'degC')
        T_cnd = Q_(unknowns[1], 'degC')

        logger.info(
            f"Iteration {i + 1}: "
            f"Try with: T_evp = {T_evp:~P.3f}, T_cnd = {T_cnd:~P.3f}"
        )

        self.compressor.T_evp = T_evp  # saturation temperature at compressor inlet
        self.compressor.T_cnd = T_cnd  # saturation temperature at compressor outlet
        cmp_rfg_m_dot = self.compressor.m_dot.to('kg / hr')

        dm = self._get_deviation(cmp_rfg_m_dot, i)  # deviation on mass balance
        try:
            de = (                                      # deviation on energy balance
                self.evaporator.Q_dot.to('kW').m
                + self.compressor.W_dot.to('kW').m
                - self.condenser.Q_dot.to('kW').m
            )
        except:
            de = float('inf')

        dm_scale = 1.0   # kg/hr
        de_scale = 0.05  # kW

        counter[0] += 1
        return np.array([dm/dm_scale, de/de_scale], dtype=float)

    @time_it
    def balance_by_speed(
        self,
        evp_air_in: HumidAir,
        evp_air_m_dot: Quantity,
        cnd_air_in: HumidAir,
        cnd_air_m_dot: Quantity,
        T_evp: Quantity,
        T_cnd: Quantity
    ) -> Output:
        """
        Determines the compressor speed for which the mass flow rate of
        refrigerant displaced by the compressor balances the mass flow rate of
        refrigerant let through by the expansion device, while maintaining the
        set degree of refrigerant superheating at the given evaporation and
        condensing temperature.

        Parameters
        ----------
        evp_air_in:
            State of air entering evaporator.
        evp_air_m_dot:
            Mass flow rate of air through evaporator.
        cnd_air_in:
            State of air entering condenser.
        cnd_air_m_dot:
            Mass flow rate of air through condenser.
        T_evp:
            Evaporation temperature at which the mass flow rates must be
            balanced.
        T_cnd:
            Condensing temperature at which the mass flow rates must be
            balanced.

        Notes
        -----
        This method can only be used with a variable speed compressor.
        Also, the minimum and maximum compressor speed must be set when
        instantiating the `SingleStageVaporCompressionMachine` object.

        Returns
        -------
        Object of class `Output`.
        """
        c1 = isinstance(self.compressor, VariableSpeedCompressor)
        c2 = (self.n_cmp_min is not None) and (self.n_cmp_max is not None)
        if c1 and c2:
            logger.info(
                "Starting speed balancing with..."
            )
            logger.info(
                f"evp_air_in: {evp_air_in.Tdb.to('degC'):~P.2f}, "
                f"{evp_air_in.RH.to('pct'):~P.0f}"
            )
            logger.info(
                f"evp_air_m_dot: {evp_air_m_dot:~P.2f}"
            )
            logger.info(
                f"cnd_air_in: {cnd_air_in.Tdb.to('degC'):~P.2f}, "
                f"{cnd_air_in.RH.to('pct'):~P.0f}"
            )
            logger.info(
                f"cnd_air_m_dot: {cnd_air_m_dot:~P.2f}"
            )
            logger.info(
                f"T_evp = {T_evp:~P.2f}"
            )
            logger.info(
                f"T_cnd = {T_cnd:~P.2f}"
            )

            self._init(
                evp_air_in, evp_air_m_dot,
                cnd_air_in, cnd_air_m_dot,
                T_evp=T_evp, T_cnd=T_cnd
            )

            # Find balanced speed between min and max compressor speed:
            n_cmp_min = self.n_cmp_min.to('rpm').m
            n_cmp_max = self.n_cmp_max.to('rpm').m
            counter = [0]
            f_cont, f_int = self._make_obj(counter, self.__fun_balance__)
            try:
                res = optimize.minimize_scalar(
                    f_cont,
                    bounds=(n_cmp_min, n_cmp_max),
                    method="bounded",
                    options={"xatol": 1.0}
                )
                best_rpm = int(round(res.x))
                self.n_cmp = Q_(best_rpm, 'rpm')
            except Exception as err:
                if isinstance(err, ValueError):
                    logger.error(
                        "No compressor speed found to establish the "
                        "given evaporation temperature "
                        f"{self.compressor.T_evp.to('degC'):~P.3f} and condensing "
                        f"temperature {self.compressor.T_cnd.to('degC'):~P.3f} "
                        f"under given operating conditions."
                    )
                else:
                    logger.error(
                        f'Speed balancing failed due to '
                        f'"{type(err).__name__}: {err}".'
                    )
                self.output = Output(
                    evp_air_m_dot=self.evp_air_m_dot,
                    cnd_air_m_dot=self.cnd_air_m_dot,
                    evp_air_in=self.evp_air_in,
                    cnd_air_in=self.cnd_air_in,
                    n_cmp=self.n_cmp,
                    dT_sh=self.dT_sh
                )
                self.output.success = False
                return self.output
            else:
                logger.info(
                    f"Speed balancing finished after {counter[0]} iterations with "
                    f"message: {res.message}"
                )

                f_int(best_rpm)

                self.output = Output(
                    evp_air_m_dot=self.evp_air_m_dot,
                    cnd_air_m_dot=self.cnd_air_m_dot,
                    evp_air_in=self.evp_air_in,
                    cnd_air_in=self.cnd_air_in,
                    n_cmp=self.n_cmp,
                    dT_sh=self.dT_sh,
                    evp_air_out=self.evaporator.air_out,
                    cnd_air_out=self.condenser.air_out,
                    evp_Q_dot=self.evaporator.Q_dot,
                    cnd_Q_dot=self.condenser.Q_dot,
                    cmp_W_dot=self.compressor.W_dot,
                    rfg_m_dot=self.compressor.m_dot,
                    T_evp=self.evaporator.T_evp,
                    P_evp=self.evaporator.P_evp,
                    T_cnd=self.condenser.T_cnd,
                    P_cnd=self.condenser.P_cnd,
                    dT_sc=self.condenser.dT_sc,
                    suction_gas=self.evaporator.rfg_out,
                    discharge_gas=self.condenser.rfg_in,
                    liquid=self.condenser.rfg_out,
                    mixture=self.evaporator.rfg_in,
                    COP=self.condenser.Q_dot / self.compressor.W_dot,
                    EER=self.evaporator.Q_dot / self.compressor.W_dot,
                    evp_eps=self.evaporator.eps,
                    cnd_eps=self.condenser.eps,
                    evp_air_dP=self.evaporator.air_dP,
                    cnd_air_dP=self.condenser.air_dP
                )
                self.output.success = True
                return self.output
        elif not c1:
            raise ValueError(
                "The compressor is without variable speed drive."
            )
        elif not c2:
            raise ValueError(
                "The speed range of the compressor has not been specified."
            )
        return None

    def __fun_balance__(self, n_cmp: float, counter: list[int]) -> float:
        """
        Calculates the deviation between the refrigerant mass flow rate let
        through by the expansion device to maintain the set degree
        of refrigerant superheating, and the refrigerant mass flow rate
        generated by the compressor at the given compressor speed.

        Parameters
        ----------
        n_cmp:
            Compressor speed
        counter:
            List of a single int to count the number of calls to `__fun__`.

        Notes
        -----
        This function should only be used internally in method `balance_by_speed`.
        """
        i = counter[0]

        logger.info(
            f"Iteration {i + 1}: "
            f"Try with: {n_cmp:.3f} rpm"
        )

        self.compressor.speed = Q_(n_cmp, 'rpm')
        cmp_rfg_m_dot = self.compressor.m_dot.to('kg / hr')
        dev = self._get_deviation(cmp_rfg_m_dot, i)
        counter[0] += 1
        return abs(dev)

    @staticmethod
    def _make_obj(counter, fun_balance):
        @functools.lru_cache(maxsize=None)
        def f_int(rpm_int: int) -> float:
            return fun_balance(rpm_int, counter)

        def f_cont(rpm: float) -> float:
            rpm_i = int(round(rpm))
            return f_int(rpm_i)

        return f_cont, f_int

    def _init(
        self,
        evp_air_in: HumidAir,
        evp_air_m_dot: Quantity,
        cnd_air_in: HumidAir,
        cnd_air_m_dot: Quantity,
        n_cmp: Quantity | None = None,
        T_evp: Quantity | None = None,
        T_cnd: Quantity | None = None
    ) -> None:
        """
        Sets the operating conditions of the machine, when calling method `rate`
        or method `balance_by_speed`.
        """
        self.evp_air_in = evp_air_in
        self.evp_air_m_dot = evp_air_m_dot
        self.cnd_air_in = cnd_air_in
        self.cnd_air_m_dot = cnd_air_m_dot
        self.compressor.dT_sh = self.dT_sh
        if isinstance(self.compressor, VariableSpeedCompressor) and n_cmp is not None:
            self.n_cmp = self.compressor.speed = n_cmp
        if T_evp is not None:
            self.compressor.T_evp = T_evp
        if T_cnd is not None:
            self.compressor.T_cnd = T_cnd

    def _guess_initial_T_evp(self, T_evp_ini: Quantity | None) -> float:
        """
        If the user does not provide an initial guess for the evaporation
        temperature to method `rate`, an initial value of 20 K below the
        temperature of the air entering the evaporator is taken.
        """
        if T_evp_ini is None:
            dT = Q_(20, 'K')
            T_evp = self.evp_air_in.Tdb - dT
            return T_evp.to('degC').m
        else:
            return T_evp_ini.to('degC').m

    def _guess_initial_T_cnd(self, T_cnd_ini: Quantity | None) -> float:
        """
        If the user does not provide an initial guess for the condensation
        temperature to method `rate`, an initial value of 20 K above the
        temperature of the air entering the condenser is taken.
        """
        if T_cnd_ini is None:
            dT = Q_(20, 'K')
            T_cnd = self.cnd_air_in.Tdb + dT
            return T_cnd.to('degC').m
        else:
            return T_cnd_ini.to('degC').m

    def _get_physical_bounds(self) -> tuple[tuple[float, float], ...]:
        """
        Returns the lower and upper limit of the evaporation and condensation
        temperature for use with the minimization algorithm in method `rate`.
        The lower limit of the evaporation temperature is set to negative
        infinity. Its upper limit is set to the temperature of the air entering
        the evaporator.
        The lower limit of the condensation temperature is set to the
        temperature of the air entering the condenser. Its upper limit is set
        1 K below the critical point temperature of the refrigerant.
        """
        T_crit = self.refrigerant.critical_point.T
        T_cnd_min = self.cnd_air_in.Tdb.to('degC').m
        T_cnd_max = (T_crit - Q_(1, 'K')).to('degC').m
        T_evp_min = -np.inf  # °C
        T_evp_max = self.evp_air_in.Tdb.to('degC').m
        return (
            (T_evp_min, T_evp_max),
            (T_cnd_min, T_cnd_max)
        )

    def _check_mass_balance(self) -> tuple[Quantity, Quantity]:
        """
        Returns the absolute error and relative error between the mass flow
        rate of refrigerant according to the compressor model (computed by the
        correlation in the compressor model with the resulting evaporation and
        condensation temperature) and the mass flow rate of refrigerant
        according to the evaporator model (computed by the evaporator model in
        order to maintain the degree of superheat set on the expansion device).
        The relative error is the ratio of the absolute error to the mass flow
        rate of refrigerant according to the compressor model. It expresses
        the deviation of the refrigerant mass flow rate according to the
        evaporator model with respect to the refrigerant mass flow rate according
        to the compressor model as a fraction (percentage).
        """
        abs_err = abs(self.evaporator.rfg_m_dot - self.compressor.m_dot)
        rel_err = abs_err / abs(self.compressor.m_dot)
        return abs_err, rel_err.to('pct')

    def _check_energy_balance(self) -> tuple[Quantity, Quantity]:
        """
        Returns the absolute error and relative error between the
        refrigeration capacity (heat absorption rate) according to the
        condensing unit model (computed by taking the difference between the
        heat rejection rate of the condenser and the compressor power) and
        the refrigeration capacity according to the evaporator model.
        The relative error is the ratio of the absolute error to the
        refrigeration capacity according to the compressor model. It expresses
        the deviation of the refrigeration capacity according to the evaporator
        model with respect to the refrigeration capacity according to the
        compressor model as a fraction (percentage).
        """
        Q_evp_cmp_model = self.condenser.Q_dot - self.compressor.W_dot
        Q_evp_evp_model = self.evaporator.Q_dot
        abs_err = abs(Q_evp_evp_model - Q_evp_cmp_model)
        rel_err = abs_err / abs(Q_evp_cmp_model)
        return abs_err, rel_err.to('pct')

    def _get_deviation(
        self,
        cmp_rfg_m_dot: Quantity,
        i: int,
        logger_on: bool = True
    ) -> float:
        """
        Calculates the mass flow rate of refrigerant let through by the 
        expansion device and returns the deviation with the mass flow rate of
        refrigerant displaced by the compressor.

        Returns
        -------
        The value (float) of the deviation expressed in units kg/h.
        """
        if logger_on:
            logger.info(
                f"Iteration {i + 1}: "
                "Refrigerant mass flow rate from compressor = "
                f"{cmp_rfg_m_dot:~P.3f}"
            )

        if not self._refrigeration_lines:
            cnd_rfg_in = self.compressor.discharge_gas
        else:
            self.discharge_line.m_dot = cmp_rfg_m_dot
            self.discharge_line.rfg_state = self.compressor.discharge_gas
            dP_dis = self.discharge_line.pressure_drop
            P_dis = self.compressor.discharge_gas.P
            T_dis = self.compressor.discharge_gas.T
            cnd_rfg_in = self.refrigerant(T=T_dis, P=P_dis - dP_dis)

        if logger_on:
            logger.info(
                f"Iteration {i + 1}: "
                "Refrigerant entering condenser with "
                f"T = {cnd_rfg_in.T.to('degC'):~P.3f}, "
                f"P = {cnd_rfg_in.P.to('bar'):~P.3f}"
            )

        try:
            self.condenser.solve(
                air_in=self.cnd_air_in,
                air_m_dot=self.cnd_air_m_dot,
                rfg_in=cnd_rfg_in,
                rfg_m_dot=cmp_rfg_m_dot
            )
        except CondenserError as err:
            logger.error(
                f"Iteration {i + 1}: "
                f"{type(err).__name__}: {err}"
            )
            # raise err
            return float('inf')
            # --> keep the minimize-function running and tell it that the
            # deviation is infinity for the current combination of evaporation
            # and condensing temperature.

        if logger_on:
            logger.info(
                f"Iteration {i + 1}: "
                "Refrigerant leaving condenser with "
                f"T = {self.condenser.rfg_out.T.to('degC'):~P.3f}, "
                f"P = {self.condenser.rfg_out.P.to('bar'):~P.3f}"
            )

        if not self._refrigeration_lines:
            evp_rfg_in = self.refrigerant(
                P=self.compressor.P_evp,
                h=self.condenser.rfg_out.h,
            )
        else:
            self.liquid_line.m_dot = cmp_rfg_m_dot
            self.liquid_line.rfg_state = self.condenser.rfg_out
            dP_liq = self.liquid_line.pressure_drop
            P_liq = self.condenser.rfg_out.P
            T_liq = self.condenser.rfg_out.T
            h_liq = self.refrigerant(P=P_liq - dP_liq, T=T_liq).h
            self.suction_line.m_dot = cmp_rfg_m_dot
            self.suction_line.rfg_state = self.compressor.suction_gas
            dP_suc = self.suction_line.pressure_drop
            P_suc = self.compressor.suction_gas.P
            evp_rfg_in = self.refrigerant(P=P_suc + dP_suc, h=h_liq)

        if logger_on:
            logger.info(
                f"Iteration {i + 1}: "
                "Refrigerant entering evaporator with "
                f"T = {evp_rfg_in.T.to('degC'):~P.3f}, "
                f"P = {evp_rfg_in.P.to('bar'):~P.3f}, "
                f"x = {evp_rfg_in.x.to('pct'):~P.3f}"
            )

        try:
            self.evaporator.solve(
                air_in=self.evp_air_in,
                air_m_dot=self.evp_air_m_dot,
                rfg_in=evp_rfg_in,
                dT_sh=self.dT_sh,
            )
        except EvaporatorError as err:
            logger.error(
                f"Iteration {i + 1}: "
                f"{type(err).__name__}: {err}"
            )
            # raise err
            return float('inf')  # --> keep the minimize-function running

        if logger_on:
            logger.info(
                f"Iteration {i + 1}: "
                "Refrigerant leaving evaporator with "
                f"T = {self.evaporator.rfg_out.T.to('degC'):~P.3f}, "
                f"P = {self.evaporator.rfg_out.P.to('bar'):~P.3f}"
            )

        evp_rfg_m_dot = self.evaporator.rfg_m_dot.to('kg / hr')

        if logger_on:
            logger.info(
                f"Iteration {i + 1}: "
                "Refrigerant mass flow rate through evaporator = "
                f"{evp_rfg_m_dot.to('kg / hr'):~P.3f}"
            )

        dev = evp_rfg_m_dot - cmp_rfg_m_dot

        if logger_on:
            logger.info(
                f"Iteration {i + 1}: "
                "Deviation between refrigerant mass flow rates = "
                f"{dev:~P.3f}"
            )

        return dev.magnitude

    def test_rate(
        self,
        evp_air_in: HumidAir,
        evp_air_m_dot: Quantity,
        cnd_air_in: HumidAir,
        cnd_air_m_dot: Quantity,
        n_cmp: Quantity | None = None,
        T_evp: Quantity | None = None,
        T_cnd: Quantity | None = None
    ) -> float:
        """
        For the given operating conditions, returns the absolute value of the
        deviation between the mass flow rate let through by the expansion device
        (to maintain the set degree of refrigerant superheating at the
        evaporator outlet) and the mass flow rate of refrigerant displaced by
        the compressor.
        """
        self._init(
            evp_air_in, evp_air_m_dot,
            cnd_air_in, cnd_air_m_dot,
            n_cmp
        )
        counter = [0]
        T_evp = self._guess_initial_T_evp(T_evp)
        T_cnd = self._guess_initial_T_cnd(T_cnd)
        unknowns = np.array([T_evp, T_cnd])
        dev = self.__fun_rate__(unknowns, counter)
        return dev
