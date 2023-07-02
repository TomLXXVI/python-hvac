# import numpy as np
# from scipy.optimize import fsolve
from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import HumidAir, Fluid, FluidState, CP_HUMID_AIR
from hvac.heat_transfer.heat_exchanger.fin_tube import air_evaporator
from hvac.heat_transfer.heat_exchanger.fin_tube import air_condenser
from hvac.vapor_compression.real_compressor import FixedSpeedCompressor, VariableSpeedCompressor

Q_ = Quantity

logger = ModuleLogger.get_logger(__name__)

Evaporator = air_evaporator.rating.PFT_CO_EVP
Condenser = air_condenser.rating.PFT_CO_CND
Compressor = FixedSpeedCompressor | VariableSpeedCompressor


class SingleStageVaporCompressionMachine:

    def __init__(
        self,
        evaporator: Evaporator,
        condenser: Condenser,
        compressor: Compressor,
        Refrigerant: Fluid
    ) -> None:
        self.evaporator = evaporator
        self.condenser = condenser
        self.compressor = compressor
        self.Refrigerant = Refrigerant
        self.evp_m_dot_air: Quantity | None = None
        self.cnd_m_dot_air: Quantity | None = None
        self.evp_air_in: HumidAir | None = None
        self.cnd_air_in: HumidAir | None = None
        self.dT_super: Quantity | None = None
        self.dT_sub: Quantity | None = None
        self.evp_rfg_sat_vap: FluidState | None = None
        self.cnd_rfg_sat_liq: FluidState | None = None

    def set_operating_conditions(
        self,
        evp_m_dot_air: Quantity,
        cnd_m_dot_air: Quantity,
        evp_air_in: HumidAir,
        cnd_air_in: HumidAir,
        dT_super: Quantity,
        dT_sub: Quantity,
        n_cmp: Quantity | None = None
    ) -> None:
        self.evp_m_dot_air = evp_m_dot_air
        self.cnd_m_dot_air = cnd_m_dot_air
        self.evp_air_in = evp_air_in
        self.cnd_air_in = cnd_air_in
        self.dT_super = dT_super.to('K')
        self.dT_sub = dT_sub.to('K')
        if isinstance(self.compressor, VariableSpeedCompressor):
            self.compressor.speed = n_cmp

    def _calculate_T_evp(self, evp_Q_dot: Quantity) -> Quantity:
        r = self.evaporator(
            air_in=self.evp_air_in,
            m_dot_air=self.evp_m_dot_air,
            rfg_in=self.compressor.mixture,
            dT_rfg_sh=self.dT_super,
            m_dot_rfg_ini=self.compressor.m_dot
        )
        air_sat_out = HumidAir(
            h=self.evp_air_in.h - evp_Q_dot / (r.eps * self.evp_m_dot_air),
            RH=Q_(100, 'pct')
        )
        return air_sat_out.Tdb

    def _calculate_T_cnd(self, cnd_Q_dot: Quantity) -> Quantity:
        r = self.condenser(
            m_dot_air=self.cnd_m_dot_air,
            m_dot_rfg=self.compressor.m_dot,
            air_in=self.cnd_air_in,
            rfg_in=self.compressor.discharge_gas
        )
        C_air = self.cnd_m_dot_air * CP_HUMID_AIR
        T_cnd = self.cnd_air_in.Tdb + cnd_Q_dot / (r.eps * C_air)
        return T_cnd

    def _find_steady_operating_point(self, evp_T_ini: Quantity, cnd_T_ini: Quantity) -> None:
        """Find evaporator and condenser temperature under the given operating
        conditions."""
        T_evp, T_cnd = evp_T_ini, cnd_T_ini
        i = 0
        i_max = 20
        tol_T = Q_(0.1, 'K')
        while i < i_max:
            logger.debug(
                f"iteration {i + 1}: "
                f"Try T_evp = {T_evp.to('degC'):~P.3f}; "
                f"T_cnd = {T_cnd.to('degC'):~P.3f}"
            )
            self.compressor.Te = T_evp
            self.compressor.Tc = T_cnd
            evp_Q_dot = self.compressor.Qc_dot
            cnd_Q_dot = self.compressor.Qh_dot
            T_evp_new = self._calculate_T_evp(evp_Q_dot)
            T_cnd_new = self._calculate_T_cnd(cnd_Q_dot)
            dev_T_evp = abs(T_evp_new.to('K') - T_evp)
            dev_T_cnd = abs(T_cnd_new.to('K') - T_cnd)
            if (dev_T_evp <= tol_T) and (dev_T_cnd <= tol_T):
                self.compressor.Te = T_evp_new
                self.compressor.Tc = T_cnd_new
                break
            T_evp = T_evp_new
            T_cnd = T_cnd_new
            i += 1
        else:
            raise ValueError(
                f"no steady-state operating point found after {i_max} iterations"
            )

    def simulate(self) -> None:
        """Run the calculations to find evaporator temperature and condenser
         temperature at the current set working conditions. Call this method after
         first setting `lt_air_in`, `ht_air_in` and/or `rpm`."""
        if self.evp_air_in is not None:
            evp_T_ini = self.evp_air_in.Tdb - Q_(20, 'K')
        else:
            evp_T_ini = None
        if self.cnd_air_in is not None:
            cnd_T_ini = self.cnd_air_in.Tdb + Q_(20, 'K')
        else:
            cnd_T_ini = None
        if evp_T_ini is not None and cnd_T_ini is not None:
            self._find_steady_operating_point(evp_T_ini, cnd_T_ini)

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
