"""Simulation of a heat pump (single-stage vapor compression machine) in
cooling operation. Aim of the simulation is to find the steady-state operating
point of the heat pump under given operating conditions.
"""
import time
from pathlib import Path
from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import Fluid, HumidAir
from hvac.vapor_compression import VariableSpeedCompressor
from hvac.vapor_compression.machine_bis import SingleStageVaporCompressionMachine
from hvac.heat_transfer.heat_exchanger.fin_tube import air_evaporator, air_condenser

logger = ModuleLogger.get_logger('hvac.fluids.fluid')
logger.setLevel(ModuleLogger.ERROR)

Q_ = Quantity
R134a = Fluid('R134a')
Evaporator = air_evaporator.rating.PFT_CO_EVP
Condenser = air_condenser.rating.PFT_CO_CND

compressor = VariableSpeedCompressor(
    coeff_file=Path("VTZ038-G_R134a.csv"),
    refrigerant_type=R134a,
    units={'m_dot': 'kg / hr', 'speed': '1 / s'}
)

evaporator = Evaporator(
    L1=Q_(0.844, 'm'),
    L3=Q_(0.211, 'm'),
    N_r=6,
    S_t=Q_(25.4, 'mm'),
    S_l=Q_(22.0, 'mm'),
    D_i=Q_(8.422, 'mm'),
    D_o=Q_(10.2, 'mm'),
    t_f=Q_(0.3302, 'mm'),
    N_f=1 / Q_(3.175, 'mm')
)

condenser = Condenser(
    L1=Q_(0.550, 'm'),
    L3=Q_(0.275, 'm'),
    N_r=10,
    S_t=Q_(25.4, 'mm'),
    S_l=Q_(22.0, 'mm'),
    D_i=Q_(8.422, 'mm'),
    D_o=Q_(10.2, 'mm'),
    t_f=Q_(0.3302, 'mm'),
    N_f=1 / Q_(3.175, 'mm')
)

machine = SingleStageVaporCompressionMachine(
    evaporator,
    condenser,
    compressor,
    R134a
)
machine.set_operating_conditions(
    evp_m_dot_air=Q_(1500, 'kg / hr'),
    cnd_m_dot_air=Q_(2454.820, 'kg / hr'),
    evp_air_in=HumidAir(Tdb=Q_(24, 'degC'), RH=Q_(50, 'pct')),
    cnd_air_in=HumidAir(Tdb=Q_(35.0, 'degC'), RH=Q_(30, 'pct')),
    dT_super=Q_(10, 'K'),
    n_cmp=Q_(3932, '1 / min')
)

t_start = time.time()
machine.simulate()
t_end = time.time()
print(f"execution time: {t_end - t_start}")

print(
    f"heat absorption rate = {machine.Qc_dot.to('kW'):~P.3f}\n"
    f"heat rejection rate = {machine.Qh_dot.to('kW'):~P.3f}\n"
    f"compressor power = {machine.Wc_dot.to('kW'):~P.3f}\n"
    f"COP = {machine.COP.to('frac'):~P.2f}\n"
    f"refrigerant mass flow rate = {machine.m_dot.to('kg / hr'):~P.3f}\n"
    f"evaporation temperature = {machine.Te.to('degC'):~P.2f}\n"
    f"condensing temperature = {machine.Tc.to('degC'):~P.2f}\n"
)

print(
    f"suction gas temperature = {machine.suction_gas.T.to('degC'):~P.2f}\n"
    f"discharge gas temperature = {machine.discharge_gas.T.to('degC'):~P.2f}\n"
    f"liquid temperature = {machine.liquid.T.to('degC'):~P.2f}\n"
    f"liquid/vapor mixture temperature = {machine.mixture.T.to('degC'):~P.2f}\n"
    f"subcooling degree = {machine.sub_cooling.to('K'):~P.2f}\n"
)

print(
    f"evaporator air out = {evaporator.air_out.Tdb.to('degC'):~P.2f} DB, "
    f"{evaporator.air_out.RH.to('pct'):~P.1f} RH\n"
    f"condenser air out = {condenser.air_out.Tdb.to('degC'):~P.2f} DB, "
    f"{condenser.air_out.RH.to('pct'):~P.1f} RH\n"
)

print(
    f"evaporator air-side pressure drop = {evaporator.dP_air.to('Pa'):~P.1f}\n"
    f"condenser air-side pressure drop = {condenser.dP_air.to('Pa'):~P.1f}"
)
