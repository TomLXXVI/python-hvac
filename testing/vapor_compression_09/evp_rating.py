from pathlib import Path
from hvac import Quantity
from hvac.fluids import Fluid, HumidAir
from hvac.heat_exchanger.fintube.continuous_fin import (
    PlainFinTubeCounterFlowAirEvaporator,
    PlainFinTubeCounterFlowAirCondenser
)
from hvac.vapor_compression import (
    VariableSpeedCompressor,
    SingleStageVaporCompressionMachine
)

Q_ = Quantity
R410a = Fluid('R410a')


evp_air_in = HumidAir(Tdb=Q_(28.2, 'degC'), W=Q_(10.6, 'g / kg'))
evp_air_out = HumidAir(Tdb=Q_(14.0, 'degC'), W=Q_(9.1, 'g / kg'))
evp_air_m_dot = Q_(10_762.4, 'kg / hr')
evp_Q_dot = evp_air_m_dot * (evp_air_in.h - evp_air_out.h)
print(evp_Q_dot.to('kW'))

T_evp = Q_(6, 'degC')
evp_rfg_sv = R410a(T=T_evp, x=Q_(1, 'frac'))
P_evp = evp_rfg_sv.P

cnd_air_in = HumidAir(Tdb=Q_(32, 'degC'), Twb=Q_(21, 'degC'))
cnd_dT_air = Q_(20, 'K')

T_cnd = Q_(60, 'degC')
cnd_rfg_sl = R410a(T=T_cnd, x=Q_(0, 'frac'))
P_cnd = cnd_rfg_sl.P
dT_sc = Q_(10, 'K')
cnd_rfg_out = R410a(T=T_cnd.to('K') - dT_sc, P=P_cnd)

evp_rfg_in = R410a(h=cnd_rfg_out.h, P=P_evp)

dT_sh = Q_(10, 'K')
evp_rfg_out = R410a(T=T_evp.to('K') + dT_sh, P=P_evp)

rfg_m_dot = evp_Q_dot / (evp_rfg_out.h - evp_rfg_in.h)
print(rfg_m_dot.to('kg / hr'))

evp_v_fa = Q_(2, 'm / s')
evp_air_V_dot = evp_air_m_dot / evp_air_in.rho
evp_A_fa = evp_air_V_dot / evp_v_fa
evp_H_fa = Q_(1, 'm')
evp_W_fa = evp_A_fa / evp_H_fa

evaporator = PlainFinTubeCounterFlowAirEvaporator(
    W_fro=evp_W_fa,               # width of frontal area
    H_fro=evp_H_fa,               # height of frontal area
    N_rows=4,                     # number of rows
    S_trv=Q_(22.42, 'mm'),        # vertical distance between tubes
    S_lon=Q_(22.27, 'mm'),        # horizontal distance between tubes
    D_int=Q_(8.422, 'mm'),        # inner tube diameter
    D_ext=Q_(10.2, 'mm'),         # outer tube diameter
    t_fin=Q_(0.3302, 'mm'),       # fin thickness
    N_fin=1 / Q_(2.8, 'mm'),      # fin density
    k_fin=Q_(237, 'W / (m * K)')  # conductivity of fin material
)

rfg_m_dot = evaporator.solve(
    air_in=evp_air_in,
    air_m_dot=evp_air_m_dot,
    rfg_in=evp_rfg_in,
    dT_sh=dT_sh,
    rfg_m_dot_ini=rfg_m_dot
)
print(rfg_m_dot)
print(evaporator.air_out)
print(evaporator.Q_dot.to('kW'))

compressor = VariableSpeedCompressor(
    coeff_file=Path('./compressor_data/VZH170CGM_R410a.csv'),
    refrigerant_type=R410a,
    dT_sh=dT_sh,
    dT_sc=dT_sc,
    units={'m_dot': 'kg / hr', 'speed': '1 / s'}
)

compressor.Te = evaporator.T_evp
compressor.Tc = T_cnd
n_cmp = compressor.get_compressor_speed(m_dot=evaporator.rfg_m_dot)
print(n_cmp.to('1 / min'))
compressor.speed = n_cmp

cnd_Q_dot = evaporator.Q_dot + compressor.Wc_dot
cnd_air_out = HumidAir(
    Tdb=cnd_air_in.Tdb.to('K') + cnd_dT_air,
    W=cnd_air_in.W
)
cnd_air_m_dot = cnd_Q_dot / (cnd_air_out.h - cnd_air_in.h)
print(cnd_air_m_dot.to('kg / hr'))

cnd_v_fa = Q_(1.5, 'm / s')
cnd_air_V_dot = cnd_air_m_dot / cnd_air_in.rho
cnd_A_fa = cnd_air_V_dot / cnd_v_fa
cnd_H_fa = Q_(1, 'm')
cnd_W_fa = cnd_A_fa / cnd_H_fa

condenser = PlainFinTubeCounterFlowAirCondenser(
    W_fro=cnd_W_fa,               # width of frontal area
    H_fro=cnd_H_fa,               # height of frontal area
    N_rows=6,                     # number of rows
    S_trv=Q_(22.42, 'mm'),        # vertical distance between tubes
    S_lon=Q_(22.27, 'mm'),        # horizontal distance between tubes
    D_int=Q_(8.422, 'mm'),        # inner tube diameter
    D_ext=Q_(10.2, 'mm'),         # outer tube diameter
    t_fin=Q_(0.3302, 'mm'),       # fin thickness
    N_fin=1 / Q_(2.8, 'mm'),      # fin density
    k_fin=Q_(237, 'W / (m * K)')  # conductivity of fin material
)

cnd_air_out, cnd_rfg_out = condenser.solve(
    air_in=cnd_air_in,
    air_m_dot=cnd_air_m_dot,
    rfg_in=compressor.discharge_gas,
    rfg_m_dot=rfg_m_dot
)
print(cnd_air_out)
print(condenser.T_cnd.to('degC'))


vcm = SingleStageVaporCompressionMachine(
    evaporator, condenser, compressor, R410a,
    dT_sh=dT_sh,
    n_cmp_min=Q_(1500, '1 / min'),
    n_cmp_max=Q_(6000, '1 / min')
)


output = vcm.balance_by_speed(
    evp_air_in=evp_air_in,
    evp_air_m_dot=evp_air_m_dot,
    cnd_air_in=cnd_air_in,
    cnd_air_m_dot=cnd_air_m_dot,
    T_evp=T_evp,
    T_cnd=T_cnd
)

print(output.to_text())
