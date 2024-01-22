from hvac import Quantity
from hvac.fluids import HumidAir, Fluid
from hvac.heat_exchanger.fintube.continuous_fin import PlainFinTubeCounterFlowAirEvaporator


Q_ = Quantity
R134a = Fluid('R134a')


evaporator = PlainFinTubeCounterFlowAirEvaporator(
    W_fro=Q_(0.731, 'm'),
    H_fro=Q_(0.244, 'm'),
    N_rows=3,
    S_trv=Q_(25.4, 'mm'),
    S_lon=Q_(22.0, 'mm'),
    D_int=Q_(8.422, 'mm'),
    D_ext=Q_(10.2, 'mm'),
    t_fin=Q_(0.3302, 'mm'),
    N_fin=1 / Q_(3.175, 'mm')
)


m_dot_rfg = evaporator.solve(
    air_in=HumidAir(Tdb=Q_(24.0, 'degC'), RH=Q_(50, 'pct')),
    air_m_dot=Q_(1500, 'kg / hr'),
    rfg_in=R134a(T=Q_(5.0, 'degC'), x=Q_(22.788, 'pct')),
    dT_sh=Q_(10, 'K')
)


air_Q_dot = evaporator.air_m_dot * (evaporator.air_in.h - evaporator.air_out.h)
rfg_Q_dot = evaporator.rfg_m_dot * (evaporator.rfg_out.h - evaporator.rfg_in.h)
print(air_Q_dot.to('kW'), rfg_Q_dot.to('kW'))
print(evaporator.P_evp.to('bar'))
