from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import HumidAir, Fluid
from hvac.heat_exchanger.fintube.continuous_fin import PlainFinTubeCounterFlowAirCondenser
from hvac.heat_exchanger.fintube.continuous_fin.air_condenser import logger

logger.setLevel(ModuleLogger.DEBUG)

Q_ = Quantity
R134a = Fluid('R134a')


condenser = PlainFinTubeCounterFlowAirCondenser(
    W_fro=Q_(1.003, 'm'),
    H_fro=Q_(0.400, 'm'),
    N_rows=5,
    S_trv=Q_(25.4, 'mm'),
    S_lon=Q_(22.0, 'mm'),
    D_int=Q_(8.422, 'mm'),
    D_ext=Q_(10.2, 'mm'),
    t_fin=Q_(0.3302, 'mm'),
    N_fin=1 / Q_(3.175, 'mm')
)


air_out, rfg_out = condenser.solve(
    air_in=HumidAir(Tdb=Q_(35.0, 'degC'), RH=Q_(30, 'pct')),
    air_m_dot=Q_(3216.315, 'kg / hr'),
    rfg_in=R134a(T=Q_(84.585, 'degC'), P=Q_(12.669, 'bar')),
    rfg_m_dot=Q_(171.731, 'kg / hr')
)

air_Q_dot = condenser.air_m_dot * (condenser.air_out.h - condenser.air_in.h)
rfg_Q_dot = condenser.rfg_m_dot * (condenser.rfg_in.h - condenser.rfg_out.h)
print(air_Q_dot.to('kW'), rfg_Q_dot.to('kW'))
