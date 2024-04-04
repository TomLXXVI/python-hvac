"""
ESTIMATION OF ENERGY USE IN HEATING OPERATION (part 1):
Introducing the `HeatingLoad` class and the `HeatPump` class.
"""
from hvac import Quantity
from hvac.energy_estimation.heat_pump import HeatPump
from hvac.energy_estimation.load import HeatingLoad

Q_ = Quantity

# ------------------------------------------------------------------------------
# HEATING LOAD
# ------------------------------------------------------------------------------
# Define the heating load of a building:
load = HeatingLoad(
    T_int=Q_(20, 'degC'),
    T_ext_min=Q_(-8, 'degC'),
    H_trm=Q_(4 / 28, 'kW / K'),
    V_dot_ven=Q_(100, 'm**3 / hr'),
    Q_dot_ihg=Q_(500, 'W'),
    eta_sys=Q_(80, 'pct'),
    time_segment=None
)

# Set the outdoor temperature and the number of hours this temperature is
# present:
load.T_ext = Q_(5, 'degC')
load.num_hours = 10


print(
    "fixed heating load characteristics:"
)
print(
    f"ventilation heat loss coefficient = {load.H_ven.to('kW / K'):~P.3f}",
    f"total heat loss coefficient = {load.H_tot.to('kW / K'):~P.3f}",
    f"building's balance temperature = {load.T_bal.to('degC'):~P.1f}",
    sep='\n', end='\n\n'
)
print(
    "heating load characteristics that depend on the set outdoor temperature:"
)
print(
    f"current outdoor temperature = {load.T_ext.to('degC'):~P.1f}",
    f"heating demanded by load = {load.Q_dot_out.to('kW'):~P.3f}",
    f"heating delivered by system = {load.Q_dot_in.to('kW'):~P.3f}",
    f"transmission heat loss = {load.Q_dot_trm.to('kW'):~P.3f}",
    f"ventilation heat loss = {load.Q_dot_ven.to('kW'):~P.3f}",
    f"delivered heating energy = {load.Q_in.to('kWh'):~P.3f}",
    sep='\n', end='\n\n'
)

# ------------------------------------------------------------------------------
# HEAT PUMP
# ------------------------------------------------------------------------------

# Define the heat pump used to heat the building, passing the heat pump
# performance table which tabulates the heating capacity and input power of the
# heat pump at different outdoor air temperatures, while the condenser-side
# temperature remains fixed (in this specific example water at 35 °C):
heat_pump = HeatPump(file_path="hp_atlantic_alfea_extensa_AI5_W35.csv")

# Attach the building's heating load to the heat pump:
heat_pump.load = load

# Set the outdoor temperature and the number of hours this temperature is
# present:
heat_pump.T_ext = Q_(5, 'degC')
heat_pump.num_hours = 10

# Get the balance point of the heat pump for the given load:
T_ext_bal, Q_dot_bal = heat_pump.balance_point


print(
    "heat pump characteristics that depend on the set outdoor temperature "
    "only:"
)
print(
    f"heating capacity = {heat_pump.Q_dot.to('kW'):~P.3f}",
    f"input power = {heat_pump.W_dot.to('kW'):~P.3f}",
    f"COP = {heat_pump.COP:~P.2f}",
    sep='\n', end='\n\n'
)
print(
    "heat pump characteristics that also depend on the set load:"
)
print(
    f"balance point temperature = {T_ext_bal.to('degC'):~P.1f}",
    f"heating capacity available at balance point = {Q_dot_bal.to('kW'):~P.1f}",
    f"PLR = {heat_pump.PLR:.3f}",
    f"PLF = {heat_pump.PLF:.3f}",
    f"COP = {heat_pump.COP_pl:~P.2f}",
    f"heat output = {heat_pump.Q_dot_pl.to('kW'):~P.3f}",
    f"input power = {heat_pump.W_dot_pl.to('kW'):~P.3f}",
    f"auxiliary heat = {heat_pump.Q_dot_aux.to('kW'):~P.3f}",
    f"heat pump energy consumption = {heat_pump.W.to('kWh'):~P.3f}",
    f"auxiliary heat energy = {heat_pump.Q_aux.to('kWh'):~P.3f}",
    sep='\n', end='\n\n'
)
