"""
DESIGN OF CEILING SUPPLY WITH CIRCULAR DIFFUSER

For a room with given dimensions, room air temperature, and room load, determine
the required supply air volume flow rate, supply air temperature, and size of
the supply opening.
"""
from hvac import Quantity
from hvac.air_diffusion.design.ceiling_supply import (
    RoomInfo,
    design_circular_ceiling_supply,
    VAVSupply
)

Q_ = Quantity

room_info = RoomInfo(
    L=Q_(4, 'm'),
    B=Q_(4, 'm'),
    H=Q_(3, 'm'),
    Z=Q_(1.8, 'm'),
    T_r=Q_(24, 'degC'),
    Q_dot=Q_(800, 'W')
)

ccs = design_circular_ceiling_supply(
    room_info,
    dT_o=Q_(-5, 'K'),
    d_o=Q_(0.6, 'm'),
)
print(ccs.output)


vav_supply = VAVSupply(
    room_info,
    A_o=ccs.output.A_o,
    dT_o=ccs.output.dT_o,
)
print(
    f"minimum room load that can be handled "
    f"with d_To {vav_supply.dT_o.to('K'):~P}: "
    f"{vav_supply.Q_dot_min.to('W'):~P.0f})",
    f"mean room air velocity at design load: {vav_supply.v_r.to('m / s'):~P.2f}",
    f"maximum local velocity at design load: {vav_supply.v_m.to('m / s'):~P.2f}",
    sep='\n'
)

v_m = Q_(0.35, 'm / s')
v_r_min = Q_(0.1, 'm / s')
v_r_max = Q_(0.25, 'm / s')

Q_dot_v_m, Q_dot_v_r_min, Q_dot_v_r_max = vav_supply.get_room_load_limits(v_m, v_r_min, v_r_max)

print(
    f"Q_dot that corresponds with v_m {v_m:~P}: {Q_dot_v_m.to('W'):~P.1f}",
    f"Q_dot that corresponds with v_r_min {v_r_min:~P}: {Q_dot_v_r_min.to('W'):~P.1f}",
    f"Q_dot that corresponds with v_r_max {v_r_max:~P}: {Q_dot_v_r_max.to('W'):~P.1f}",
    sep='\n'
)
