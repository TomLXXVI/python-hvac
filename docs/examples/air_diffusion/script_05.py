"""
DESIGN OF CEILING SUPPLY WITH LINEAR DIFFUSER

For a room with given dimensions, room air temperature, and room load, determine
the required supply air volume flow rate, supply air temperature, and size of
the supply opening.
"""
from hvac import Quantity
from hvac.air_diffusion.design.ceiling_supply import RoomInfo, design_linear_ceiling_supply

Q_ = Quantity

room_info = RoomInfo(
    L=Q_(4, 'm'),
    B=Q_(4, 'm'),
    H=Q_(3, 'm'),
    Z=Q_(1.8, 'm'),
    T_r=Q_(24, 'degC'),
    Q_dot=Q_(2400, 'W')
)

lcs = design_linear_ceiling_supply(
    room_info,
    dT_o=Q_(-15, 'K'),
    position='end'
)
print(lcs.output)
