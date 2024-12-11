"""
DESIGN OF SIDE-WALL AIR SUPPLY

For a room with given dimensions, room air temperature, and room load, determine
the required supply air volume flow rate, supply air temperature, and size and
aspect ratio of the supply opening, such that the critical room Archimedes
number and mean room air speed are met.
"""
from hvac import Quantity
from hvac.air_diffusion.design.side_wall_supply import RoomInfo, design_side_wall_supply


Q_ = Quantity

room_info = RoomInfo(
    L=Q_(8, 'm'),
    B=Q_(4, 'm'),
    H=Q_(3.5, 'm'),
    Z=Q_(1.8, 'm'),
    T_r=Q_(26, 'degC'),
    Q_dot=Q_(1.280, 'kW')
)

sws = design_side_wall_supply(
    room_info,
    h_o=Q_(400, 'mm'),
    d=Q_(200, 'mm'),
    K1=6.3,
    Ar_r_crit=5500
)
print(sws.output)
