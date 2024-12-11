"""
DESIGN OF SILL SUPPLY

For a room with given dimensions, room air temperature, and room load, determine
the required supply air volume flow rate, supply air temperature, and size of
the supply opening.
"""
from hvac import Quantity
from hvac.air_diffusion.design.sill_supply import RoomInfo, design_sill_supply


Q_ = Quantity


room_info = RoomInfo(
    L=Q_(8, 'm'),
    B=Q_(4, 'm'),
    H=Q_(3.5, 'm'),
    Z=Q_(1.8, 'm'),
    T_r=Q_(24, 'degC'),
    Q_dot=Q_(1.280, 'kW')
)

ss = design_sill_supply(
    room_info,
    b=Q_(1.2, 'm'),
    H_sc=Q_(2.5, 'm'),
    dT_o=Q_(-10, 'K')  # initial guess
)
print(ss.output)
