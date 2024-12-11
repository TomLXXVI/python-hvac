"""Example of a quick method for preliminary calculation of the design
extraction air volume flow rate of a kitchen hood based on capture velocity
according to EN 16282-1.
"""
from hvac import Quantity
from hvac.kitchen_ventilation import (
    ApplianceDutyCategory,
    hood_capture_velocity_method
)

Q_ = Quantity

V_dot_ext = hood_capture_velocity_method(
    hood_perimeter=Q_(2 * (2.5 + 1.25), 'm'),
    free_height=Q_(1.2, 'm'),
    appliance_duty_category=ApplianceDutyCategory.LIGHT
)
print(V_dot_ext)
