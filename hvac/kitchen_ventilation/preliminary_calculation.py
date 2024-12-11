"""Implementation of quick methods according to standard EN 16282-1 for the
preliminary calculation of the design extraction air volume flow rate in
commercial kitchens.
"""
from hvac import Quantity
from enum import Enum

Q_ = Quantity


class KitchenActivity(Enum):
    """Relates the kind of kitchen activity to the required extraction airflow
    per unit of kitchen floor area.
    """
    OVERALL_KITCHEN = Q_(90, 'm ** 3 / (h * m ** 2)')
    ROASTING_GRILLING_AND_BAKING = Q_(120, 'm ** 3 / (h * m ** 2)')
    DISHWASHING_AREA = Q_(120, 'm ** 3 / (h * m ** 2)')


def hourly_air_change_rate_method(
    floor_area: Quantity,
    kitchen_type: KitchenActivity
) -> Quantity:
    """Returns a rough estimation for the extraction airflow in a kitchen
    depending on the floor area of the kitchen.

    Parameters
    ----------
    floor_area:
        Floor area of the kitchen.
    kitchen_type:
        Specifies the kind of kitchen activity in the considered area. See
        enum `KitchenType`.

    Notes
    -----
    To prevent spreading of kitchen odours to areas outside the kitchen, a
    slight negative pressure in the kitchen is recommended. However, the
    extraction airflow shall not exceed the supply airflow by more than 10 %.
    """
    A = floor_area.to('m ** 2')
    V_dot_ext = kitchen_type.value * A
    return V_dot_ext


class ApplianceDutyCategory(Enum):
    """Relates the duty category of appliances to an airflow velocity for
    capturing and removing the particulate emission.

    LIGHT:
        Steaming ovens, boiling pans, bains maries, and stock-pot stoves.
    MEDIUM:
        Deep fat fryers, bratt pans, solid and open top ranges and griddles.
    HEAVY:
        Chargrills, mesquite and specialist broiler units.
    """
    LIGHT = Q_(0.15, 'm / s')
    MEDIUM = Q_(0.225, 'm / s')
    HEAVY = Q_(0.3, 'm / s')


def hood_capture_velocity_method(
    hood_perimeter: Quantity,
    free_height: Quantity,
    appliance_duty_category: ApplianceDutyCategory
) -> Quantity:
    """Returns the required extraction airflow for the hood based on the size
    of the hood and the load category of the appliance, which determines the
    necessary capture velocity.

    Parameters
    ----------
    hood_perimeter:
        Unobstructed perimeter of the hood.
    free_height:
        Free height between the lower edge of the hood and the top edge of the
        cooking appliance. In case of an appliance with a vertical opening, e.g.
        salamander or oven, the free height is measured from the middle of the
        opening height.
    appliance_duty_category:
        The duty category of the appliance: light, medium or heavy. See
        enum `ApplianceDutyCategory`.
    """
    U = hood_perimeter.to('m')
    h_d = free_height.to('m')
    seconds_to_hours = Q_(3600, 's / h')
    V_dot_ext = appliance_duty_category.value * seconds_to_hours * U * h_d
    return V_dot_ext
