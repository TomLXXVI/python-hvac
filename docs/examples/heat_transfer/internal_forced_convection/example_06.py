"""
EXAMPLE 6
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 9.7
A pump is used to fill an open tank with water. The inlet to the pump pulls
water directly from a reservoir at T_o = 70 °F and p_o = 1 atm. The pump
discharge is connected to a length of new cast iron pipe with inner diameter
D = 1.5 inch and length L = 500 ft. The exit of the pipe discharges into the
tank at approximately the same elevation as the reservoir. The volume of the
tank is V = 600 gal. The pump behavior is characterized by the pump curve in
[Figure 1], which shows the pressure rise by the pump as a function of the
flow rate that it provides. Determine the time required to fill the tank with
this pump.
"""
import warnings
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from hvac.heat_transfer.forced_convection.internal_flow import CircularTube

warnings.filterwarnings('ignore', category=CoolPropWarning)

Q_ = Quantity

Water = Fluid('Water')
water = Water(T=Q_(70, 'degF'), P=Q_(1, 'atm'))

# New cast iron pipe:
pipe = CircularTube(
    Di=Q_(1.5, 'inch'),
    L=Q_(500, 'ft'),
    fluid=water,
    e=Q_(0.5, 'mm')
)

# Pump:
V_dot_oc = Q_(100, 'gal / min')  # maximum flow rate when open circuited
dp_dh = Q_(25, 'psi')  # shut-off pump pressure


def pump_curve(V_dot: Quantity) -> Quantity:
    """Returns the pump's pressure rise for a given volume flow rate."""
    V_dot_frac = V_dot.to('gal / min') / V_dot_oc
    dp_frac = 1 - 0.0353 * V_dot_frac - 0.9647 * (V_dot_frac ** 2)
    dp_pump = dp_frac * dp_dh
    return dp_pump


def system_curve(V_dot: Quantity) -> Quantity:
    """Returns the pressure loss in the pipe for a given volume flow rate."""
    pipe.V_dot = V_dot
    dp_loss = pipe.pressure_drop()
    return dp_loss.to('psi')


def working_point() -> tuple[Quantity, Quantity]:
    """Returns the working point (volume flow rate, pressure rise) on the
    pump's curve.
    """
    # Numerical pump equation:
    def _pump_curve(V_dot: float) -> float:
        V_dot = Q_(V_dot, 'gal / min')
        dp_pump = pump_curve(V_dot)
        return dp_pump.to('psi').m

    # Numerical system equation:
    def _system_curve(V_dot: float) -> float:
        V_dot = Q_(V_dot, 'gal / min')
        dp_sys = system_curve(V_dot)
        return dp_sys.to('psi').m

    # The working point of the pump is where the pressure rise by the pump
    # equals the pressure loss in the system. We will use a root finding
    # algorithm from the third-party library Scipy to find the volume flow rate
    # for which the difference between the pump's pressure rise and the system's
    # pressure loss becomes zero.
    def _working_point(V_dot: float) -> float:
        dp_pump = _pump_curve(V_dot)
        dp_sys = _system_curve(V_dot)
        return dp_pump - dp_sys

    # Find this volume flow rate within the pump's volume flow rate range:
    V_dot = root_scalar(_working_point, bracket=[1.e-5, V_dot_oc.m]).root

    # Get the pump's pressure rise that corresponds with this volume flow rate:
    dp = _pump_curve(V_dot)
    return Q_(V_dot, 'gal / min'), Q_(dp, 'psi')


if __name__ == '__main__':
    # Determine the working point of the pump:
    V_dot_wp, dp_wp = working_point()
    print(
        f"V_dot pump = {V_dot_wp}",
        f"dP pump = {dp_wp}",
        sep='\n'
    )

    # Time required to fill the tank:
    V_tank = Q_(600, 'gal')
    print(f"fill time = {(V_tank / V_dot_wp).to('s')}")
