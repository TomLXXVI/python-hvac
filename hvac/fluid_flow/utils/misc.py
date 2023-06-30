from hvac import Quantity
from hvac.fluids import Fluid, FluidState

Water = Fluid('Water')
STANDARD_PRESSURE = Quantity(101_325, 'Pa')


def design_volume_flow_rate(Q_design: Quantity, T_avg: Quantity, deltaT: Quantity) -> Quantity:
    """
    Calculates the required network volume flow rate of water to meet the network heating/cooling load.

    Parameters
    ----------
    Q_design: Quantity
        The network thermal power load.
    T_avg: Quantity
        Average water temperature (average of supply and return temperature).
    deltaT: Quantity
        Temperature difference between supply and return water.

    Returns
    -------
    Quantity
        The network volume flow rate.
    """
    w = Water(T=T_avg, P=STANDARD_PRESSURE)
    rho = w.rho.to('kg / m ** 3').magnitude
    c = w.cp.to('J / (kg * K)').magnitude
    deltaT = deltaT.to('K').magnitude
    Q = Q_design.to('W').magnitude
    V = Q / (rho * c * deltaT)
    return Quantity(V, 'm ** 3 / s')


def head_to_pressure(
    H: Quantity,
    fluid: FluidState = Water(T=Quantity(15.0, 'degC'), P=STANDARD_PRESSURE),
    g: Quantity = Quantity(9.81, 'm / s ** 2')
) -> Quantity:
    f = fluid
    rho = f.rho.to('kg / m ** 3').magnitude
    g = g.to('m / s ** 2').magnitude
    H = H.to('m').magnitude
    P = rho * g * H
    return Quantity(P, 'Pa')


def pressure_to_head(
    P: Quantity,
    fluid: FluidState = Water(T=Quantity(15.0, 'degC'), P=STANDARD_PRESSURE),
    g: Quantity = Quantity(9.81, 'm / s ** 2')
) -> Quantity:
    f = fluid
    rho = f.rho.to('kg / m ** 3').magnitude
    g = g.to('m / s ** 2').magnitude
    P = P.to('Pa').magnitude
    H = P / (rho * g)
    return Quantity(H, 'm')
