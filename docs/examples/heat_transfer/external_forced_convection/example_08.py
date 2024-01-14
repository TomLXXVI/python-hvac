"""
EXAMPLE 8
---------
From: Nellis G. F., & Klein S. A. (2021)
INTRODUCTION TO ENGINEERING HEAT TRANSFER.
Cambridge University Press.

Example 8.9
A copper sphere is falling through a bath of Therminol 66. The diameter of the
sphere is D = 20 mm and the therminol temperature is T_inf = 20 °C. Initially,
the sphere is falling at a rate of u = 0.5 m/s and has a temperature of
T = 200 °C. Determine the rate of change of the temperature and the rate of
change of the velocity of the sphere.
"""
import warnings
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from hvac.charts import LineChart
from hvac.heat_transfer.forced_convection.external_flow.sphere import Sphere

warnings.filterwarnings('ignore', category=CoolPropWarning)


Q_ = Quantity


Therminol66 = Fluid('T66', backend='INCOMP')
T_gly = Q_(20, 'degC').to('K')

c_Cu = Q_(0.385, 'J / (g * K)')
rho_Cu = Q_(8.96, 'g / cm ** 3')
sphere = Sphere(D=Q_(20, 'mm'), fluid=None)

g = Q_(9.81, 'm / s ** 2')  # gravity constant
dt = Q_(0.01, 's')  # time step


def fun(T_sph: Quantity, u_sph: Quantity) -> tuple[Quantity, Quantity]:
    """Returns the rate of change of the sphere's temperature and the
    rate of change of the sphere's velocity at a time moment t.

    Parameters
    ----------
    T_sph:
        Temperature of the sphere at time moment t.
    u_sph:
        Speed of the sphere at time moment t.
    """
    T_flm = (T_sph + T_gly) / 2
    t66 = Therminol66(T=T_flm, P=Q_(101.325, 'kPa'))
    sphere.fluid = t66
    sphere.u_inf = u_sph
    h_avg = sphere.avg_heat_trf_coeff()
    F_D = sphere.drag_force()
    Q = h_avg * sphere.area * (T_sph - T_gly)  # heat loss of sphere
    m_sph = rho_Cu * sphere.volume  # mass of sphere
    C_sph = m_sph * c_Cu  # thermal capacity of sphere
    der_T_sph = -Q / C_sph  # dT/dt sphere from heat balance
    der_u_sph = g - (F_D / m_sph)  # du/dt sphere from force balance
    return der_T_sph.to('K / s'), der_u_sph.to('m / s ** 2')


if __name__ == '__main__':
    # Initial conditions:
    T_sph_ini = Q_(200, 'degC').to('K')
    u_sph_ini = Q_(0.5, 'm / s')
    x_sph_ini = Q_(0, 'm')

    T_sph_arr = [T_sph_ini]
    u_sph_arr = [u_sph_ini]
    x_sph_arr = [x_sph_ini]
    for i in range(200):
        der_T_sph, der_u_sph = fun(T_sph_ini, u_sph_ini)
        dT_sph = der_T_sph * dt  # temperature change after time step dt
        du_sph = der_u_sph * dt  # speed change after time step dt
        T_sph = T_sph_ini + dT_sph  # temperature of the sphere at t + dt
        u_sph = u_sph_ini + du_sph  # speed of the sphere at t + dt
        dx_sph = (u_sph_ini + u_sph) / 2 * dt  # displacement of the sphere during dt
        x_sph = x_sph_ini + dx_sph
        T_sph_arr.append(T_sph)
        u_sph_arr.append(u_sph)
        x_sph_arr.append(x_sph)
        T_sph_ini = T_sph
        u_sph_ini = u_sph
        x_sph_ini = x_sph

    # Plot speed and temperature of the sphere as function of its position:
    chart = LineChart()
    chart.add_y2_axis()
    chart.add_xy_data(
        label='velocity',
        x1_values=[x_sph.to('m').m for x_sph in x_sph_arr],
        y1_values=[u_sph.to('m / s').m for u_sph in u_sph_arr]
    )
    chart.add_xy_data(
        label='temperature',
        x1_values=[x_sph.to('m').m for x_sph in x_sph_arr],
        y2_values=[T_sph.to('degC').m for T_sph in T_sph_arr],
        style_props={'color': 'orange'}
    )
    chart.x1.add_title('position, m')
    chart.y1.add_title('velocity, m/s')
    chart.y2.add_title('temperature, °C')
    chart.show()
