"""WINTER PEAK DESIGN CALCULATIONS OF A SINGLE-ZONE AIR SYSTEM

This script follows the example taken from the book "Heating and Cooling of
Buildings" by T. Agami Reddy, Jan F. Kreider, Peter S. Curtiss and Ari
Rabl (3rd Edition, 2017), Chapter 19: Example 19.4.
"""
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning.single_zone.heating_design import AircoSystem


Q_ = Quantity

# AIR-CONDITIONING SYSTEM SETUP
# =============================

airco_system = AircoSystem(
    # Zone air setpoint:
    zone_air=HumidAir(Tdb=Q_(72, 'degF'), RH=Q_(50, 'pct')),
    # Winter peak design condition of outdoor air:
    outdoor_air=HumidAir(Tdb=Q_(40, 'degF'), RH=Q_(40, 'pct')),
    # Heating zone load at winter peak design conditions (from the heating load
    # calculation of the building):
    Q_dot_zone=Q_(-150_000, 'Btu / hr'),
    SHR_zone=Q_(0.8, 'frac'),
    # Required minimum outdoor air volume flow rate for ventilation:
    V_dot_vent_ntp=Q_(1000, 'ft ** 3 / min'),
    # Required supply air mass flow rate to the zone at summer peak design
    # conditions (from the cooling load calculation of the building):
    m_dot_supply_cooling=Q_(17_140, 'lb / hr'),
    # Allowable maximum supply air temperature (to avoid stratification):
    T_supply_max=Q_(105, 'degF'),
    # Temperature of steam for supply air humidification:
    T_steam=Q_(200, 'degF'),
    # Set the units for displaying the output:
    units={
        'm_dot': ('lb / hr', 3),
        'V_dot': ('ft ** 3 / min', 3),
        'Q_dot': ('Btu / hr', 3),
        'T': ('degF', 3),
        'W': ('lb / lb', 6),
        'RH': ('pct', 3)
    }
)

# RESULTS
# =======
output = airco_system.design()
print(output)

psy_chart = airco_system.psychrometric_chart
psy_chart.show()
