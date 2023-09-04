"""Example of a cooling coil load calculation in a single-zone air cooling
system.

Based on the cooling load calculation of a single-zone building, the following
elements need to be determined:
- the required capacity of the cooling coil
- the sensible heat ratio of the cooling coil
- the state of air entering and leaving the cooling coil
- the required mass flow rate of air through the cooling coil
- the required volume flow rate of air through the cooling coil referred to
NTP conditions for selection of the air supply fan
"""
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning.single_zone.cooling_design import AircoSystem


Q_ = Quantity


# DESIGN CONDITIONS
# =================

# Ventilation requirements
# ------------------------
# Ventilation air requirement per person:
V_p_vent = Q_(22, 'm ** 3 / hr')
# Design occupation of the zone (number of persons in the zone):
n_p = 200
# Required ventilation air volume flow rate for the zone:
V_vent = n_p * V_p_vent


# Cooling loads
# -------------
# Sensible cooling load of the zone at summer peak design conditions:
Q_sen_zone = Q_(30.514, 'kW')
# Latent cooling load of the zone:
Q_lat_zone = Q_(12.967, 'kW')
# Total cooling load:
Q_zone = Q_sen_zone + Q_lat_zone
# Sensible heat ratio of the zone:
SHR_zone = Q_sen_zone / Q_zone

# State of zone air at summer peak design conditions
# --------------------------------------------------
zone_air = HumidAir(Tdb=Q_(26, 'degC'), W=Q_(8.0, 'g / kg'))

# State of outdoor air under summer peak design conditions
# --------------------------------------------------------
outdoor_air = HumidAir(Tdb=Q_(31, 'degC'), W=Q_(10.0, 'g / kg'))


# AIR-CONDITIONING SYSTEM SETUP
# =============================
airco_system = AircoSystem(
    zone_air=zone_air,
    outdoor_air=outdoor_air,
    Q_dot_zone=Q_zone,
    SHR_zone=SHR_zone,
    V_dot_vent_ntp=V_vent,
    T_supply=Q_(16, 'degC'),  # designer's choice
    eps_hr_h=Q_(22, 'pct'),   # enthalpy-effectiveness of heat recovery wheel
    eps_hr_W=Q_(11, 'pct')    # humidity-effectiveness of heat recovery wheel
)

output = airco_system.design()

# RESULTS
# =======
# Display the results needed for sizing the air cooling coil:
print(output)

# Show the air conditioning processes on the psychrometric chart:
psy_chart = airco_system.psychrometric_chart
psy_chart.show()
