"""Cooling Coil Load Calculation

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
from hvac.air_conditioning.single_zone.design import AircoSystem


Q_ = Quantity


# DESIGN CONDITIONS

# Ventilation requirements:

# Ventilation air requirement per person:
V_p_vent = Q_(22, 'm ** 3 / hr')

# Design occupation of the zone (number of persons in the zone):
n_p = 200

# Required ventilation air volume flow rate for the zone:
V_vent = n_p * V_p_vent


# Cooling loads:

# Sensible cooling load of the zone at climatic design conditions:
Q_sen_zone = Q_(30.514, 'kW')

# Latent cooling load of the zone:
Q_lat_zone = Q_(12.967, 'kW')

# Total cooling load:
Q_zone = Q_sen_zone + Q_lat_zone

# Sensible heat ratio of the zone:
SHR_zone = Q_sen_zone / Q_zone

# State of zone air at design conditions:
zone_air = HumidAir(Tdb=Q_(26, 'degC'), W=Q_(8.0, 'g / kg'))

# State of outdoor air at climatic design conditions:
outdoor_air = HumidAir(Tdb=Q_(31, 'degC'), W=Q_(10.0, 'g / kg'))

# Set up air conditioning system:
airco_system = AircoSystem(
    zone_air=zone_air,
    outdoor_air=outdoor_air,
    Q_zone=Q_zone,
    SHR_zone=SHR_zone,
    V_vent=V_vent,
    T_supply=Q_(13, 'degC'),  # designer's choice
    eps_hr_h=Q_(22, 'pct'),
    eps_hr_W=Q_(11, 'pct')
)

# Get the results needed for sizing the air cooling coil:
print(
    f"cooling coil load @ design = "
    f"{airco_system.Q_cc.to('kW'):~P.3f}",
    f"SHR cooling coil @ design = "
    f"{airco_system.SHR_cc.to('pct'):~P.1f}",
    f"mixed air @ design = "
    f"{airco_system.mixed_air.Tdb.to('degC'):~P.1f} DB, "
    f"{airco_system.mixed_air.W.to('g / kg'):~P.1f}",
    f"cooled air @ design = "
    f"{airco_system.cooled_air.Tdb.to('degC'):~P.1f} DB, "
    f"{airco_system.cooled_air.W.to('g / kg'):~P.1f}",
    f"mass flow rate supply air @ design = "
    f"{airco_system.m_supply.to('kg / hr'):~P.1f}",
    f"volume flow rate supply air @ NTP = "
    f"{airco_system.V_supply_ntp.to('m ** 3 / hr'):~P.1f}",
    sep='\n'
)
