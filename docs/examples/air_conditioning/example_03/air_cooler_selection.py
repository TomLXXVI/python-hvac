"""AIR COOLER SELECTION AND SIZING

In this example, an air conditioning system is designed for cooling a single
zone using the performance data of an air cooler selected from a catalog.

The aim is to determine:
(1) the required supply airflow to the zone, and
(2) the required face area of the air cooler.

The air face velocity remains equal to the air face velocity specified in the
datasheet of the selected air cooler. Also, the same conditions on the
refrigerant side of the air cooler need to be maintained.
In this way, it can be assumed that the heat and mass transfer effectiveness of
the air cooler remain unchanged.
Then, it is possible to determine the state of air that leaves the air cooler
for any state of mixed air that will enter the air cooler.

Knowing the state of air leaving the air cooler and being supplied to the zone,
the required mass (and volume) flow rate of this supply air can then be determined
from a sensible heat-load balance of the zone, in order to maintain the desired
zone air temperature. The resulting zone air humidity will, however, depend on
the latent cooling load of the zone.

Once the volume flow rate of supply air through the air cooler has been
determined, the required face area of the air cooler can be determined with
the known value of the air face velocity.

As the required supply air mass flow rate is still unknown at the beginning, the
mass flow rate of recirculation air will also be unknown. Therefore, to be able
to determine the state of mixed air at the entry of the air cooler, the
ventilation air mass flow rate needs to be specified as a fraction of the still
unknown supply air mass flow rate. At the end of the calculations, the resulting
ventilation flow rate can be determined, which can then be checked against the
ventilation requirements of the zone. If the ventilation requirement is not
fulfilled, we can try with a higher fraction.
"""
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning.single_zone.cooling_coil_sizing import DxAirCoolingCoilSizer


Q_ = Quantity


# AIR COOLER SPECIFICATIONS FROM DATASHEET
# ========================================
# The air cooler selected from the catalog has a cooling capacity of 56 kW under
# the following conditions:
# - state of entering air: Tdb 27.3 °C, RH 80 %
# - state of leaving air: Tdb 18.5 °C, RH 90 %
# - refrigerant evaporation temperature: 5 °C
# - volume flow rate of air at NTP (dry air at 20 °C and 101.325 kPa): 7500 m³/h
# - air face velocity: 2.37 m/s
# - air-side pressure drop: 50 Pa

class Specs:
    Q_cc = Q_(56, 'kW')
    air_in = HumidAir(
        Tdb=Q_(27.3, 'degC'),
        RH=Q_(80.0, 'pct')
    )
    air_out = HumidAir(
        Tdb=Q_(18.5, 'degC'),
        RH=Q_(98, 'pct')
    )
    T_rfg = Q_(5.0, 'degC')
    V_ntp = Q_(7500, 'm ** 3 / hr')
    v_fa_ntp = Q_(2.37, 'm / s')
    dP_ntp = Q_(50, 'Pa')


air_cooler = DxAirCoolingCoilSizer.WetDxAirCooler(
    air_in_ref=Specs.air_in,
    air_out_ref=Specs.air_out,
    T_rfg=Specs.T_rfg,
    v_fa_ntp=Specs.v_fa_ntp,
    dP_ntp=Specs.dP_ntp
)


# SUMMARY OF THE ZONE'S COOLING-LOAD CALCULATION
# ==============================================
# From the cooling load calculation of the zone, it follows that under summer
# peak design conditions the sensible cooling load of the zone is 31.126 kW and
# the latent cooling load is 17.936 kW.

Q_sen_zone = Q_(31.126, 'kW')
Q_lat_zone = Q_(17.936, 'kW')
Q_zone = Q_sen_zone + Q_lat_zone
SHR_zone = Q_sen_zone / Q_zone

# The cooling load calculation was based on this zone air state:
zone_air = HumidAir(Tdb=Q_(26.0, 'degC'), RH=Q_(50.0, 'pct'))

# At peak summer design conditions, the outdoor air state is:
outdoor_air = HumidAir(Tdb=Q_(26.7, 'degC'), Twb=Q_(19.2, 'degC'))


# CONFIGURE THE AC SYSTEM
# =======================
# To configure the AC system, we need to:
# - specify the state of zone air for which the AC system will be designed
# - specify the state of outdoor air at design conditions
# - specify the total cooling load of the zone
# - specify the sensible heat ratio of the total zone load
# - specify the fraction of outdoor ventilation air in the total supply air flow
#   to the zone
# - configure the air-cooling coil using the known specs from the datasheet

# Assume a fraction of the supply air mass flow rate that is outdoor ventilation
# air:
f_vent = Q_(40, 'pct')


airco_system = DxAirCoolingCoilSizer(
    zone_air=zone_air,
    outdoor_air=outdoor_air,
    Q_zone=Q_zone,
    SHR_zone=SHR_zone,
    f_vent=f_vent,
    air_cooler=air_cooler
)

# RESULTS
# =======
print(airco_system.info())
