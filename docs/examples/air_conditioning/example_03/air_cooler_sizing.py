from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning.single_zone.air_cooling import AircoSystem

Q_ = Quantity


# AIR COOLER SPECIFICATIONS FROM DATASHEET
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


# AIRCO-SYSTEM DESIGN WITH KNOWN AIR COOLER
Q_sen_zone = Q_(30.514, 'kW')
Q_lat_zone = Q_(12.967, 'kW')
Q_zone = Q_sen_zone + Q_lat_zone
SHR_zone = Q_sen_zone / Q_zone

zone_air = HumidAir(Tdb=Q_(26, 'degC'), W=Q_(8.0, 'g / kg'))
outdoor_air = HumidAir(Tdb=Q_(31, 'degC'), W=Q_(10.0, 'g / kg'))


airco_system = AircoSystem(
    zone_air=zone_air,
    outdoor_air=outdoor_air,
    Q_zone=Q_zone,
    SHR_zone=SHR_zone,
    f_vent=Q_(36, 'pct'),
    air_cooler=AircoSystem.WetDXAirCooler(
        air_in_ref=Specs.air_in,
        air_out_ref=Specs.air_out,
        T_rfg=Specs.T_rfg,
        v_fa_ntp=Specs.v_fa_ntp,
        dP_ntp=Specs.dP_ntp
    )
)
print(airco_system.info())
