"""Air Cooler Selection

In this example, an air conditioning system is designed for the cooling of a
single zone using performance data of an air cooler that was selected from a
catalog.

The aim of the exercise is to determine:
(1) the required supply air mass and volume flow rate to the zone, and
(2) the required face area of the air cooler with known heat exchanger core
geometry, such that the air face velocity remains equal to the air face velocity
specified in the datasheet.
Also, the same conditions on the refrigerant side of the air cooler are
maintained. In this way, it can be assumed that the heat and mass transfer
effectiveness of the air cooler remain unchanged. Then, it is possible to
determine the state of air leaving the air cooler for any state of mixed air
entering the air cooler.

Knowing the state of air leaving the air cooler and being supplied to the zone,
the required mass flow rate of supply air can be determined from a sensible
heat-load balance of the zone, in order to maintain the desired zone air
temperature. The resulting zone air humidity will, however, depend on the latent
load of the zone.

Once the volume flow rate of supply air through the air cooler has been
determined, the required face area of the air cooler can be determined with
the known value of the air face velocity.

As the required supply air mass flow rate is unknown at the beginning, the
mass flow rate of recirculation air will also be unknown. Therefore, to be able
to determine the state of mixed air at the entry of the air cooler, the
ventilation air mass flow rate needs to be specified as a fraction of the still
unknown supply air mass flow rate. At the end of the calculation, the resulting
ventilation flow rate can be determined, which can then be checked against the
ventilation requirements of the zone. If the ventilation requirement is not
fulfilled, we can try with another fraction.
"""
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning.single_zone.air_coil_sizing import AircoSystem

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


# AIRCO-SYSTEM DESIGN WITH KNOWN AIR COOLER SPECS FROM DATASHEET
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
