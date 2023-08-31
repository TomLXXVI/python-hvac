"""Design calculations for a single-zone air system on the winter peak design
day (heating and humidifying the supply air to the zone).

This script follows the example taken from the book "Heating and Cooling of
Buildings" by T. Agami Reddy, Jan F. Kreider, Peter S. Curtiss and Ari
Rabl (3rd Edition, 2017), Chapter 19: Example 19.4.
"""
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning.single_zone.heating_design import AircoSystem
from hvac.charts import PsychrometricChart, StatePoint


Q_ = Quantity


airco_system = AircoSystem(
    zone_air=HumidAir(Tdb=Q_(72, 'degF'), RH=Q_(50, 'pct')),
    outdoor_air=HumidAir(Tdb=Q_(40, 'degF'), RH=Q_(40, 'pct')),
    Q_dot_zone=Q_(-150_000, 'Btu / hr'),  # don't forget the minus sign for heat losses!
    SHR_zone=Q_(0.8, 'frac'),
    V_dot_vent_ntp=Q_(1000, 'ft ** 3 / min'),
    m_dot_supply_cooling=Q_(17_140, 'lb / hr'),  # design value for cooling operation
    T_supply_max=Q_(105, 'degF'),  # to avoid stratification
    T_steam=Q_(200, 'degF')  # temperature of injected steam
)

print(
    "required steam mass flow rate for humidification = "
    f"{airco_system.m_dot_steam.to('lb / hr'):~P.2f}",
    "required heating capacity of preheat-coil = "
    f"{airco_system.Q_dot_preheat.to('Btu / hr')}",
    "required NTP volume flow rate of supply air = "
    f"{airco_system.V_dot_supply_ntp.to('ft ** 3 / min'):~P.2f}",
    "mixed air = "
    f"{airco_system.mixed_air.Tdb.to('degF'):~P.2f} DB, "
    f"{airco_system.mixed_air.W.to('lb / lb'):~P.5f} AH",
    "supply air = "
    f"{airco_system.supply_air.Tdb.to('degF'):~P.2f} DB, "
    f"{airco_system.supply_air.W.to('lb / lb'):~P.5f} AH",
    "preheated air = "
    f"{airco_system.preheated_air.Tdb.to('degF'):~P.2f} DB, "
    f"{airco_system.preheated_air.W.to('lb / lb'):~P.5f} AH",
    sep='\n'
)

psy_chart = PsychrometricChart()
psy_chart.plot_process(
    name='mixing',
    start_point=StatePoint(
        airco_system.outdoor_air.Tdb,
        airco_system.outdoor_air.W
    ),
    end_point=StatePoint(
        airco_system.return_air.Tdb,
        airco_system.return_air.W
    ),
    mix_point=StatePoint(
        airco_system.mixed_air.Tdb,
        airco_system.mixed_air.W
    )
)
psy_chart.plot_process(
    name='heating',
    start_point=StatePoint(
        airco_system.mixed_air.Tdb,
        airco_system.mixed_air.W
    ),
    end_point=StatePoint(
        airco_system.preheated_air.Tdb,
        airco_system.preheated_air.W
    )
)
psy_chart.plot_process(
    name='humidification',
    start_point=StatePoint(
        airco_system.preheated_air.Tdb,
        airco_system.preheated_air.W
    ),
    end_point=StatePoint(
        airco_system.supply_air.Tdb,
        airco_system.supply_air.W
    )
)
psy_chart.plot_process(
    name='zone',
    start_point=StatePoint(
        airco_system.supply_air.Tdb,
        airco_system.supply_air.W
    ),
    end_point=StatePoint(
        airco_system.zone_air.Tdb,
        airco_system.zone_air.W
    )
)
psy_chart.show()
