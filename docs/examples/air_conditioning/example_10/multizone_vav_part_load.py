"""
ANALYSIS OF A MULTI-ZONE VAV SYSTEM UNDER SUMMER PART-LOAD CONDITIONS
-------------------------------------------------------------------

To analyze the summer part-load operation of a multi-zone VAV system, we
use the module `part_load` in subpackage `hvac.air_conditioning.multi_zone`.

In this script, a part-load analysis is demonstrated of system having two zones,
zone A and zone B.

Inside the `__init__()`-method of class `PeakSummerDesign`, we first run the
design calculations under peak summer design conditions. See the line:
`self.results_summer = self.vav_system.design_summer(...)`

This `PeakSummerDesign` object is then passed to the `__init__()`-method of
class `PartLoadSummerAnalysis`, together with a specification about the current
part-load conditions: the state of outdoor air and the sensible and latent zone
loads, which are expressed as fractions of the peak summer design loads. On
instantiation of a `PartLoadSummerAnalysis` object, the part-load analysis
routine is called:
`self.results_summer = self.vav_system.part_load_summer(...)`

`self.results_summer` is a dictionary that contains only the main system
quantities. Inside the `main()` function we access the results through
object attributes.

By encapsulating the part-load analysis inside the `__init__()`-method of a
class, we can easily create a series of "part-load objects" for a range of
different part-load conditions, e.g., determined from TMY data and a thermal
model of the building to calculate the sensible and latent zone loads.
"""
import warnings

from hvac.fluids import CoolPropWarning
warnings.filterwarnings('ignore', category=CoolPropWarning)

from hvac import Quantity
from hvac.fluids import HumidAir, Fluid
from hvac.air_conditioning.multi_zone import Zone, Season
import hvac.air_conditioning.multi_zone.design as design
import hvac.air_conditioning.multi_zone.part_load as analysis


Q_ = Quantity
Water = Fluid('Water')


class PeakSummerDesign:
    """Peak-summer design calculations."""

    def __init__(self):
        # Create zone A:
        # From the design cooling and heating load calculations of the zone, we
        # have the cooling / heating load of the zone under design conditions
        # at the desired state of the zone air.
        self.zone_A = Zone(
            name='zone A',
            summer=Season(
                Q_sen=Q_(24, 'kW'),
                Q_lat=Q_(6, 'kW'),
                zone_air=HumidAir(Tdb=Q_(26, 'degC'), RH=Q_(50, 'pct'))
            ),
            winter=Season(
                Q_sen=Q_(-11.25, 'kW'),
                Q_lat=Q_(0, 'kW'),
                zone_air=HumidAir(Tdb=Q_(22, 'degC'), RH=Q_(50, 'pct'))
            )
        )
        # Create zone B:
        self.zone_B = Zone(
            name='zone B',
            summer=Season(
                Q_sen=Q_(15, 'kW'),
                Q_lat=Q_(5, 'kW'),
                zone_air=HumidAir(Tdb=Q_(26, 'degC'), RH=Q_(50, 'pct'))
            ),
            winter=Season(
                Q_sen=Q_(-5, 'kW'),
                Q_lat=Q_(0, 'kW'),
                zone_air=HumidAir(Tdb=Q_(22, 'degC'), RH=Q_(50, 'pct'))
            )
        )
        # Create the VAV system:
        self.vav_system = design.VAVSystem(
            zones=[self.zone_A, self.zone_B],
            outdoor_air_summer=HumidAir(Tdb=Q_(28.9, 'degC'), Twb=Q_(20.1, 'degC')),
            outdoor_air_winter=HumidAir(Tdb=Q_(-5.7, 'degC'), Twb=Q_(-6.7, 'degC')),
            V_vent=Q_(1225, 'm**3 / hr'),
        )
        # Design the VAV system for peak-summer operation:
        self.results_summer = self.vav_system.design_summer(
            eta_fan_sup=Q_(60, 'pct'),
            dP_fan_sup=Q_(746, 'Pa')
        )


class PartLoadSummerAnalysis:

    def __init__(
        self,
        peak_summer_design: PeakSummerDesign,
        outdoor_air: HumidAir,
        zone_A_frac_Q_sen: float,
        zone_A_frac_Q_lat: float,
        zone_B_frac_Q_sen: float,
        zone_B_frac_Q_lat: float
    ) -> None:
        self.peak_summer_design = peak_summer_design
        self.zone_A = Zone(
            name='zone A',
            summer=Season(
                Q_sen=zone_A_frac_Q_sen * self.peak_summer_design.zone_A.summer.Q_sen,
                Q_lat=zone_A_frac_Q_lat * self.peak_summer_design.zone_A.summer.Q_lat,
                zone_air=self.peak_summer_design.zone_A.summer.zone_air,
                m_supply_des=self.peak_summer_design.zone_A.summer.m_supply
            )
        )
        self.zone_B = Zone(
            name='zone B',
            summer=Season(
                Q_sen=zone_B_frac_Q_sen * self.peak_summer_design.zone_B.summer.Q_sen,
                Q_lat=zone_B_frac_Q_lat * self.peak_summer_design.zone_B.summer.Q_lat,
                zone_air=self.peak_summer_design.zone_B.summer.zone_air,
                m_supply_des=self.peak_summer_design.zone_B.summer.m_supply
            )
        )
        self.vav_system = analysis.VAVSystem(
            zones=[self.zone_A, self.zone_B],
            T_supply_des=self.peak_summer_design.vav_system.summer.T_supply,
            outdoor_air_summer=outdoor_air,
            V_vent=self.peak_summer_design.vav_system.summer.V_vent
        )
        self.results_summer = self.vav_system.part_load_summer(
            eta_fan_sup=Q_(60, 'pct'),
            dP_fan_sup=Q_(746, 'Pa')
        )


def main():
    peak_summer_design = PeakSummerDesign()
    ptl_obj = PartLoadSummerAnalysis(
        peak_summer_design=peak_summer_design,
        outdoor_air=HumidAir(Tdb=Q_(26, 'degC'), RH=Q_(50, 'pct')),
        zone_A_frac_Q_sen=0.75,
        zone_A_frac_Q_lat=1.00,
        zone_B_frac_Q_sen=0.50,
        zone_B_frac_Q_lat=1.00
    )
    Q_cc_tot_des = ptl_obj.peak_summer_design.vav_system.summer.cooling_coil.Q
    V_sup_des = ptl_obj.peak_summer_design.vav_system.summer.V_supply
    Q_cc_tot_ptl = ptl_obj.vav_system.summer.cooling_coil.Q
    V_sup_ptl = ptl_obj.vav_system.summer.V_supply
    print(
        "cooling coil load under peak summer design conditions = "
        f"{Q_cc_tot_des.to('kW'):~P.3f}"
    )
    print(
        "cooling coil load under current part-load conditions = "
        f"{Q_cc_tot_ptl.to('kW'):~P.3f}"
    )
    print(
        f"supply air volume flow rate under peak summer design conditions = "
        f"{V_sup_des.to('m**3/hr'):~P.3f}"
    )
    print(
        f"supply air volume flow rate under current part-load conditions = "
        f"{V_sup_ptl.to('m**3/hr'):~P.3f}"
    )
    for zone in ptl_obj.vav_system.zones:
        print(
            f"{zone.name}: ",
            f"supply air volume flow rate = {zone.summer.V_supply.to('m**3/hr'):~P.3f}",
            f"reheat-coil load = {zone.reheat_coil.Q_sen.to('W'):~P.3f}",
            sep='\n'
        )


if __name__ == '__main__':
    main()
