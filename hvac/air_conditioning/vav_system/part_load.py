from typing import Optional, List, Dict, Tuple
from hvac import Quantity
from ..core import AirConditioningProcess, Fan, AirStream, AdiabaticMixing
from hvac.fluids import HumidAir
from .design import Zone


Q_ = Quantity


class VAVSystem:

    class Summer:

        def __init__(self, system: 'VAVSystem', T_supply: Quantity, outdoor_air: HumidAir, V_vent: Quantity):
            self.system = system
            self.T_supply = T_supply
            self.outdoor_air = outdoor_air
            self.m_vent = V_vent * outdoor_air.rho
            self.m_supply: Optional[Quantity] = None
            self.V_supply: Optional[Quantity] = None
            self.supply_air: Optional[HumidAir] = None
            self.m_return: Optional[Quantity] = None
            self.V_return: Optional[Quantity] = None
            self.return_air: Optional[HumidAir] = None
            self.recirculated_air: Optional[HumidAir] = None
            self.mixed_air: Optional[HumidAir] = None
            self.cooled_air: Optional[HumidAir] = None
            self.cooling_coil: Optional[AirConditioningProcess] = None
            self.Q_reheat: Optional[Quantity] = None

        def determine_m_supply(self):
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    T_ai=self.T_supply,
                    T_ao=zone.summer.zone_air.Tdb,
                    Q_sen=zone.summer.Q_sen
                )
                zone.summer.m_supply = max(0.6 * zone.summer.m_supply, p.m_da + zone.summer.m_exhaust)
            self.m_supply = sum(zone.summer.m_supply for zone in self.system.zones)

        def determine_supply_air(self):
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    air_out=zone.summer.zone_air,
                    Q_sen=zone.summer.Q_sen,
                    Q_lat=zone.summer.Q_lat,
                    m_da=zone.summer.m_supply
                )
                zone.summer.supply_air = p.air_in
            W_avg = sum(z.summer.supply_air.W for z in self.system.zones) / len(self.system.zones)
            self.supply_air = HumidAir(Tdb=self.T_supply, W=W_avg)
            self.V_supply = self.m_supply * self.supply_air.v

        def determine_cooling_coil_exit_air(
                self,
                supply_fan_pressure: Optional[Quantity] = None,
                supply_fan_efficiency: Optional[Quantity] = None,
                supply_duct_heat_gain: Optional[Quantity] = None
        ):
            if supply_fan_pressure is not None:
                fan = Fan(
                    air_out=self.supply_air,
                    fan_efficiency=supply_fan_efficiency if supply_fan_efficiency is not None else Q_(100, 'pct'),
                    fan_pressure=supply_fan_pressure
                )
                dT_supply_fan = fan.air_out.Tdb - fan.air_in.Tdb
            else:
                dT_supply_fan = Q_(0.0, 'K')

            if supply_duct_heat_gain is not None:
                supply_duct = AirConditioningProcess(
                    T_ao=self.supply_air.Tdb,
                    m_da=self.m_supply,
                    Q_sen=supply_duct_heat_gain
                )
                dT_supply_duct = supply_duct.T_ao - supply_duct.T_ai
            else:
                dT_supply_duct = Q_(0.0, 'K')

            T_cold = self.T_supply - dT_supply_fan - dT_supply_duct
            self.cooled_air = HumidAir(Tdb=T_cold, RH=self.supply_air.RH)

        def determine_return_air(
                self,
                return_fan_pressure: Optional[Quantity] = None,
                return_fan_efficiency: Optional[Quantity] = None,
                return_duct_heat_gain: Optional[Quantity] = None
        ):
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    air_in=self.supply_air,
                    T_ao=zone.summer.zone_air.Tdb,
                    SHR=zone.summer.Q_sen / (zone.summer.Q_sen + zone.summer.Q_lat)
                )
                zone.summer.return_air = p.air_out

            self.m_return = sum(zone.summer.m_return for zone in self.system.zones)
            h_return = sum(z.summer.m_return * z.summer.return_air.h for z in self.system.zones) / self.m_return
            W_return = sum(z.summer.m_return * z.summer.return_air.W for z in self.system.zones) / self.m_return
            self.return_air = HumidAir(h=h_return, W=W_return)
            self.V_return = self.m_return * self.return_air.v

            if return_fan_pressure is not None:
                fan = Fan(
                    air_out=self.return_air,
                    fan_efficiency=return_fan_efficiency if return_fan_efficiency is not None else Q_(100, 'pct'),
                    fan_pressure=return_fan_pressure
                )
                dT_return_fan = fan.air_out.Tdb - fan.air_in.Tdb
            else:
                dT_return_fan = Q_(0.0, 'K')

            if return_duct_heat_gain is not None:
                return_duct = AirConditioningProcess(
                    T_ai=self.return_air.Tdb,
                    m_da=self.m_return,
                    Q_sen=return_duct_heat_gain
                )
                dT_return_duct = return_duct.T_ao - return_duct.T_ai
            else:
                dT_return_duct = Q_(0.0, 'K')

            T_recirculated = self.return_air.Tdb + dT_return_fan + dT_return_duct
            self.recirculated_air = HumidAir(Tdb=T_recirculated, RH=self.return_air.RH)

        def determine_mixed_air(self):
            m_recirculated = self.m_return - self.m_vent
            m_exhaust = sum(z.summer.m_exhaust for z in self.system.zones)
            mixing_chamber = AdiabaticMixing(
                in1=AirStream(state=self.recirculated_air, m_da=m_recirculated),
                in2=AirStream(state=self.outdoor_air, m_da=self.m_vent + m_exhaust),
                out=AirStream(m_da=self.m_supply)
            )
            self.mixed_air = mixing_chamber.stream_out.state

        def determine_cooling_coil(self):
            self.cooling_coil = AirConditioningProcess(
                air_in=self.mixed_air,
                air_out=self.cooled_air,
                m_da=self.m_supply,
                h_w=Q_(0.0, 'J / kg')
            )

        def determine_reheat_coils(self):
            for zone in self.system.zones:
                zone.reheat_coil = AirConditioningProcess(
                    T_ai=self.T_supply,
                    T_ao=zone.summer.supply_air.Tdb,
                    m_da=zone.summer.m_supply
                )
            self.Q_reheat = sum(z.reheat_coil.Q_sen for z in self.system.zones)

    def __init__(
            self,
            zones: List[Zone],
            T_supply: Quantity,
            outdoor_air: HumidAir,
            V_vent: Quantity
    ):
        self.zones = zones
        self.summer = VAVSystem.Summer(self, T_supply, outdoor_air, V_vent)

    def part_load_summer(self, **kwargs) -> Dict[str, Quantity]:
        eta_supply_fan = kwargs.get('supply_fan_efficiency')
        dP_supply_fan = kwargs.get('supply_fan_pressure')
        Q_supply_duct = kwargs.get('supply_duct_heat_gain')
        eta_return_fan = kwargs.get('return_fan_efficiency')
        dP_return_fan = kwargs.get('return_fan_pressure')
        Q_return_duct = kwargs.get('return_duct_heat_gain')

        self.summer.determine_m_supply()
        self.summer.determine_supply_air()
        self.summer.determine_cooling_coil_exit_air(dP_supply_fan, eta_supply_fan, Q_supply_duct)
        self.summer.determine_return_air(dP_return_fan, eta_return_fan, Q_return_duct)
        self.summer.determine_mixed_air()
        self.summer.determine_cooling_coil()
        self.summer.determine_reheat_coils()

        results = {
            'cooling coil load': self.summer.cooling_coil.Q,
            'sensible cooling coil load': self.summer.cooling_coil.Q_sen,
            'latent cooling coil load': self.summer.cooling_coil.Q_lat,
            'total reheat coil load': self.summer.Q_reheat,
            'supply air volume flow rate': self.summer.V_supply,
            'return air volume flow rate': self.summer.V_return,
            'system supply air temperature': self.summer.supply_air.Tdb,
            'system return air temperature': self.summer.return_air.Tdb
        }
        return results

    @staticmethod
    def show_results_markdown(results: Dict[str, Quantity], units: Dict[str, Tuple[str, int]]) -> List[str]:
        results_as_string = []
        for k, v in results.items():
            if v.dimensionality == Q_(1, 'W').dimensionality:
                result_str = f"{k}: <b>{v.to(units['Q'][0]):~P.{units['Q'][1]}f}</b>"
            elif v.dimensionality == Q_(1, 'm ** 3 / s').dimensionality:
                result_str = f"{k}: <b>{v.to(units['V'][0]):~P.{units['V'][1]}f}</b>"
            elif v.dimensionality == Q_(1, 'K').dimensionality:
                result_str = f"{k}: <b>{v.to(units['K'][0]):~P.{units['K'][1]}f}</b>"
            else:
                result_str = f"{k}: <b>{v}</b>"
            results_as_string.append(result_str)
        return results_as_string
