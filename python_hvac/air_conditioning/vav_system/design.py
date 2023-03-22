from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass, field
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning import (
    AirConditioningProcess,
    Fan,
    AirStream,
    AdiabaticMixing,
    SpaceConditionLine
)

Q_ = Quantity


@dataclass
class Season:
    Q_sen: Quantity
    Q_lat: Quantity
    zone_air: HumidAir
    m_exhaust: Quantity = Q_(0.0, 'kg / s')
    m_supply: Quantity = Q_(float('nan'), 'kg / s')
    supply_air: Optional[HumidAir] = field(init=False, default=None)
    return_air: Optional[HumidAir] = field(init=False, default=None)

    @property
    def m_return(self) -> Quantity:
        return self.m_supply - self.m_exhaust

    @property
    def V_supply(self) -> Quantity:
        return self.m_supply * self.supply_air.v


@dataclass
class Zone:
    name: str
    summer: Optional[Season] = None
    winter: Optional[Season] = None
    reheat_coil: Optional[AirConditioningProcess] = field(init=False, default=None)


class VAVSystem:

    class Summer:
        dT_supply = Q_(12.0, 'K')

        def __init__(self, outdoor_air: HumidAir, V_vent: Quantity, system: 'VAVSystem'):
            self.outdoor_air = outdoor_air
            self.V_vent = V_vent
            self.m_vent = V_vent * outdoor_air.rho
            self.system = system
            self.T_supply: Quantity = Q_(float('nan'), 'degC')
            self.supply_air: Optional[HumidAir] = None
            self.m_supply: Quantity = Q_(float('nan'), 'kg /s')
            self.V_supply: Quantity = Q_(float('nan'), 'kg /s')
            self.T_cold: Quantity = Q_(float('nan'), 'degC')
            self.cooled_air: Optional[HumidAir] = None
            self.m_return: Quantity = Q_(float('nan'), 'kg /s')
            self.V_return: Quantity = Q_(float('nan'), 'kg /s')
            self.return_air: Optional[HumidAir] = None
            self.recirculated_air: Optional[HumidAir] = None
            self.mixed_air: Optional[HumidAir] = None
            self.cooling_coil: Optional[AirConditioningProcess] = None
            self.m_supply_part_load: Quantity = Q_(float('nan'), 'kg /s')
            self.V_supply_part_load: Quantity = Q_(float('nan'), 'kg /s')

        def determine_supply_air(self):
            T_zone_avg = sum(zone.summer.zone_air.Tdb for zone in self.system.zones) / len(self.system.zones)
            W_zone_avg = sum(zone.summer.zone_air.W for zone in self.system.zones) / len(self.system.zones)
            zone_air_avg = HumidAir(Tdb=T_zone_avg, W=W_zone_avg)
            Q_sen_tot = sum(zone.summer.Q_sen for zone in self.system.zones)
            Q_lat_tot = sum(zone.summer.Q_lat for zone in self.system.zones)
            SCL = SpaceConditionLine(zone_air_avg, Q_sen_tot, Q_lat_tot)
            self.T_supply = T_zone_avg - self.dT_supply
            W_supply = SCL.W_ai(self.T_supply)
            self.supply_air = HumidAir(Tdb=self.T_supply, W=W_supply)
            for zone in self.system.zones:
                zone.summer.supply_air = self.supply_air

        def determine_m_supply(self):
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    T_ai=self.T_supply,
                    T_ao=zone.summer.zone_air.Tdb,
                    Q_sen=zone.summer.Q_sen,
                )
                zone.summer.m_supply = p.m_da
                zone.summer.m_supply += zone.summer.m_exhaust
            self.m_supply = sum(zone.summer.m_supply for zone in self.system.zones)
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

            self.T_cold = self.T_supply - dT_supply_fan - dT_supply_duct
            self.cooled_air = HumidAir(Tdb=self.T_cold, W=self.supply_air.W)
        
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
                    air_in=self.return_air,
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

    class Winter:
        T_supply_max = Q_(40.0, 'degC')
        RH_supply = Q_(0.0, 'pct')

        def __init__(self, outdoor_air: HumidAir, V_vent: Quantity, system: 'VAVSystem'):
            self.system = system
            self.outdoor_air = outdoor_air
            self.V_vent = V_vent
            self.m_vent = V_vent * self.outdoor_air.rho
            self.preheat_coil: Optional[AirConditioningProcess] = None
            self.Q_preheat_peak_load: Quantity = Q_(float('nan'), 'W')
            self.m_supply: Quantity = Q_(float('nan'), 'kg / s')
            self.V_supply: Quantity = Q_(float('nan'), 'kg / s')
            self.T_supply: Quantity = Q_(float('nan'), 'degC')
            self.supply_air: Optional[HumidAir] = None
            self.m_return: Quantity = Q_(float('nan'), 'kg /s')
            self.V_return: Quantity = Q_(float('nan'), 'kg / s')
            self.return_air: Optional[HumidAir] = None
            self.recirculated_air: Optional[HumidAir] = None
            self.m_recirculated: Quantity = Q_(float('nan'), 'kg /s')
            self.mixed_air: Optional[HumidAir] = None
            self.preheated_air: Optional[HumidAir] = None
            self.T_cold: Quantity = Q_(float('nan'), 'degC')
            self.cooled_air: Optional[HumidAir] = None
            self.cooling_coil: Optional[AirConditioningProcess] = None
            self.Q_reheat: Quantity = Q_(float('nan'), 'W')

        def determine_preheat_peak_load(self):
            preheat_coil = AirConditioningProcess(
                T_ai=self.outdoor_air.Tdb,
                T_ao=self.system.summer.T_cold,
                m_da=self.m_vent
            )
            self.Q_preheat_peak_load = preheat_coil.Q_sen

        def determine_m_supply(self):
            for zone in self.system.zones:
                if zone.winter.Q_sen < 0.0:  # zone requires heating
                    p = AirConditioningProcess(
                        T_ai=self.T_supply_max,
                        T_ao=zone.winter.zone_air.Tdb,
                        Q_sen=zone.winter.Q_sen
                    )
                else:  # zone requires cooling
                    T_supply = zone.winter.zone_air.Tdb - self.system.summer.dT_supply
                    p = AirConditioningProcess(
                        T_ai=T_supply,
                        T_ao=zone.winter.zone_air.Tdb,
                        Q_sen=zone.winter.Q_sen
                    )
                zone.winter.m_supply = p.m_da + zone.winter.m_exhaust
                zone.winter.m_supply = max(zone.winter.m_supply, 0.6 * zone.summer.m_supply)
            self.m_supply = sum(z.winter.m_supply for z in self.system.zones)

        def determine_T_supply(self):
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    T_ao=zone.winter.zone_air.Tdb,
                    m_da=zone.winter.m_supply,
                    Q_sen=zone.winter.Q_sen
                )
                zone.winter.supply_air = HumidAir(Tdb=p.T_ai, RH=self.RH_supply)
            self.T_supply = min(z.winter.supply_air.Tdb for z in self.system.zones)
            self.supply_air = HumidAir(Tdb=self.T_supply, RH=self.RH_supply)
            self.V_supply = self.m_supply * self.supply_air.v

        def determine_return_air(
                self,
                return_fan_pressure: Optional[Quantity] = None,
                return_fan_efficiency: Optional[Quantity] = None,
                return_duct_heat_gain: Optional[Quantity] = None
        ):
            self.m_return = sum(z.winter.m_return for z in self.system.zones)
            h_return = sum(z.winter.zone_air.h * z.winter.m_return for z in self.system.zones) / self.m_return
            W_return = sum(z.winter.zone_air.W * z.winter.m_return for z in self.system.zones) / self.m_return
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
                    T_ao=self.return_air.Tdb,
                    m_da=self.m_return,
                    Q_sen=return_duct_heat_gain
                )
                dT_return_duct = return_duct.T_ao - return_duct.T_ai
            else:
                dT_return_duct = Q_(0.0, 'K')

            T_recirculated = self.return_air.Tdb + dT_return_fan + dT_return_duct
            self.recirculated_air = HumidAir(Tdb=T_recirculated, RH=self.return_air.RH)

        def determine_mixed_air(self):
            self.m_recirculated = self.m_return - self.m_vent
            m_exhaust = sum(z.winter.m_exhaust for z in self.system.zones)
            mixing_chamber = AdiabaticMixing(
                in1=AirStream(state=self.recirculated_air, m_da=self.m_recirculated),
                in2=AirStream(state=self.outdoor_air, m_da=self.m_vent + m_exhaust),
                out=AirStream(m_da=self.m_supply)
            )
            self.mixed_air = mixing_chamber.stream_out.state

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

            self.T_cold = self.T_supply - dT_supply_fan - dT_supply_duct
            self.cooled_air = HumidAir(Tdb=self.T_cold, RH=self.RH_supply)
            if self.cooled_air.Tdb > self.mixed_air.Tdb:
                # the required cooling coil leaving temperature is greater than the mixed air temperature at the cooling
                # coil entrance ==> preheat the mixed air to be equal to the set cooled air temperature.
                self.preheat_coil = AirConditioningProcess(
                    air_in=self.mixed_air,
                    air_out=self.cooled_air,
                    m_da=self.m_supply
                )
                self.preheated_air = self.cooled_air
            else:
                self.preheated_air = self.mixed_air

        def determine_cooling_coil(self):
            self.cooling_coil = AirConditioningProcess(
                air_in=self.preheated_air,
                air_out=self.cooled_air,
                m_da=self.m_supply,
                h_w=Q_(0.0, 'J / kg')
            )

        def determine_reheat_coils(self):
            for zone in self.system.zones:
                zone.reheat_coil = AirConditioningProcess(
                    T_ai=self.T_supply,
                    T_ao=zone.winter.supply_air.Tdb,
                    m_da=zone.winter.m_supply
                )
            self.Q_reheat = sum(z.reheat_coil.Q_sen for z in self.system.zones)

    def __init__(
            self,
            zones: List[Zone],
            outdoor_air_summer: HumidAir,
            outdoor_air_winter: HumidAir,
            V_vent: Quantity
    ):
        self.zones = zones
        self.summer = VAVSystem.Summer(outdoor_air_summer, V_vent, self)
        self.winter = VAVSystem.Winter(outdoor_air_winter, V_vent, self)

    def design_summer(self, **kwargs) -> Dict[str, Quantity]:
        dT_supply = kwargs.get('dT_supply')

        eta_supply_fan = kwargs.get('supply_fan_efficiency')
        dP_supply_fan = kwargs.get('supply_fan_pressure')
        Q_supply_duct = kwargs.get('supply_duct_heat_gain')

        eta_return_fan = kwargs.get('return_fan_efficiency')
        dP_return_fan = kwargs.get('return_fan_pressure')
        Q_return_duct = kwargs.get('return_duct_heat_gain')

        if dT_supply is not None: self.summer.dT_supply = dT_supply

        self.summer.determine_supply_air()
        self.summer.determine_m_supply()
        self.summer.determine_cooling_coil_exit_air(dP_supply_fan, eta_supply_fan, Q_supply_duct)
        self.summer.determine_return_air(dP_return_fan, eta_return_fan, Q_return_duct)
        self.summer.determine_mixed_air()
        self.summer.determine_cooling_coil()

        results = {
            'cooling coil load': self.summer.cooling_coil.Q,
            'sensible cooling coil load': self.summer.cooling_coil.Q_sen,
            'latent cooling coil load': self.summer.cooling_coil.Q_lat,
            'supply air volume flow rate': self.summer.V_supply,
            'return air volume flow rate': self.summer.V_return,
            'system supply air temperature': self.summer.supply_air.Tdb,
            'system return air temperature': self.summer.return_air.Tdb
        }
        return results

    def design_winter(self, **kwargs) -> Dict[str, Quantity]:
        T_supply_max = kwargs.get('T_supply_max')
        RH_supply = kwargs.get('RH_supply')

        eta_supply_fan = kwargs.get('supply_fan_efficiency')
        dP_supply_fan = kwargs.get('supply_fan_pressure')
        Q_supply_duct = kwargs.get('supply_duct_heat_gain')

        eta_return_fan = kwargs.get('return_fan_efficiency')
        dP_return_fan = kwargs.get('return_fan_pressure')
        Q_return_duct = kwargs.get('return_duct_heat_gain')

        if T_supply_max is not None: self.winter.T_supply_max = T_supply_max
        if RH_supply is not None: self.winter.RH_supply = RH_supply

        self.winter.determine_preheat_peak_load()
        self.winter.determine_m_supply()
        self.winter.determine_T_supply()
        self.winter.determine_return_air(dP_return_fan, eta_return_fan, Q_return_duct)
        self.winter.determine_mixed_air()
        self.winter.determine_cooling_coil_exit_air(dP_supply_fan, eta_supply_fan, Q_supply_duct)
        self.winter.determine_cooling_coil()
        self.winter.determine_reheat_coils()

        results = {
            'preheat coil peak load': self.winter.Q_preheat_peak_load,
            'preheat coil load': self.winter.preheat_coil.Q_sen,
            'total reheat coil load': self.winter.Q_reheat,
            'cooling coil load': self.winter.cooling_coil.Q,
            'sensible cooling coil load': self.winter.cooling_coil.Q_sen,
            'latent cooling coil load': self.winter.cooling_coil.Q_lat,
            'supply air volume flow rate': self.winter.V_supply,
            'return air volume flow rate': self.winter.V_return,
            'system supply air temperature': self.winter.supply_air.Tdb,
            'system return air temperature': self.winter.return_air.Tdb
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
