from typing import List, Optional, Dict
from dataclasses import dataclass, field
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState
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
    """Dataclass that holds the design day data and calculation results of a
    zone for, either the cooling season (summer), or the heating season
    (winter).

    Attributes
    ----------
    Q_sen:
        Sensible cooling or heating load of the zone.
    Q_lat:
        Latent load of the zone.
    zone_air:
        Desired zone air state for the design of the system.
    m_exhaust: default 0.0 kg/s
        Mass flow rate of air that is locally exhausted from the zone.
    m_supply:
        Required mass flow rate of supply air to the zone (calculation result).
    supply_air:
        Required state of supply air (calculation result).
    return_air:
        Resulting state of return air (calculation result).
    m_return:
        Resulting mass flow rate of return air (calculation result)
    V_supply:
        Required volume flow rate of supply air to the zone (calculation
        result).
    """
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
    """Dataclass that represents a temperature zone of a building.

    Attributes
    ----------
    name:
        Identifier for the zone.
    summer: optional
        `Season` object that holds the summer design data and calculation
        results of the zone.
    winter: optional
        `Season` object that holds the winter design data and calculation
        results of the zone.
    reheat_coil: optional.
        Reference to the sensible heating process in the reheat coil of the
        zone.
    """
    name: str
    summer: Optional[Season] = None
    winter: Optional[Season] = None
    reheat_coil: Optional[AirConditioningProcess] = field(
        init=False,
        default=None
    )


class VAVSystem:
    """Class that incorporates the design routines of a VAV system."""

    class Summer:
        """Class for designing the VAV system for the peak summer day.

        Attributes
        ----------
        dT_sup:
            Allowable temperature difference between the zone air and the supply
            air to ensure proper mixing.
        outdoor_air:
            State of outdoor air on the summer peak design day.
        V_vent:
            Required ventilation air volume flow rate of the building.
        T_supply:
            Temperature of supply air to the zones. It is calculated by taking
            the average zone air temperature minus `dT_sup`, which is by
            default set at 12 K.
        supply_air:
            State of supply air to the zones.
        m_supply:
            Mass flow rate of supply air to the zones.
        V_supply:
            Volume flow rate of supply air to the zones (referred to the
            state of the supply air downstream of the supply fan).
        T_cold:
            Temperature of air leaving the cooling coil.
        cooled_air:
            State of air leaving the cooling coil.
        m_return:
            Mass flow rate of return air from the zones.
        V_return:
            Volume flow rate of return air from the zones (referred to the
            state of the return air).
        return_air:
            State of return air from the zones.
        recirculated_air:
            State of air being recirculated to the mixing chamber of the system.
        mixed_air:
            State of air resulting from adiabatic mixing in the mixing chamber.
        cooling_coil:
            Refers to the cooling and dehumidification process in the cooling
            coil (instance of `AirConditioningProcess`).
        m_supply_part_load:
            Mass flow rate of supply air to the zones at part load.
        V_supply_part_load:
            Volume flow rate of supply air to the zones at part load.
        """
        dT_sup = Q_(12.0, 'K')

        def __init__(
            self,
            outdoor_air: HumidAir,
            V_vent: Quantity,
            system: 'VAVSystem'
        ) -> None:
            """Creates the `Summer` design stage of the VAV system.

            Parameters
            ----------
            outdoor_air:
                State of outdoor air on the summer peak design day.
            V_vent:
                Required ventilation air volume flow rate of the building.
            system:
                Refers to the `VAVSystem` parent class.
            """
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

        def determine_supply_air(self) -> None:
            """Calculates the wanted state of system supply air to all zones.
            The system supply air temperature is determined by taking the
            average of the desired zone air temperatures from which a
            predetermined temperature difference (that still ensures proper
            mixing) is subtracted.
            The average humidity ratio of all zones is calculated. The average
            zone air temperature and humidity ratio determine the average state
            of the zone air in the building.
            The totals of the sensible and the latent zone loads are calculated.
            Applying the space condition line equation to the average zone air,
            delivers the humidity ratio of the system supply air that goes with
            the system supply air temperature.
            """
            # Average zone air temperature:
            T_zone_avg = sum(
                zone.summer.zone_air.Tdb
                for zone in self.system.zones
            ) / len(self.system.zones)
            # Average zone air humidity ratio:
            W_zone_avg = sum(
                zone.summer.zone_air.W
                for zone in self.system.zones
            ) / len(self.system.zones)
            # Average state of zone air:
            zone_air_avg = HumidAir(Tdb=T_zone_avg, W=W_zone_avg)
            # Total sensible cooling load of all zones:
            Q_sen_tot = sum(
                zone.summer.Q_sen
                for zone in self.system.zones
            )
            # Total latent load of all zones:
            Q_lat_tot = sum(
                zone.summer.Q_lat
                for zone in self.system.zones
            )
            # Average space condition line of all zones:
            scl = SpaceConditionLine(zone_air_avg, Q_sen_tot, Q_lat_tot)
            # Selection of supply air temperature to the zones:
            self.T_supply = T_zone_avg - self.dT_sup
            # Resulting humidity ratio determined from the average space
            # condition line of the zones:
            W_supply = scl.W_ai(self.T_supply)
            # State of supply air to the zones:
            self.supply_air = HumidAir(Tdb=self.T_supply, W=W_supply)
            # Set the state of supply air of each zone in the VAV system (reheat
            # coils are off):
            for zone in self.system.zones:
                zone.summer.supply_air = self.supply_air

        def determine_m_supply(self) -> None:
            """Calculates the required mass flow rate of supply air to each zone
            using the sensible heat balance of the zone and the system supply
            air temperature previously determined.
            Any local exhaust air mass flow rate is added to get at the total
            mass flow rate of air that must be supplied to each zone.
            By taking the sum, the total mass and volume flow rate of system
            supply air to the zones is determined.
            """
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    T_ai=self.T_supply,
                    T_ao=zone.summer.zone_air.Tdb,
                    Q_sen=zone.summer.Q_sen,
                )
                zone.summer.m_supply = p.m_da
                # Any air that is exhausted locally in a zone must also be
                # supplied to the zone:
                zone.summer.m_supply += zone.summer.m_exhaust

            self.m_supply = sum(
                zone.summer.m_supply
                for zone in self.system.zones
            )
            self.V_supply = self.m_supply * self.supply_air.v

        def determine_cooled_air(
            self,
            dP_fan_sup: Optional[Quantity] = None,
            eta_fan_sup: Optional[Quantity] = None,
            Q_duct_sup: Optional[Quantity] = None
        ) -> None:
            """Determines the required state of air leaving the cooling coil of
            the VAV system, possibly taking the temperature rise of air through
            the supply fan and any heat gain along the supply duct into account.

            Parameters
            ----------
            dP_fan_sup: optional, default None
                The required pressure gain of the supply fan needed to deliver
                the required mass flow rate of supply air to the zones,
                determined by a pressure loss calculation of the air supply duct
                system.
            eta_fan_sup: optional, default None
                The efficiency of the supply fan at the required operating
                point of the fan.
            Q_duct_sup: optional, default None
                Any heat added from the environment to the supply air along the
                supply duct between the cooling coil and the zones.
            """
            # Temperature rise due to fan heating:
            if dP_fan_sup is not None:
                fan = Fan(
                    air_out=self.supply_air,
                    eta_fan=(
                        eta_fan_sup
                        if eta_fan_sup is not None
                        else Q_(100, 'pct')
                    ),
                    dP_fan=dP_fan_sup
                )
                dT_supply_fan = fan.air_out.Tdb - fan.air_in.Tdb
            else:
                dT_supply_fan = Q_(0.0, 'K')
            # Temperature rise due to duct heat gain:
            if Q_duct_sup is not None:
                supply_duct = AirConditioningProcess(
                    T_ao=self.supply_air.Tdb,
                    m_da=self.m_supply,
                    Q_sen=Q_duct_sup
                )
                dT_supply_duct = supply_duct.T_ao - supply_duct.T_ai
            else:
                dT_supply_duct = Q_(0.0, 'K')
            # Determine required air state at the cooling coil exit:
            self.T_cold = self.T_supply - dT_supply_fan - dT_supply_duct
            self.cooled_air = HumidAir(Tdb=self.T_cold, W=self.supply_air.W)
        
        def determine_return_air(
            self,
            dP_fan_ret: Optional[Quantity] = None,
            eta_fan_ret: Optional[Quantity] = None,
            Q_duct_ret: Optional[Quantity] = None
        ) -> None:
            """Determines the global state of air returning from the zones and
            the state of air that is recirculated to the cooling coil, taking
            possible increase of the return air temperature into account due to
            additional heating in the return fan (if present) or due to
            additional heat taken up from the environment along the return air
            duct system.

            Parameters
            ----------
            dP_fan_ret: optional, default None
                The required pressure gain of the return fan needed to deliver
                the required mass flow rate of return air to the zones,
                determined by a pressure loss calculation of the air return duct
                system.
            eta_fan_ret: optional, default None
                The efficiency of the return fan at the required operating
                point of the fan.
            Q_duct_ret: optional, default None
                Any heat added from the environment to the return air along the
                return duct between the zones and the recirculation section of
                the VAV system.
            """
            # Determine the state of return air of each zone:
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    air_in=self.supply_air,
                    T_ao=zone.summer.zone_air.Tdb,
                    SHR=(
                        zone.summer.Q_sen /
                        (zone.summer.Q_sen + zone.summer.Q_lat)
                    )
                )
                zone.summer.return_air = p.air_out
            # Determine mass flow rate of return air from all zones:
            self.m_return = sum(
                zone.summer.m_return
                for zone in self.system.zones
            )
            # Determine global state of return air:
            h_return = sum(
                z.summer.m_return * z.summer.return_air.h
                for z in self.system.zones
            ) / self.m_return
            W_return = sum(
                z.summer.m_return * z.summer.return_air.W
                for z in self.system.zones
            ) / self.m_return
            self.return_air = HumidAir(h=h_return, W=W_return)
            # Determine volume flow rate of return air from all zones:
            self.V_return = self.m_return * self.return_air.v
            # Determine the temperature rise due to fan heating:
            if dP_fan_ret is not None:
                fan = Fan(
                    air_in=self.return_air,
                    eta_fan=(
                        eta_fan_ret
                        if eta_fan_ret is not None
                        else Q_(100, 'pct')
                    ),
                    dP_fan=dP_fan_ret
                )
                dT_return_fan = fan.air_out.Tdb - fan.air_in.Tdb
            else:
                dT_return_fan = Q_(0.0, 'K')
            # Determine the additional temperature rise due to duct heat gain:
            if Q_duct_ret is not None:
                return_duct = AirConditioningProcess(
                    T_ai=self.return_air.Tdb,
                    m_da=self.m_return,
                    Q_sen=Q_duct_ret
                )
                dT_return_duct = return_duct.T_ao - return_duct.T_ai
            else:
                dT_return_duct = Q_(0.0, 'K')
            # Determine the state of return air at the recirculation section:
            T_recirculated = (
                self.return_air.Tdb
                + dT_return_fan
                + dT_return_duct
            )
            self.recirculated_air = HumidAir(
                Tdb=T_recirculated,
                RH=self.return_air.RH
            )

        def determine_mixed_air(self) -> None:
            """Determines the state of air leaving the mixing chamber and
            entering the cooling coil of the VAV system. In the mixing chamber
            recirculated return air is mixed with outdoor ventilation air.
            """
            # Determine mass flow rate of recirculated return air:
            m_recirculated = self.m_return - self.m_vent
            # Determine mass flow rate of air that is exhausted locally in the
            # zones. This air is also taken from outdoors and must be added to
            # required mass flow rate for ventilation:
            m_exhaust = sum(
                z.summer.m_exhaust
                for z in self.system.zones
            )
            # Determine the state of air leaving the mixing chamber:
            mixing_chamber = AdiabaticMixing(
                in1=AirStream(
                    state=self.recirculated_air,
                    m_da=m_recirculated
                ),
                in2=AirStream(
                    state=self.outdoor_air,
                    m_da=self.m_vent + m_exhaust
                ),
                out=AirStream(m_da=self.m_supply)
            )
            self.mixed_air = mixing_chamber.stream_out.state

        def determine_cooling_coil(self) -> None:
            """Determines the cooling coil load."""
            self.cooling_coil = AirConditioningProcess(
                air_in=self.mixed_air,
                air_out=self.cooled_air,
                m_da=self.m_supply,
                h_w=Q_(0.0, 'J / kg')
            )

    class Winter:
        """Class for designing the VAV system for the peak winter day.

        Attributes
        ----------
        T_sup_max: default 40 °C
            Maximum allowable supply air temperature to avoid stratification.
        outdoor_air:
            State of outdoor air on the winter peak design day.
        V_vent:
            Required ventilation air volume flow rate of the building.
        m_vent:
            Required ventilation air mass flow rate of the building.
        preheat_coil:
            Refers to the sensible heating process in the preheat coil of the
            VAV system (instance of `AirConditioningProcess`).
        Q_ph_peak:
            Peak load on the preheat coil for sizing the preheat coil.
        m_supply:
            Required mass flow rate of supply air to the zones.
        V_supply:
            Required volume flow rate of supply air to the zones.
        T_supply:
            Selected supply air temperature.
        supply_air:
            State of supply air to the zones.
        m_return:
            Total mass flow rate of return air from the zones.
        V_return:
            Total volume flow rate of return air from the zones.
        return_air:
            Global state of return air from the zones.
        recirculated_air:
            State of air at the recirculation section of the VAV system.
        m_recirculated:
            Mass flow rate of air recirculated to the mixing chamber of the
            VAV system.
        mixed_air:
            State of air leaving the mixing chamber and entering the preheat
            coil.
        preheated_air:
            State of air leaving the preheat coil.
        T_cold:
            Required temperature of air leaving the cooling coil.
        cooled_air:
            State of air leaving the cooling coil.
        cooling_coil:
            Refers to the cooling and dehumidification process in the cooling
            coil of the VAV system (instance of `AirConditioningProcess`).
        Q_rh_tot:
            Total heating load of the reheat coils in the zones.
        humidifier:
            Refers to the humidification process in the humidifier of the VAV
            system (instance of `AirConditioningProcess`).
        humidified air:
            State of air at the humidifier exit = cooling coil entry.
        """
        T_sup_max = Q_(40.0, 'degC')

        def __init__(
            self,
            outdoor_air: HumidAir,
            steam: FluidState,
            V_vent: Quantity,
            system: 'VAVSystem'
        ) -> None:
            """Creates the `Winter` design stage of the VAV system.

            Parameters
            ----------
            outdoor_air:
                State of outdoor air on the winter peak design day.
            steam:
                State of steam injected in the air humidifier to humidify the
                supply air to the zones.
            V_vent:
                Required ventilation air volume flow rate of the building.
            system:
                Refers to the `VAVSystem` parent class.
            """
            self.system = system
            self.outdoor_air = outdoor_air
            self.V_vent = V_vent
            self.steam = steam

            self.m_vent = V_vent * self.outdoor_air.rho
            self.preheat_coil: Optional[AirConditioningProcess] = None
            self.Q_ph_peak: Quantity = Q_(float('nan'), 'W')
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
            self.Q_rh_tot: Quantity = Q_(float('nan'), 'W')
            self.humidifier: Optional[AirConditioningProcess] = None
            self.humidified_air: Optional[HumidAir] = None

        def determine_preheat_peak_load(self) -> None:
            """Determines the peak load of the preheat coil, being the heat rate
            needed to heat only the ventilation air mass flow rate to the
            required temperature of air leaving the cooling coil.
            """
            preheat_coil = AirConditioningProcess(
                T_ai=self.outdoor_air.Tdb,
                T_ao=self.system.summer.T_cold,
                m_da=self.m_vent
            )
            self.Q_ph_peak = max(
                Q_(0.0, 'W'),
                preheat_coil.Q_sen
            )

        def determine_m_supply(self) -> None:
            """Determines the required mass flow rate of supply air to each zone.
            Any air that is exhausted locally to the outdoor in a zone is also
            added the required mass flow rate of supply air of this zone.
            """
            for zone in self.system.zones:
                if zone.winter.Q_sen < 0.0:
                    # Zone requires heating.
                    # To determine the supply air mass flow rate, the maximum
                    # allowable supply air temperature `T_sup_max`is
                    # used:
                    p = AirConditioningProcess(
                        T_ai=self.T_sup_max,
                        T_ao=zone.winter.zone_air.Tdb,
                        Q_sen=zone.winter.Q_sen
                    )
                else:
                    # Zone requires cooling.
                    # Supply air temperature to the zone:
                    T_supply = zone.winter.zone_air.Tdb - self.system.summer.dT_sup
                    p = AirConditioningProcess(
                        T_ai=T_supply,
                        T_ao=zone.winter.zone_air.Tdb,
                        Q_sen=zone.winter.Q_sen
                    )
                zone.winter.m_supply = p.m_da + zone.winter.m_exhaust
                # The mass flow rate of supply air to a zone cannot be reduced
                # below 60 % of the peak summer flow rate:
                zone.winter.m_supply = max(
                    zone.winter.m_supply,
                    0.6 * zone.summer.m_supply
                )
            # Total mass flow rate of supply air to the zones:
            self.m_supply = sum(
                z.winter.m_supply
                for z in self.system.zones
            )

        def determine_supply_air(self) -> None:
            """Determines the required supply air temperature to each zone,
            after the supply air mass flow rates to each zone have been
            determined first. The smallest required supply air temperature to a
            zone will also be the global supply air temperature of the system.
            """
            # Determine the required state of supply air to each zone based on
            # a sensible and latent load balance of each zone:
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    air_out=zone.winter.zone_air,
                    m_da=zone.winter.m_supply,
                    Q_sen=zone.winter.Q_sen,
                    Q_lat=zone.winter.Q_lat
                )
                zone.winter.supply_air = p.air_in
            # Determine the average zone air state of all zones together:
            T_zone_avg = sum(
                zone.winter.zone_air.Tdb
                for zone in self.system.zones
            ) / len(self.system.zones)
            W_zone_avg = sum(
                zone.winter.zone_air.W
                for zone in self.system.zones
            ) / len(self.system.zones)
            zone_air_avg = HumidAir(Tdb=T_zone_avg, W=W_zone_avg)
            # Total sensible load of all zones:
            Q_sen_tot = sum(
                zone.winter.Q_sen
                for zone in self.system.zones
            )
            # Total latent load of all zones:
            Q_lat_tot = sum(
                zone.winter.Q_lat
                for zone in self.system.zones
            )
            # Global space condition line of all zones together:
            SCL = SpaceConditionLine(zone_air_avg, Q_sen_tot, Q_lat_tot)
            # Determine the system supply air temperature, being the smallest
            # supply air temperature required by the zones:
            self.T_supply = min(
                z.winter.supply_air.Tdb
                for z in self.system.zones
            )
            # Determine the humidity ratio of the system supply air based on the
            # global space condition line:
            W_supply = SCL.W_ai(self.T_supply)
            # Determine the state of the system supply air:
            self.supply_air = HumidAir(
                Tdb=self.T_supply,
                W=W_supply
            )
            # Now that the state of supply air has been determined, the volume
            # flow rate of supply air can also be determined:
            self.V_supply = self.m_supply * self.supply_air.v
            # Determine the actual state of supply air to each zone (as the
            # humidity ratio of the system supply air cannot change in the
            # reheat coils of the zones):
            for zone in self.system.zones:
                zone.winter.supply_air = HumidAir(
                    Tdb=zone.winter.supply_air.Tdb,
                    W=self.supply_air.W
                )

        def determine_return_air(
            self,
            dP_fan_ret: Optional[Quantity] = None,
            eta_fan_ret: Optional[Quantity] = None,
            Q_duct_ret: Optional[Quantity] = None
        ) -> None:
            """Determines the global state of return air from all zones and the
            state of air at the recirculation section of the system, taking any
            additional heating into account due to fan heating or duct heat gain.

            Parameters
            ----------
            dP_fan_ret: optional, default None
                The required pressure gain of the return fan needed to deliver
                the required mass flow rate of return air to the zones,
                determined by a pressure loss calculation of the air return duct
                system.
            eta_fan_ret: optional, default None
                The efficiency of the return fan at the required operating
                point of the fan.
            Q_duct_ret: optional, default None
                Any heat added from the environment to the return air along the
                return duct between the zones and the recirculation section of
                the VAV system.
            """
            # Determine the state of return air of each zone:
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    air_in=zone.winter.supply_air,
                    T_ao=zone.winter.zone_air.Tdb,
                    SHR=(
                        zone.winter.Q_sen /
                        (zone.winter.Q_sen + zone.winter.Q_lat)
                    )
                )
                zone.winter.return_air = p.air_out
            # Determine total mass flow rate of return air:
            self.m_return = sum(
                z.winter.m_return
                for z in self.system.zones
            )
            # Determine the state of system return air:
            h_return = sum(
                z.winter.return_air.h * z.winter.m_return
                for z in self.system.zones
            ) / self.m_return
            W_return = sum(
                z.winter.return_air.W * z.winter.m_return
                for z in self.system.zones
            ) / self.m_return
            self.return_air = HumidAir(h=h_return, W=W_return)
            # Determine total volume flow rate of return air:
            self.V_return = self.m_return * self.return_air.v
            # Determine the temperature rise of return air due to fan heating:
            if dP_fan_ret is not None:
                fan = Fan(
                    air_out=self.return_air,
                    eta_fan=(
                        eta_fan_ret
                        if eta_fan_ret is not None
                        else Q_(100, 'pct')
                    ),
                    dP_fan=dP_fan_ret
                )
                dT_return_fan = fan.air_out.Tdb - fan.air_in.Tdb
            else:
                dT_return_fan = Q_(0.0, 'K')
            # Determine the additional temperature rise due to duct heat gain:
            if Q_duct_ret is not None:
                return_duct = AirConditioningProcess(
                    T_ao=self.return_air.Tdb,
                    m_da=self.m_return,
                    Q_sen=Q_duct_ret
                )
                dT_return_duct = return_duct.T_ao - return_duct.T_ai
            else:
                dT_return_duct = Q_(0.0, 'K')
            # Determine the state of air at the recirculation section of the
            # system:
            T_recirculated = (
                self.return_air.Tdb
                + dT_return_fan
                + dT_return_duct
            )
            self.recirculated_air = HumidAir(
                Tdb=T_recirculated,
                RH=self.return_air.RH
            )

        def determine_mixed_air(self) -> None:
            """Determines the state of air leaving the mixing chamber of the
            VAV system.
            """
            # Mass flow rate of air that is recirculated:
            self.m_recirculated = self.m_return - self.m_vent
            # Total mass flow rate of air that is locally exhausted in the zones.
            # This mass flow rate must also be supplied from outdoors and
            # therefore added to the ventilation flow rate:
            m_exhaust = sum(z.winter.m_exhaust for z in self.system.zones)
            # Adiabatic mixing of outdoor ventilation air and recirculated air:
            mixing_chamber = AdiabaticMixing(
                in1=AirStream(
                    state=self.recirculated_air,
                    m_da=self.m_recirculated
                ),
                in2=AirStream(
                    state=self.outdoor_air,
                    m_da=self.m_vent + m_exhaust
                ),
                out=AirStream(m_da=self.m_supply)
            )
            self.mixed_air = mixing_chamber.stream_out.state

        def determine_cooled_air(
            self,
            dP_fan_sup: Optional[Quantity] = None,
            eta_fan_sup: Optional[Quantity] = None,
            Q_duct_sup: Optional[Quantity] = None
        ) -> None:
            """Determines the required state of air leaving the system's cooling
            coil, based on the required state of the supply air to the zones and
            possibly taking any fan and/or duct heat gain into account.

            Parameters
            ----------
            dP_fan_sup: optional, default None
                The required pressure gain of the supply fan needed to deliver
                the required mass flow rate of supply air to the zones,
                determined by a pressure loss calculation of the air supply duct
                system.
            eta_fan_sup: optional, default None
                The efficiency of the supply fan at the required operating
                point of the fan.
            Q_duct_sup: optional, default None
                Any heat added from the environment to the supply air along the
                supply duct between the cooling coil and the zones.
            """
            # Determine the temperature rise of supply air due to fan heating:
            if dP_fan_sup is not None:
                fan = Fan(
                    air_out=self.supply_air,
                    eta_fan=(
                        eta_fan_sup
                        if eta_fan_sup is not None
                        else Q_(100, 'pct')
                    ),
                    dP_fan=dP_fan_sup
                )
                dT_supply_fan = fan.air_out.Tdb - fan.air_in.Tdb
            else:
                dT_supply_fan = Q_(0.0, 'K')
            # Determine the additional temperature rise due to duct heat gain:
            if Q_duct_sup is not None:
                supply_duct = AirConditioningProcess(
                    T_ao=self.supply_air.Tdb,
                    m_da=self.m_supply,
                    Q_sen=Q_duct_sup
                )
                dT_supply_duct = supply_duct.T_ao - supply_duct.T_ai
            else:
                dT_supply_duct = Q_(0.0, 'K')
            # Determine the required state of air at the cooling coil exit:
            self.T_cold = (
                self.T_supply
                - dT_supply_fan
                - dT_supply_duct
            )
            self.cooled_air = HumidAir(
                Tdb=self.T_cold,
                W=self.supply_air.W
            )

        def determine_preheated_air(self) -> None:
            """Determines the required state of air leaving the preheat coil.
            If the required temperature at the cooling coil exit is greater
            than the mixed air temperature at the preheat coil entry, the air is
            heated by the preheat coil to the required air temperature at the
            cooling coil exit.
            """
            if self.cooled_air.Tdb > self.mixed_air.Tdb:
                self.preheat_coil = AirConditioningProcess(
                    air_in=self.mixed_air,
                    air_out=HumidAir(
                        Tdb=self.cooled_air.Tdb,
                        W=self.mixed_air.W
                    ),
                    m_da=self.m_supply
                )
                self.preheated_air = self.preheat_coil.air_out
            else:
                self.preheat_coil = None
                self.preheated_air = self.mixed_air

        def determine_humidified_air(self) -> None:
            """Determines the required state of air leaving the humidifier.
            If the required humidity ratio at the cooling coil exit is greater
            than the humidity ratio of mixed air, the air is humidified by the
            humidifier to the required state at the cooling coil exit.
            """
            if self.cooled_air.W > self.mixed_air.W:
                self.humidifier = AirConditioningProcess(
                    air_in=self.preheated_air,
                    air_out=HumidAir(
                        Tdb=self.preheated_air.Tdb,
                        W=self.cooled_air.W
                    ),
                    m_da=self.m_supply,
                    Q=Q_(0.0, 'W'),
                    h_w=self.steam.h
                )
                self.humidified_air = self.humidifier.air_out
            else:
                self.humidifier = None
                self.humidified_air = self.mixed_air

        def determine_cooling_coil(self) -> None:
            """Determines the cooling coil load."""
            self.cooling_coil = AirConditioningProcess(
                air_in=self.humidified_air,
                air_out=self.cooled_air,
                m_da=self.m_supply,
                h_w=Q_(0.0, 'J / kg')
            )

        def determine_reheat_coils(self) -> None:
            """Determines the load of the reheat coil in each zone and the
            total load on all reheat coils.
            """
            for zone in self.system.zones:
                zone.reheat_coil = AirConditioningProcess(
                    T_ai=self.T_supply,
                    T_ao=zone.winter.supply_air.Tdb,
                    m_da=zone.winter.m_supply
                )
            self.Q_rh_tot = sum(
                z.reheat_coil.Q_sen
                for z in self.system.zones
            )

    def __init__(
        self,
        zones: List[Zone],
        outdoor_air_summer: HumidAir,
        outdoor_air_winter: HumidAir,
        steam_winter: FluidState,
        V_vent: Quantity
    ) -> None:
        """Creates a `VAVSystem` instance.

        Parameters
        ----------
        zones:
            list of `Zone` objects that represent the temperature zones of the
            building served by the VAV system.
        outdoor_air_summer:
            State of outdoor air on the peak summer design day.
        outdoor_air_winter:
            State of outdoor air on the peak winter design day.
        steam_winter:
            State of steam for humidification of supply air.
        V_vent:
            Minimum required outdoor air volume flow rate for ventilating the
            building.
        """
        self.zones = zones
        self.summer = VAVSystem.Summer(outdoor_air_summer, V_vent, self)
        self.winter = VAVSystem.Winter(outdoor_air_winter, steam_winter, V_vent, self)

    def design_summer(self, **kwargs) -> Dict[str, Quantity]:
        """Designs the VAV system for peak summer conditions.

        Optional keyword arguments
        --------------------------
        dT_sup:
            Allowable temperature difference between zone air and supply air
            to the zones.
        eta_fan_sup:
            Efficiency of the supply fan.
        dP_fan_sup:
            Feed pressure the supply fan needs to create to produce the required
            flow rate of supply air.
        Q_duct_sup:
            Heat transfer rate from the environment to the supply air along the
            air supply duct system.
        eta_fan_ret:
            Efficiency of the return fan.
        dP_fan_ret:
            Feed pressure the return fan needs to create to produce the required
            flow rate of return air.
        Q_duct_ret:
            Heat transfer rate from the environment to the return air along the
            air return duct system.

        Returns
        -------
        Dictionary with calculation results:
        - 'Q_cc_tot'
        - 'Q_cc_sen'
        - 'Q_cc_lat'
        - 'V_sup'
        - 'V_ret'
        - 'T_sup'
        - 'T_ret'
        """
        dT_supply = kwargs.get('dT_sup')

        eta_supply_fan = kwargs.get('eta_fan_sup')
        dP_supply_fan = kwargs.get('dP_fan_sup')
        Q_supply_duct = kwargs.get('Q_duct_sup')

        eta_return_fan = kwargs.get('eta_fan_ret')
        dP_return_fan = kwargs.get('dP_fan_ret')
        Q_return_duct = kwargs.get('Q_duct_ret')

        if dT_supply is not None:
            self.summer.dT_sup = dT_supply

        self.summer.determine_supply_air()
        self.summer.determine_m_supply()
        self.summer.determine_cooled_air(
            dP_supply_fan, eta_supply_fan,
            Q_supply_duct
        )
        self.summer.determine_return_air(
            dP_return_fan, eta_return_fan,
            Q_return_duct
        )
        self.summer.determine_mixed_air()
        self.summer.determine_cooling_coil()

        results = {
            'Q_cc_tot': self.summer.cooling_coil.Q,
            'Q_cc_sen': self.summer.cooling_coil.Q_sen,
            'Q_cc_lat': self.summer.cooling_coil.Q_lat,
            'V_sup': self.summer.V_supply,
            'V_ret': self.summer.V_return,
            'T_sup': self.summer.supply_air.Tdb,
            'T_ret': self.summer.return_air.Tdb
        }
        return results

    def design_winter(self, **kwargs) -> Dict[str, Quantity]:
        """Designs the VAV system for peak winter conditions.

        Optional keyword arguments
        --------------------------
        T_sup_max:
            Permissible maximum temperature of supply air to a zone.
        eta_fan_sup:
            Efficiency of the supply fan.
        dP_fan_sup:
            Feed pressure the supply fan needs to create to produce the required
            flow rate of supply air.
        Q_duct_sup:
            Heat transfer rate from the environment to the supply air along the
            air supply duct system.
        eta_fan_ret:
            Efficiency of the return fan.
        dP_fan_ret:
            Feed pressure the return fan needs to create to produce the required
            flow rate of return air.
        Q_duct_ret:
            Heat transfer rate from the environment to the return air along the
            air return duct system.

        Returns
        -------
        Dictionary with calculation results:
        - 'Q_ph_peak'
        - 'Q_ph'
        - 'm_steam'
        - 'Q_cc_tot'
        - 'Q_cc_sen'
        - 'Q_cc_lat'
        - 'Q_rh_tot'
        - 'V_sup'
        - 'V_ret'
        - 'T_sup'
        - 'T_ret'
        """
        T_sup_max = kwargs.get('T_sup_max')

        eta_supply_fan = kwargs.get('eta_fan_sup')
        dP_supply_fan = kwargs.get('dP_fan_sup')
        Q_supply_duct = kwargs.get('Q_duct_sup')

        eta_return_fan = kwargs.get('eta_fan_ret')
        dP_return_fan = kwargs.get('dP_fan_ret')
        Q_return_duct = kwargs.get('Q_duct_ret')

        if T_sup_max is not None: self.winter.T_sup_max = T_sup_max

        self.winter.determine_preheat_peak_load()
        self.winter.determine_m_supply()
        self.winter.determine_supply_air()
        self.winter.determine_return_air(
            dP_return_fan, eta_return_fan,
            Q_return_duct
        )
        self.winter.determine_mixed_air()
        self.winter.determine_cooled_air(
            dP_supply_fan, eta_supply_fan,
            Q_supply_duct
        )
        self.winter.determine_preheated_air()
        self.winter.determine_humidified_air()
        self.winter.determine_cooling_coil()
        self.winter.determine_reheat_coils()

        results = {
            'Q_ph_peak': self.winter.Q_ph_peak,
            'Q_ph': (
                self.winter.preheat_coil.Q_sen
                if self.winter.preheat_coil is not None
                else Q_(0.0, 'kW')
            ),
            'm_steam': (
                self.winter.humidifier.m_w
                if self.winter.humidifier is not None
                else Q_(0.0, 'kg / s')
            ),
            'Q_cc_tot': self.winter.cooling_coil.Q,
            'Q_cc_sen': self.winter.cooling_coil.Q_sen,
            'Q_cc_lat': self.winter.cooling_coil.Q_lat,
            'Q_rh_tot': self.winter.Q_rh_tot,
            'V_sup': self.winter.V_supply,
            'V_ret': self.winter.V_return,
            'T_sup': self.winter.supply_air.Tdb,
            'T_ret': self.winter.return_air.Tdb
        }
        return results
