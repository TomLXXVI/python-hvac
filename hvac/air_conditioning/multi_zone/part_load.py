from typing import Optional, List, Dict
from hvac import Quantity
from ..core import AirConditioningProcess, Fan, AirStream, AdiabaticMixing
from hvac.fluids import HumidAir
from .design import Zone


Q_ = Quantity


class VAVSystem:
    """Class that implements the analysis of a VAV system under part-load
    conditions. Only part-load analysis during summer (cooling) has been
    implemented.
    """

    class Summer:
        """Represents VAV system operation under summer part-load conditions.

        Attributes
        ----------
        T_supply_des
        outdoor_air
        m_vent
        m_supply
        V_supply
        supply_air
        m_return
        V_return
        return_air
        recirculated_air
        mixed_air
        cooled_air
        cooling_coil
        Q_reheat
        """
        def __init__(
            self,
            system: 'VAVSystem',
            T_supply_des: Quantity,
            outdoor_air: HumidAir,
            V_vent: Quantity
        ) -> None:
            """Creates the `Summer` part of the `VAVSystem` object.

            Parameters
            ----------
            system:
                Reference to the parent `VAVSystem` object.
            T_supply_des:
                Supply air temperature to the zones (design value).
            outdoor_air:
                State of outdoor air under part-load conditions.
            V_vent:
                Minimum required outdoor air volume flow rate for ventilation
                of the building (design value).
            """
            self.system = system
            self.T_supply_des = T_supply_des
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

        def determine_m_supply(self) -> None:
            """Determines the mass flow rate of supply air to each zone
             based on the fixed system supply air temperature and a sensible
             heat balance of each zone. Then, determines the total mass flow rate
             the supply fan needs to transport to all zones.
            """
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    T_ai=self.T_supply_des,  # fixed system supply air temperature
                    T_ao=zone.summer.zone_air.Tdb,
                    Q_sen=zone.summer.Q_sen
                )
                # Supply air mass flow rate to a zone at part-load cannot be
                # reduced below 60 % of the full-load design value:
                zone.summer.m_supply = max(
                    0.6 * zone.summer.m_supply_des,  # --> full-load design value
                    p.m_da + zone.summer.m_exhaust
                )
            # Total mass flow rate of system supply air:
            self.m_supply = sum(
                zone.summer.m_supply
                for zone in self.system.zones
            )

        def determine_supply_air(self) -> None:
            """Determines the required state of supply air to each zone to offset
            both the sensible and latent load of the zone, so that the desired
            state of zone air is maintained. Then, determines the global state
            of the system's supply air by keeping its temperature fixed to the
            design value given by the user and taking the average of the
            required humidity ratios in the zones.
            """
            # Determine the state of supply air that each zone would require:
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    air_out=zone.summer.zone_air,
                    Q_sen=zone.summer.Q_sen,
                    Q_lat=zone.summer.Q_lat,
                    m_da=zone.summer.m_supply
                )
                zone.summer.supply_air = p.air_in
            # Determine the required humidity ratio of the system air supplied to
            # all zones by averaging the zone air humidities:
            W_supply_avg = sum(
                z.summer.supply_air.W * z.summer.m_supply
                for z in self.system.zones
            ) / self.m_supply
            # Determine the required state of system supply air to all zones:
            self.supply_air = HumidAir(
                Tdb=self.T_supply_des,  # --> fixed system value
                W=W_supply_avg
            )
            # Determine the volume flow rate of supply air to all the zones:
            self.V_supply = self.m_supply * self.supply_air.v

        def determine_cooled_air(
            self,
            dP_fan_sup: Optional[Quantity] = None,
            eta_fan_sup: Optional[Quantity] = None,
            Q_duct_sup: Optional[Quantity] = None
        ) -> None:
            """Determines the required state of air leaving the cooling coil,
            taking possible fan and duct heat gains into account between the
            cooling coil and the zones.
            """
            # Determine the temperature rise of the supply air due to fan
            # heating:
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
            # Determine the temperature rise of the supply air due to duct heat
            # gain:
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
            T_cold = self.T_supply_des - dT_supply_fan - dT_supply_duct
            self.cooled_air = HumidAir(
                Tdb=T_cold,
                W=self.supply_air.W
                #  For simplicity, it is assumed that the operating conditions
                #  of the cooling coil will be such that the calculated
                #  required state of the air leaving the cooling coil
                #  is always achieved.
            )

        def determine_return_air(
            self,
            dP_fan_ret: Optional[Quantity] = None,
            eta_fan_ret: Optional[Quantity] = None,
            Q_duct_ret: Optional[Quantity] = None
        ) -> None:
            """Determines the actual state of return air from each zone, being
            supplied with air of which the temperature is controlled by the
            zone's reheat coil. The humidity level of the air supplied to a zone
            cannot be controlled.
            """
            # For each zone, determine the actual state of zone air which is
            # also the air being returned to the air handling unit:
            for zone in self.system.zones:
                p = AirConditioningProcess(
                    air_in=HumidAir(
                        Tdb=zone.summer.supply_air.Tdb,  # downstream of reheat-coil
                        W=self.supply_air.W  # remains fixed
                    ),
                    T_ao=zone.summer.zone_air.Tdb,
                    SHR=(
                        zone.summer.Q_sen
                        / (zone.summer.Q_sen + zone.summer.Q_lat)
                    )
                )
                zone.summer.return_air = p.air_out
            # Determine the total mass flow rate of air returned to the air
            # handling unit:
            self.m_return = sum(
                zone.summer.m_return
                for zone in self.system.zones
            )
            # Determine the state of the system return air:
            h_return = sum(
                z.summer.m_return * z.summer.return_air.h
                for z in self.system.zones
            ) / self.m_return
            W_return = sum(
                z.summer.m_return * z.summer.return_air.W
                for z in self.system.zones
            ) / self.m_return
            self.return_air = HumidAir(h=h_return, W=W_return)
            # Determine the volume flow rate of return air:
            self.V_return = self.m_return * self.return_air.v
            # Determine temperature rise of return air due to fan heating:
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
            # Determine temperature rise of return air due to duct heat gain:
            if Q_duct_ret is not None:
                return_duct = AirConditioningProcess(
                    T_ai=self.return_air.Tdb,
                    m_da=self.m_return,
                    Q_sen=Q_duct_ret
                )
                dT_return_duct = return_duct.T_ao - return_duct.T_ai
            else:
                dT_return_duct = Q_(0.0, 'K')
            # Determine the state of air at the recirculation section:
            T_recirculated = (
                self.return_air.Tdb
                + dT_return_fan
                + dT_return_duct
            )
            self.recirculated_air = HumidAir(
                Tdb=T_recirculated,
                RH=self.return_air.RH
            )

        def determine_mixed_air(self):
            """Determines the state of air leaving the mixing chamber."""
            m_recirculated = self.m_return - self.m_vent
            m_exhaust = sum(
                z.summer.m_exhaust
                for z in self.system.zones
            )
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

        def determine_cooling_coil(self):
            """Determines the cooling coil load."""
            self.cooling_coil = AirConditioningProcess(
                air_in=self.mixed_air,
                air_out=self.cooled_air,
                m_da=self.m_supply,
                h_w=Q_(0.0, 'J / kg')
            )

        def determine_reheat_coils(self):
            """Determines the load on the reheat-coils of the zones."""
            for zone in self.system.zones:
                zone.reheat_coil = AirConditioningProcess(
                    T_ai=self.T_supply_des,
                    T_ao=zone.summer.supply_air.Tdb,
                    m_da=zone.summer.m_supply
                )
            self.Q_reheat = sum(
                z.reheat_coil.Q_sen
                for z in self.system.zones
            )

    def __init__(
        self,
        zones: List[Zone],
        T_supply_des: Quantity,
        outdoor_air_summer: HumidAir,
        V_vent: Quantity
    ) -> None:
        """Creates a `VAVSystem` instance.

        Parameters
        ----------
        zones:
            List of `Zone` objects that represent the temperature zones of the
            building served by the VAV system.
        T_supply_des:
            The design value of the supply air temperature to the zones.
        outdoor_air_summer:
            State of outdoor air under part-load conditions.
        V_vent:
            Minimum required outdoor air volume flow rate for ventilating the
            building (design value).
        """
        self.zones = zones
        self.summer = VAVSystem.Summer(self, T_supply_des, outdoor_air_summer, V_vent)

    def part_load_summer(self, **kwargs) -> Dict[str, Quantity]:
        """Determine VAV system operation under summer part-load conditions.

        Optional keyword arguments
        --------------------------
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
        - 'Q_rh_tot'
        - 'V_sup'
        - 'V_ret'
        - 'T_sup'
        - 'T_ret'
        """
        eta_supply_fan = kwargs.get('eta_fan_sup')
        dP_supply_fan = kwargs.get('dP_fan_sup')
        Q_supply_duct = kwargs.get('Q_duct_sup')
        eta_return_fan = kwargs.get('eta_fan_ret')
        dP_return_fan = kwargs.get('dP_fan_ret')
        Q_return_duct = kwargs.get('Q_duct_ret')

        self.summer.determine_m_supply()
        self.summer.determine_supply_air()
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
        self.summer.determine_reheat_coils()

        results = {
            'Q_cc_tot': self.summer.cooling_coil.Q,
            'Q_cc_sen': self.summer.cooling_coil.Q_sen,
            'Q_cc_lat': self.summer.cooling_coil.Q_lat,
            'Q_rh_tot': self.summer.Q_reheat,
            'V_sup': self.summer.V_supply,
            'V_ret': self.summer.V_return,
            'T_sup': self.summer.supply_air.Tdb,
            'T_ret': self.summer.return_air.Tdb
        }
        return results
