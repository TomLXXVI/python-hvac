"""
Steady-state part-load analysis of a single-zone VAV air-cooling system equipped
with a DX air-cooling coil and optionally an air-heating coil.
"""
from collections.abc import Iterator
import warnings
from dataclasses import dataclass
import pandas as pd
from hvac import Quantity
from hvac.fluids import HumidAir, Fluid, CoolPropWarning, CP_HUMID_AIR
import hvac.heat_transfer.heat_exchanger.fin_tube.core as core
from hvac.air_conditioning import AirConditioningProcess, AirStream, AdiabaticMixing
from hvac.heat_transfer.heat_exchanger.fin_tube.air_evaporator.dx_coil import DXAirCoolingCoil
from hvac.charts.psychrometric_chart import PsychrometricChart, StatePoint
from .cooling_design import Output


warnings.filterwarnings('ignore', category=CoolPropWarning)


Q_ = Quantity
HexCore = core.plain_fin_tube.PlainFinTubeHeatExchangerCore


# Define standard (NTP) air:
Air = Fluid('Air')
standard_air = Air(
    T=Q_(20, 'degC'),
    P=Q_(101_325, 'Pa')
)


@dataclass
class DesignData:
    """Groups the design data of the single-zone VAV air-cooling system.

    Attributes
    ----------
    Q_zone_sen: Quantity
        Sensible zone load at design conditions.
    Q_zone_lat: Quantity
        Latent zone load at design conditions.
    V_dot_vent_ntp:
        Minimum required volume flow rate of outdoor ventilation air referred
        to NTP-standard air (dry air at 20 °C and 101.325 kPa).
    V_dot_supply_ntp:
        Required supply air volume flow rate at design conditions referred
        to NTP-standard air.
    outdoor_air:
        State of outdoor air at design conditions.
    zone_air:
        Required state of zone air.
    supply_air:
        Required state of supply air at design conditions.
    ----------------------------------------------------------------------------
    m_dot_vent:
        Outdoor ventilation air mass flow rate. Internally calculated.
    m_dot_supply:
        Supply air mass flow rate. Internally calculated.
    m_dot_recir:
        Recirculation air mass flow rate. Internally calculated.
    V_dot_recir_ntp:
        Recirculation air volume flow rate. Internally calculated.
    """
    Q_zone_sen: Quantity
    Q_zone_lat: Quantity
    V_dot_vent_ntp: Quantity
    V_dot_supply_ntp: Quantity
    outdoor_air: HumidAir
    zone_air: HumidAir
    supply_air: HumidAir

    def __post_init__(self):
        self.m_dot_vent = self.V_dot_vent_ntp * standard_air.rho
        self.m_dot_supply = self.V_dot_supply_ntp * standard_air.rho
        self.m_dot_recir = self.m_dot_supply - self.m_dot_vent
        self.V_dot_recir_ntp = self.m_dot_recir / standard_air.rho


class VAVSingleZoneAirCoolingSystem:
    """Model for analyzing steady-state part-load operation of a single-zone
    VAV air-cooling system equipped with a DX air-cooling coil and optionally an
    air-heating coil.
    """
    def __init__(
        self,
        design_data: DesignData,
        dx_coil_hex_core: HexCore,
        dx_coil_Rfg: Fluid,
        dx_coil_T_evp: Quantity,
        dx_coil_x_rfg: Quantity = Q_(0.25, 'frac'),
        T_zone: Quantity | None = None,
        T_cool: Quantity | None = None,
        heating_coil_present: bool = True,
        units: dict[str, tuple[str, int]] | None = None
    ) -> None:
        """
        Creates a `VAVSingleZoneAirCoolingSystem` object.

        To analyze the steady-state operation of the single-zone VAV air cooling
        system, the refrigerant-side operating parameters of the DX air-cooling
        coil are needed in order to be able to determine the state of cooled air
        at the air-cooling coil outlet.

        Parameters
        ----------
        design_data:
            `DesignData` instance that groups the design data of the AC system
            (see class `DesignData` in this module).
        dx_coil_hex_core:
            Object that represents the heat exchanger core of the DX-coil.
        dx_coil_Rfg:
            Type of refrigerant used in the DX air-cooling coil.
        dx_coil_T_evp:
            Evaporation temperature of refrigerant in the DX air-cooling coil.
        dx_coil_x_rfg: optional
            Vapor quality of refrigerant at the inlet of the DX air-cooling
            coil.
        T_zone: optional
            Setpoint of zone air temperature. If None, the zone air temperature
            in `design_data` will be used.
        T_cool: optional
            Setpoint of air temperature at cooling coil outlet. If None, the
            supply air temperature in `design_data` will be used.
        heating_coil_present: default True
            Indicates if a heating coil is present downstream of the cooling
            coil or not.
        units: optional
            Units to be used for displaying the results. Keys are:
            'm_dot' for mass flow rate, 'V_dot' for volume flow rate, 'Q_dot'
            for heat transfer rate, 'T' for temperature, 'W' for absolute
            humidity, 'RH' for relative humidity, and 'SHR' for the sensible heat
            ratio. Values are 2-tuples: the first element is the unit
            (a string), the second element is the number of decimals to be
            displayed (an int).
        """
        self.design_data = design_data
        self.hex_core = dx_coil_hex_core
        self.T_evp = dx_coil_T_evp
        self.P_evp = dx_coil_Rfg(T=self.T_evp, x=Q_(1, 'frac')).P
        self.rfg_in = dx_coil_Rfg(P=self.P_evp, x=dx_coil_x_rfg)
        self.T_zone = T_zone if T_zone is not None else self.design_data.zone_air.Tdb
        self.T_cool = T_cool if T_cool is not None else self.design_data.supply_air.Tdb
        self.heating_coil_present = heating_coil_present
        self.units = units or {}

        self.outdoor_air: HumidAir | None = None
        self.rel_Q_zone_sen: Quantity | None = None
        self.rel_Q_zone_lat: Quantity | None = None

        self.m_dot_supply: Quantity | None = None
        self.m_dot_recir: Quantity | None = None
        self.m_dot_vent: Quantity | None = None
        self.cooled_air: HumidAir | None = None
        self.supply_air: HumidAir | None = None
        self.return_air: HumidAir | None = None
        self.mixed_air: HumidAir | None = None
        self.Q_dot_cc: Quantity | None = None
        self.SHR_cc: Quantity | None = None
        self.Q_dot_hc: Quantity | None = None

    def analyze(
        self,
        outdoor_air: HumidAir,
        rel_Q_zone_sen: Quantity,
        rel_Q_zone_lat: Quantity,
        i_max: int = 50,
        W_tol: Quantity = Q_(0.01, 'g / kg')
    ) -> Output:
        """For the operating parameters given, returns the needed mass flow
        rate of supply air, the mass flow rate of recirculation air, the mass
        flow rate of outdoor ventilation air, the states of supply air, return
        air and mixed air in the AC system, and the cooling coil load.

        Parameters
        ----------
        outdoor_air:
            Outdoor air state.
        rel_Q_zone_sen:
            Fraction of the sensible zone load at design conditions (either
            a fraction between 0 and 1, or a percentage between 0 and 100 %).
        rel_Q_zone_lat:
            Fraction of the latent zone load at design conditions (either
            a fraction between 0 and 1, or a percentage between 0 and 100 %).
        i_max: optional, default 50
            Maximum number of iterations.
        W_tol: optional, default 0.01 g/kg
            Tolerance for the supply air absolute humidity.

        Returns
        -------
        `Output` object (see class `Output` in module `cooling_design`).
        """
        self.outdoor_air = outdoor_air
        self.rel_Q_zone_sen = rel_Q_zone_sen
        self.rel_Q_zone_lat = rel_Q_zone_lat

        # The minimum supply air mass flow rate must not be less than the
        # minimum ventilation air requirement or not be less than 60 % of the
        # design value to ensure proper mixing of supply air with zone air:
        m_dot_supply_min = max(
            self.design_data.m_dot_vent,
            0.60 * self.design_data.m_dot_supply
        )

        # Initial guess of the supply air to the zone:
        supply_air = HumidAir(
            Tdb=self.T_cool,
            RH=self.design_data.supply_air.RH
        )

        for i in range(i_max):
            # The zone thermostat determines the setpoint for the supply air
            # mass flow rate to the zone.
            zone = AirConditioningProcess(
                air_in=supply_air,
                T_ao=self.T_zone,
                Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
                Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
            )
            # However, the supply air mass flow rate must not become less than
            # the minimum allowable limit.
            if zone.m_da < m_dot_supply_min:
                zone = AirConditioningProcess(
                    air_in=supply_air,
                    m_da=m_dot_supply_min,
                    Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
                    Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
                )
            # Should the supply air mass flow rate now become greater than the
            # mass flow rate that is required to maintain the setpoint zone air
            # temperature, the zone air temperature will drop below the setpoint
            # value. In that case, the supply air temperature will need to be
            # raised in a heating coil located downstream of the cooling coil:
            if self.heating_coil_present and zone.air_out.Tdb < self.T_zone:
                zone = AirConditioningProcess(
                    W_ai=supply_air.W,
                    T_ao=self.T_zone,
                    m_da=zone.m_da,
                    Q_sen=self.rel_Q_zone_sen * self.design_data.Q_zone_sen,
                    Q_lat=self.rel_Q_zone_lat * self.design_data.Q_zone_lat
                )
            supply_air = zone.air_in
            return_air = zone.air_out
            m_dot_supply = zone.m_da

            # Outdoor air control (economizer function)
            if self.outdoor_air.Tdb >= return_air.Tdb:
                # If outdoor air temperature is higher than return air temperature:
                # outdoor ventilation air flow rate is restricted to the minimum
                # ventilation requirement (i.e., the design flow rate)
                # --> OA-damper in its minimum open position (NC) and RA-damper
                # in its maximum open position (NO)
                m_dot_vent = self.design_data.m_dot_vent
                m_dot_recir = m_dot_supply - m_dot_vent
                mixing_chamber = AdiabaticMixing(
                    in1=AirStream(
                        state=return_air,
                        m_da=m_dot_recir
                    ),
                    in2=AirStream(
                        state=self.outdoor_air,
                        m_da=self.design_data.m_dot_vent
                    ),
                    out=AirStream(m_da=m_dot_supply)
                )
                mixed_air = mixing_chamber.stream_out.state
            elif self.outdoor_air.Tdb < self.T_cool:
                # If outdoor air temperature is less than the setpoint supply
                # air temperature: outdoor ventilation air is mixed with return
                # air to get mixed air with a temperature equal to the setpoint
                # supply air temperature --> air mixing controller controls
                # the position of the OA- and RA-damper.
                mixing_chamber = AdiabaticMixing(
                    in1=AirStream(
                        state=return_air
                    ),
                    in2=AirStream(
                        state=self.outdoor_air
                    ),
                    out=AirStream(Tdb=self.T_cool, m_da=m_dot_supply)
                )
                m_dot_vent = mixing_chamber.stream_in2.m_da
                m_dot_recir = mixing_chamber.stream_in1.m_da
                mixed_air = mixing_chamber.stream_out.state
            else:
                # If outdoor air temperature is between the setpoint supply air
                # temperature and the return air temperature: mixing
                # cannot attain the setpoint supply air temperature, because
                # mixing outdoor air with return air would inevitably raise the
                # mixed air temperature above the outdoor air temperature.
                # It will cost more cooling energy to cool this warmer mixed air
                # to the setpoint supply air temperature. Therefore, no mixing
                # takes place --> OA-damper in its maximum open position
                # and RA-damper in its minimum open position.
                m_dot_vent = m_dot_supply
                m_dot_recir = Q_(0.0, 'kg / s')
                mixed_air = self.outdoor_air

            # The cooling coil thermostat controls the cooling coil outlet air
            # temperature.
            if mixed_air.Tdb > self.T_cool:
                # Mechanical cooling is activated.
                cooling_coil = DXAirCoolingCoil(
                    m_dot_air=m_dot_supply,
                    air_in=mixed_air,
                    T_air_out=self.T_cool,  # fixed by cooling coil thermostat
                    rfg_in=self.rfg_in,
                    hex_core=self.hex_core
                )
                cooled_air = cooling_coil.air_out
            else:
                # No mechanical cooling is needed.
                cooling_coil = None
                cooled_air = mixed_air

            # If the cooling coil outlet air temperature is less than the
            # required supply air temperature to maintain the setpoint zone air
            # temperature, heating will be needed:
            if self.heating_coil_present and cooled_air.Tdb < supply_air.Tdb:
                heating_coil = AirConditioningProcess(
                    air_in=cooled_air,
                    air_out=supply_air,
                    m_da=m_dot_supply
                )
                supply_air_new = HumidAir(Tdb=supply_air.Tdb, W=cooled_air.W)
            else:
                heating_coil = None
                supply_air_new = cooled_air

            dev = supply_air_new.W.to('g / kg') - supply_air.W.to('g / kg')
            if abs(dev.m) < W_tol.to('g / kg'):
                self.m_dot_supply = m_dot_supply
                self.m_dot_vent = m_dot_vent
                self.m_dot_recir = m_dot_recir
                self.cooled_air = cooled_air
                self.supply_air = supply_air_new
                self.return_air = return_air
                self.mixed_air = mixed_air
                if cooling_coil is None:
                    self.Q_dot_cc = Q_(0.0, 'W')
                    self.SHR_cc = Q_(float('nan'), 'frac')
                else:
                    self.Q_dot_cc = cooling_coil.Q
                    Q_dot_cc_sen = (
                        self.m_dot_supply * CP_HUMID_AIR
                        * (self.mixed_air.Tdb - self.cooled_air.Tdb)
                    )
                    self.SHR_cc = Q_dot_cc_sen.to('W') / self.Q_dot_cc.to('W')
                if heating_coil is None:
                    self.Q_dot_hc = Q_(0.0, 'W')
                else:
                    self.Q_dot_hc = heating_coil.Q_sen
                return Output(
                    self.m_dot_supply, self.m_dot_vent, self.m_dot_recir,
                    self.mixed_air, self.cooled_air, self.supply_air, self.return_air,
                    self.Q_dot_cc, self.SHR_cc, self.Q_dot_hc, self.units
                )
            supply_air = supply_air_new

    @property
    def psychrometric_chart(self) -> PsychrometricChart:
        """Returns a psychrometric chart with the AC processes drawn on it."""
        psy_chart = PsychrometricChart()
        psy_chart.plot_process(
            name='air mixing',
            start_point=StatePoint(
                self.outdoor_air.Tdb,
                self.outdoor_air.W
            ),
            end_point=StatePoint(
                self.return_air.Tdb,
                self.return_air.W
            ),
            mix_point=StatePoint(
                self.mixed_air.Tdb,
                self.mixed_air.W
            )
        )
        psy_chart.plot_process(
            name='air cooling',
            start_point=StatePoint(
                self.mixed_air.Tdb,
                self.mixed_air.W
            ),
            end_point=StatePoint(
                self.cooled_air.Tdb,
                self.cooled_air.W
            )
        )
        psy_chart.plot_process(
            name='air heating',
            start_point=StatePoint(
                self.cooled_air.Tdb,
                self.cooled_air.W
            ),
            end_point=StatePoint(
                self.supply_air.Tdb,
                self.supply_air.W
            )
        )
        psy_chart.plot_process(
            name='zone heating and humidification',
            start_point=StatePoint(
                self.supply_air.Tdb,
                self.supply_air.W
            ),
            end_point=StatePoint(
                self.return_air.Tdb,
                self.return_air.W
            )
        )
        return psy_chart


class CoolingSimData:
    """Helper class to prepare the simulation data needed to analyze the
    `VAVSingleZoneAirCoolingSystem` model on a daily basis.
    """
    def __init__(
        self,
        T_outdoor_db_rng: list[Quantity],
        T_outdoor_wb_rng: list[Quantity],
        df_Q_zone: pd.DataFrame,
        design_data: DesignData
    ) -> None:
        """Creates an instance of class `CoolingSimData`.

        Parameters
        ----------
        T_outdoor_db_rng:
            List of hourly dry-bulb outdoor air temperatures ordered from 0 to 23 h.
        T_outdoor_wb_rng:
            List of hourly wet-bulb outdoor air temperatures ordered from 0 to 23 h.
        df_Q_zone:
            Pandas-DataFrame object retrieved from a `Building`, a `BuildingEntity`,
            a `VentilationZone`, or a `Space` object in subpackage
            `hvac.cooling_load_calc` that contains the hourly sensible and latent
            cooling loads.
        design_data:
            `DesignData` object that contains data coming from the design calculations
            of the single-zone VAV air-cooling system, including the sensible and
            latent design cooling load.
        """
        self.T_outdoor_db_rng = T_outdoor_db_rng
        self.T_outdoor_wb_rng = T_outdoor_wb_rng
        self.df_dQ_zone = df_Q_zone
        self.design_data = design_data

        self.outdoor_air_rng = [
            HumidAir(Tdb=T_outdoor_db, Twb=T_outdoor_wb)
            for T_outdoor_db, T_outdoor_wb in zip(T_outdoor_db_rng, T_outdoor_wb_rng)
        ]

        self.Q_zone_sen_rng = [
            Q_(Q, 'kW')
            for Q in df_Q_zone['Q_sen_load'].values
        ]
        self.Q_zone_lat_rng = [
            Q_(Q, 'kW')
            for Q in df_Q_zone['Q_lat_load'].values
        ]
        self.rel_Q_zone_sen_rng = [
            Q_zone_sen.to('kW') / design_data.Q_zone_sen.to('kW')
            for Q_zone_sen in self.Q_zone_sen_rng
        ]
        self.rel_Q_zone_lat_rng = [
            Q_zone_lat.to('kW') / design_data.Q_zone_lat.to('kW')
            for Q_zone_lat in self.Q_zone_lat_rng
        ]

    def __call__(self, *args, **kwargs) -> Iterator[list[tuple[HumidAir, Quantity, Quantity]]]:
        """Returns an iterator that on each iteration will return a tuple of
        three elements:
        -   the state of outdoor air
        -   the corresponding relative sensible cooling load (ratio of the sensible
            part-load to the sensible design-load)
        -   the corresponding relative latent cooling load (ratio of the latent
            part-load to the latent design-load)
        These tuples can be used as arguments in the method `analyze` of class
        `VAVSingleZoneAirCoolingSystem`.
        """
        return zip(
            self.outdoor_air_rng,
            self.rel_Q_zone_sen_rng,
            self.rel_Q_zone_lat_rng
        )
