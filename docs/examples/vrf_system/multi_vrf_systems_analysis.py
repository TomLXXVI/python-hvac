"""Script to compare the performance of VRF-systems with different outdoor units
with respect to a given building.
"""
from time import time
import pandas as pd
from hvac import Quantity
from hvac.vrf_system import (
    WorkingMode,
    VRFSystem,
    VRFModel,
    Load,
    DesignValues,
    EnergyConsumption,
    find_balance_point
)
import indoor_units as iu
from hvac.energy_estimation import (
    BinTableCreator,
    TimeSegment,
    Month
)
from hvac.charts import LineChart
from hvac.logging import ModuleLogger

Q_ = Quantity

# ----------------------------------------------------------------------------------------------------------------------

# Turn of the DEBUG and WARNING log messages from `hvac.vrf_system.model` and
# `hvac.vrf_system.vrf` by raising the log level to ERROR
logger = ModuleLogger.get_logger('hvac.vrf_system.model')
logger.setLevel(ModuleLogger.ERROR)

logger = ModuleLogger.get_logger('hvac.vrf_system.vrf')
logger.setLevel(ModuleLogger.ERROR)

# ----------------------------------------------------------------------------------------------------------------------
# Set the VRF-system characteristics that are common to all systems

# Equivalent pipe length between outdoor unit and furthest indoor unit
Leq_pipe = Q_(100.0, 'm')

# Maximum height difference between outdoor and indoor units or between indoor units
h_pipe = Q_(10.0, 'm')

# We will analyze the performance of the VRF-systems during the heating season
working_mode = WorkingMode.HEATING

# Retrieve the indoor units of the VRF-system from script `indoor_units.py`
indoor_units = iu.get_indoor_units()


# ----------------------------------------------------------------------------------------------------------------------

def _create_pury_P450():
    """Creates the VFR-system MITSUBISHI PURY P450."""
    pury_P450 = VRFSystem(
        Q_rated=Q_(56, 'kW'),
        W_rated=Q_(12.64, 'kW'),
        model_size=450,
        Leq_pipe=Leq_pipe,
        h_pipe=h_pipe,
        working_mode=working_mode,
        vrf_model_heating=VRFModel.load('./models/pury-p450-heating.pickle'),
        indoor_units=indoor_units
    )
    return pury_P450


def _create_pury_P550():
    """Creates the VFR-system MITSUBISHI PURY P550."""
    pury_P550 = VRFSystem(
        Q_rated=Q_(69, 'kW'),
        W_rated=Q_(16.62, 'kW'),
        model_size=550,
        Leq_pipe=Leq_pipe,
        h_pipe=h_pipe,
        working_mode=working_mode,
        vrf_model_heating=VRFModel.load('./models/pury-p550-heating.pickle'),
        indoor_units=indoor_units
    )
    return pury_P550


def _create_pury_P600():
    """Creates the VFR-system MITSUBISHI PURY P600."""
    pury_P600 = VRFSystem(
        Q_rated=Q_(76.5, 'kW'),
        W_rated=Q_(19.12, 'kW'),
        model_size=600,
        Leq_pipe=Leq_pipe,
        h_pipe=h_pipe,
        working_mode=working_mode,
        vrf_model_heating=VRFModel.load('./models/pury-p600-heating.pickle'),
        indoor_units=indoor_units
    )
    return pury_P600


def _create_load():
    """Creates the building load characteristics based on the result of a
    heating loss calculation of the building and a constant indoor air
    setpoint temperature of the building."""
    load = Load(
        time_segment=None,
        design_values=DesignValues(
            Ql=Q_(46.2, 'kW'),  # heating load of the building at design conditions
            Tia=Q_(20.0, 'degC'),  # design indoor air temperature
            Toa=Q_(-7.0, 'degC')  # design outdoor air temperature
        ),
        Tia_avg=Q_(23.0, 'degC')  # space-load-weighted average indoor air temperature
    )
    return load
# ----------------------------------------------------------------------------------------------------------------------


class EnergyEstimator:
    """Helper class for class `VRFSystemContainer` to calculate the energy
    consumption of a VRF-system during a typical heating season and with
    respect to a given building load."""

    def __init__(self):
        # create the temperature bin table of a typical heating season
        self.bin_table_creator = BinTableCreator(
            file_path="./tmy/tmy_50.811_3.310_2005_2020.csv",
            date_time_column_name='time(UTC)',
            temperature_column_name='T2m',
            bin_limits=(Q_(-10, 'degC'), Q_(22, 'degC')),
            time_segments=[
                TimeSegment('0 - 7 h', (0, 7)),
                TimeSegment('7 - 18 h', (7, 18)),
                TimeSegment('18 - 24 h', (18, 24))
            ]
        )
        self.heating_season = [
            Month.Jan, Month.Feb, Month.Mar, Month.Apr,
            Month.Oct, Month.Nov, Month.Dec
        ]
        self.bin_table = sum(
            self.bin_table_creator.get_bin_table(m)
            for m in self.heating_season
        )

    def get_energy_consumption(
        self,
        vrf_system: VRFSystem,
        load: Load
    ) -> pd.DataFrame:
        """Returns the energy consumption of the given VRF-system during the
        typical heating season with respect to the given building load."""
        energy_consumption = EnergyConsumption(
            bin_table=self.bin_table,
            loads=[load] * 3,
            vrf_system=vrf_system
        ).estimate().round(decimals=3)
        return energy_consumption


class VRFSystemContainer:
    """Class for holding multiple VRF-systems and has methods to retrieve
    the performance parameters of all these VRF-systems for a given building
    load."""

    def __init__(
        self,
        vrf_systems: list[VRFSystem],
        Toa_limits: tuple[Quantity, Quantity],
        load: Load
    ) -> None:
        """Creates the VRF-system container.

        Parameters
        ----------
        vrf_systems:
            list of the VRF-systems in the container
        Toa_limits:
            the lower and upper limit of the outdoor air temperature range
        load:
            the building load
        """
        self.vrf_systems = vrf_systems
        self.Toa_rng = [Q_(Toa, 'degC') for Toa in range(
            int(Toa_limits[0].to('degC').m),
            int(Toa_limits[1].to('degC').m)
        )]
        self.load = load
        self.load.Tia_avg = self.vrf_systems[0].Tia_avg
        self.energy_estimator = EnergyEstimator()

    def get_outdoor_unit_capacity(self) -> tuple[list[Quantity], ...]:
        """Returns a tuple with the outdoor unit capacity in function of
        outdoor air temperature of each VRF-system in the container, preserving
        the order in which the VRF-systems were added to the container."""
        r = []
        for vrf_sys in self.vrf_systems:
            Qh_ou_rng = [
                vrf_sys.get_full_load_outdoor_unit_capacity(None, Toa)
                for Toa in self.Toa_rng
            ]
            r.append(Qh_ou_rng)
        return tuple(r)

    def get_system_capacity(self) -> tuple[list[Quantity], ...]:
        """Returns a tuple with the system capacity in function of
        outdoor air temperature of each VRF-system in the container, preserving
        the order in which the VRF-systems were added to the container."""
        r = []
        for vrf_sys in self.vrf_systems:
            Qh_ou_rng = [
                vrf_sys.get_available_capacity(Toa)
                for Toa in self.Toa_rng
            ]
            r.append(Qh_ou_rng)
        return tuple(r)

    def get_balance_point(self) -> tuple[tuple[Quantity, Quantity], ...]:
        r = []
        for vrf_sys in self.vrf_systems:
            T_bal, Q_bal = find_balance_point(
                self.load,
                vrf_sys,
                (self.Toa_rng[0], self.Toa_rng[-1])
            )
            r.append((T_bal, Q_bal))
        return tuple(r)

    def get_input_power(self) -> tuple[list[Quantity], ...]:
        """Returns a tuple with the input power in function of
        outdoor air temperature of each VRF-system in the container, preserving
        the order in which the VRF-systems were added to the container."""
        r = []
        Ql_rng = [self.load(Toa) for Toa in self.Toa_rng]
        for vrf_sys in self.vrf_systems:
            Qh_ou_rng = [
                vrf_sys.get_available_capacity(Toa)
                for Toa in self.Toa_rng
            ]
            PLR_rng = [
                min(Q_(1.0, 'frac'), Ql / Qh_ou)
                for Ql, Qh_ou in zip(Ql_rng, Qh_ou_rng)
            ]
            W_in_rng = [
                vrf_sys.get_input_power(None, Toa, PLR)
                for Toa, PLR in zip(self.Toa_rng, PLR_rng)
            ]
            r.append(W_in_rng)
        return tuple(r)
    
    def get_cop(self) -> tuple[list[Quantity], ...]:
        """Returns a tuple with the COP in function of outdoor air temperature
        of each VRF-system in the container, preserving the order in which the
        VRF-systems were added to the container."""
        r = []
        Ql_rng = [self.load(Toa) for Toa in self.Toa_rng]
        tup_Qh_vrf_rng = self.get_system_capacity()
        tup_W_in_rng = self.get_input_power()
        for Qh_vrf_rng, W_in_rng in zip(tup_Qh_vrf_rng, tup_W_in_rng):
            COP_rng = [
                min(Qh_vrf, Ql) / W_in
                for Qh_vrf, Ql, W_in in zip(Qh_vrf_rng, Ql_rng, W_in_rng)
            ]
            r.append(COP_rng)
        return tuple(r)

    def get_energy_consumption(self) -> tuple[pd.DataFrame]:
        """Returns a tuple with the energy consumption table of each VRF-system
        in the container with reference to a typical heating season, preserving
        the order in which the VRF-systems were added to the container."""
        r = []
        for vrf_sys in self.vrf_systems:
            df = self.energy_estimator.get_energy_consumption(
                vrf_sys, self.load
            )
            r.append(df)
        return tuple(r)


# ----------------------------------------------------------------------------------------------------------------------

def main():
    start_ = time()

    # create the VRF-systems we want to compare
    pury_P450 = _create_pury_P450()
    pury_P550 = _create_pury_P550()
    pury_P600 = _create_pury_P600()

    # create a container for our VRF-systems and the building load
    vrf_system_container = VRFSystemContainer(
        vrf_systems=[pury_P450, pury_P550, pury_P600],
        Toa_limits=(Q_(-20, 'degC'), Q_(16.0, 'degC')),
        load=_create_load()
    )
    # get the system capacity of our VRF-systems
    # in function of outdoor air temperature
    tup_Qh_vrf_rng = vrf_system_container.get_system_capacity()

    # get the input power of our VRF-systems
    # in function of outdoor air temperature
    tup_W_in_rng = vrf_system_container.get_input_power()

    # get the COP of our VRF-systems
    # in function of outdoor air temperature
    tup_COP_rng = vrf_system_container.get_cop()

    labels = ('PURY-P450', 'PURY-P550', 'PURY-P600')

    # get the balance point of our VRF-systems
    tup_bal_pt = vrf_system_container.get_balance_point()
    for bal_pt, lbl in zip(tup_bal_pt, labels):
        print(
            f"{lbl}: balance temperature = {bal_pt[0]:~P.1f}; "
            f"balance load = {bal_pt[1]:~P.3f}"
        )
    print()

    # draw a line chart of the system capacity of our VRF-systems
    # and the building load line in function of outdoor air temperature
    chart_1 = LineChart()
    for Qh_vrf_rng, lbl in zip(tup_Qh_vrf_rng, labels):
        chart_1.add_xy_data(
            label=lbl,
            x1_values=[Toa.to('degC').m for Toa in vrf_system_container.Toa_rng],
            y1_values=[Q.to('kW').m for Q in Qh_vrf_rng]
        )
    chart_1.add_xy_data(
        label='load line',
        x1_values=[Toa.to('degC').m for Toa in vrf_system_container.Toa_rng],
        y1_values=[
            vrf_system_container.load(Toa).to('kW').m
            for Toa in vrf_system_container.Toa_rng
        ],
        style_props={'linestyle': '--'}
    )
    chart_1.x1.add_title('outdoor temperature, °C')
    chart_1.y1.add_title('heat power, kW')
    chart_1.add_legend()
    chart_1.show()

    # draw a line chart of the input power of our VRF-systems
    # in function of outdoor air temperature
    chart_2 = LineChart()
    for W_in_rng, lbl in zip(tup_W_in_rng, labels):
        chart_2.add_xy_data(
            label=lbl,
            x1_values=[Toa.to('degC').m for Toa in vrf_system_container.Toa_rng],
            y1_values=[W.to('kW').m for W in W_in_rng]
        )
    chart_2.x1.add_title('outdoor temperature, °C')
    chart_2.y1.add_title('input power, kW')
    chart_2.add_legend()
    chart_2.show()

    # draw a line chart of the COP of our VRF-systems
    # in function of outdoor air temperature
    chart_3 = LineChart()
    for COP_rng, lbl in zip(tup_COP_rng, labels):
        chart_3.add_xy_data(
            label=lbl,
            x1_values=[Toa.to('degC').m for Toa in vrf_system_container.Toa_rng],
            y1_values=[COP.m for COP in COP_rng]
        )
    chart_3.x1.add_title('outdoor temperature, °C')
    chart_3.y1.add_title('COP')
    chart_3.add_legend()
    chart_3.show()

    # print the energy consumption tables of our VRF-systems during a
    # typical heating season
    tup_energy_consumption = vrf_system_container.get_energy_consumption()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        for df, lbl in zip(tup_energy_consumption, labels):
            print(lbl)
            print(df, end='\n\n')

    end_ = time()
    duration = end_ - start_
    print(duration)


if __name__ == '__main__':
    main()
