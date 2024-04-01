"""
ENERGY ESTIMATION FOR HEATING (part 2):
Estimation of the energy consumption of an (air-to-water) heat pump using the
bin table method.
"""
import pandas as pd
from hvac import Quantity
from hvac.energy_estimation import (
    TimeSegment,
    BinTableCreator,
    EnergyEstimator,
    HeatPump,
    HeatingLoad
)

Q_ = Quantity

# Get the temperature bin table for December based on TMY data of the geographic
# location being considered. Each day is divided into two time segments.
# (Note: the bin limits should not exceed the lower and upper limit of the
# outdoor temperature range used in the heat pump performance table. Otherwise,
# a `ValueError` exception may be raised due to exceeding interpolation limits.)
bin_table_creator = BinTableCreator(
    file_path='tmy_ghent.csv',  # file with TMY data retrieved from PVGIS
    date_time_column_name='time',
    temperature_column_name='T2m',
    bin_limits=(Q_(-10, 'degC'), Q_(35, 'degC')),
    time_segments=[
        TimeSegment('0-12 h', (0, 12)),
        TimeSegment('12-24 h', (12, 24))
    ]
)
bin_december = bin_table_creator.get_bin_table(12)


# Define the heating load in time segment 1
load_seg1 = HeatingLoad(
    T_int=Q_(20, 'degC'),
    T_ext_min=Q_(-8, 'degC'),
    H_trm=Q_(6 / 28, 'kW / K'),
    V_dot_ven=Q_(100, 'm**3 / hr'),
    Q_dot_ihg=Q_(500, 'W'),
    eta_sys=Q_(80, 'pct'),
    time_segment=None
)

# Define the heating load in time segment 2
load_seg2 = HeatingLoad(
    T_int=Q_(20, 'degC'),
    T_ext_min=Q_(-8, 'degC'),
    H_trm=Q_(6 / 28, 'kW / K'),
    V_dot_ven=Q_(200, 'm**3 / hr'),
    Q_dot_ihg=Q_(1000, 'W'),
    eta_sys=Q_(80, 'pct'),
    time_segment=None
)

# Define the heat pump:
heat_pump = HeatPump(file_path="hp_atlantic_alfea_extensa_AI5_W35.csv")


# Estimate the energy consumption of the heat pump:
energy_estimator = EnergyEstimator(
    bin_table=bin_december,
    loads=[load_seg1, load_seg2],
    heat_pump=heat_pump
)
energy_use = energy_estimator.estimate()
with pd.option_context(
    'display.max_rows', None,
    'display.max_columns', None,
    'display.width', 1000
): print(energy_use)
