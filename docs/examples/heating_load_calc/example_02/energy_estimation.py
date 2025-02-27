"""EXAMPLE 02 (part 2)
----------------------
ESTIMATION OF HEATING ENERGY CONSUMPTION USING THE BIN TABLE METHOD

In part 2 of example 02, we will estimate the heating energy consumption of the 
building using the bin table method with TMY data valid for the geographic 
location where the building is situated. The heating load characteristics are 
derived from the design heating load calculation in part 1.
"""
import pandas as pd
from hvac import Quantity
from hvac.energy_estimation import BinTableCreator, HeatingLoad, TimeSegment

Q_ = Quantity

# Create a bin table of the month January using the TMY data from the specified
# file. Divide each day of the month (including weekends) in three time
# segments.
btc = BinTableCreator(
    file_path='tmy_50.911_3.192_2005_2020.csv',
    date_time_column_name='time(UTC)',
    temperature_column_name='T2m',
    time_segments=[
        TimeSegment('0h-7h', (0, 7)),  # 07:00 not included
        TimeSegment('7h-18h', (7, 18)),
        TimeSegment('18h-24h', (7, 24))
    ]
)
bt_jan = btc.get_bin_table(1)

# Display the bin table:
with pd.option_context(
    'display.max_rows', None,
    'display.max_columns', None,
    'display.width', 800
):
    print(bt_jan)


# Define the heating load characteristics at nighttime and at daytime:
load_nighttime = HeatingLoad(
    T_int=Q_(16, 'degC'),
    H_trm=Q_(0.271, 'kW / K'),
    V_dot_ven=Q_(138.239, 'm**3 / hr'),
    eta_sys=Q_(90.0, 'pct'),
    T_ext_min=Q_(-5.7, 'degC'),
    time_segment=None
)


load_daytime = HeatingLoad(
    T_int=Q_(22, 'degC'),
    H_trm=Q_(0.271, 'kW / K'),
    V_dot_ven=Q_(138.239, 'm**3 / hr'),
    eta_sys=Q_(90.0, 'pct'),
    T_ext_min=Q_(-5.7, 'degC'),
    time_segment=None
)


# Get the energy consumption of the heating system in January by iterating over
# the temperature bins and time segments:
Q_in_tbl = []
for T_bin in bt_jan.index:
    Q_in_row = []
    for t_seg_lbl in bt_jan.columns:
        if t_seg_lbl == '7h-18h':
            # From 7h until 18h: daytime load.
            load = load_daytime
        else:
            # Otherwise: nighttime load.
            load = load_nighttime
        # Get the number of hours in the current temperature bin and time segment:
        load.num_hours = bt_jan.loc[T_bin, t_seg_lbl]
        # Get the average temperature of the temperature bin:
        load.T_ext = Q_(T_bin.mid, 'degC')
        # Get the energy consumption in the current temperature bin and time segment:
        Q_in_row.append(load.Q_in.to('kWh').m)
    Q_in_tbl.append(Q_in_row)

# Transform the energy consumption list of lists into a dataframe with the same
# temperature bin index and with the same time segment columns as the bin table
# dataframe:
Q_in_tbl = pd.DataFrame(
    data=Q_in_tbl,
    index=bt_jan.index,
    columns=bt_jan.columns
)

# Add a row with the sums per column (= monthly energy consumption in a
# temperature bin):
Q_in_tbl.loc['sum'] = Q_in_tbl.sum()

# Add a column with the sums per row (= monthly energy consumption in a
# time segment):
Q_in_tbl.loc[:, 'sum'] = Q_in_tbl.sum(axis=1)

with pd.option_context(
    'display.max_rows', None,
    'display.max_columns', None,
    'display.width', 800
):
    print(Q_in_tbl)
