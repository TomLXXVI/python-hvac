"""
EXAMPLE 1: PANEL RADIATOR
-------------------------

This script demonstrates the usage of class `PanelRadiator` in subpackage
`hvac.radiant_heat_emitter`. This class can be used to size panel radiators and
to analyze the performance of a panel radiator.
"""
import warnings
from hvac import Quantity
from hvac.radiant_emitter import PanelRadiator
from hvac.charts import LineChart

warnings.filterwarnings('ignore')

Q_ = Quantity


# Create a radiator from its nominal specs taken from the manufacturer's
# datasheet:
radiator = PanelRadiator(
    Qe_dot_nom=Q_(2000, 'W'),
    Tw_sup_nom=Q_(80, 'degC'),
    Tw_ret_nom=Q_(60, 'degC'),
    Ti_nom=Q_(20, 'degC'),
    n_exp=1.33
)
print(
    f"Nominal volume flow rate of the radiator = "
    f"{radiator.Vw_dot_nom.to('L / min'):~P.3f}"
)

# Find the heat output of the radiator when:
# - the room temperature is 22 °C
# - the volume flow rate of water through the radiator is 1.5 L/min
# - the supply water temperature is 70 °C
Qe_dot = radiator(
    Ti=Q_(22, 'degC'),
    Vw_dot=Q_(1.5, 'L/min'),
    Tw_sup=Q_(70, 'degC')
)
print(
    "Heat output @ Ti = 22 °C, Vw_dot = 1.5 L/min, Tw_sup = 70 °C: "
    f"{Qe_dot.to('kW'):~P.3f}"
)

# Find the required volume flow rate of water through the radiator when:
# - the room temperature is 20 °C
# - the supply water temperature is 80 °C
# - the needed heat output is 1.5 kW
Vw_dot = radiator(
    Ti=Q_(20, 'degC'),
    Tw_sup=Q_(80, 'degC'),
    Qe_dot=Q_(1.5, 'kW')
)
print(
    "Required water volume flow rate for heat output 1.5 kW "
    "@ Ti = 20 °C and Tw_sup = 80 °C: "
    f"{Vw_dot.to('L / min'):~P.3f}"
)

# Find the required supply water temperature through the radiator when:
# - the room temperature is 20 °C
# - the water volume flow rate is 0.6 L/min
# - the needed heat output is 1.5 kW
Tw_sup = radiator(
    Ti=Q_(20, 'degC'),
    Vw_dot=Q_(0.6, 'L / min'),
    Qe_dot=Q_(1.5, 'kW')
)
print(
    "Required supply water temperature for heat output 1.5 kW "
    "@ Ti = 20 °C and Vw_dot = 0.6 L/min: "
    f"{Tw_sup.to('degC'):~P.3f}"
)

# Plot the QV-characteristic of the radiator for a given supply water
# temperature and room temperature. From the QV-characteristic of the
# radiator one can read what volume flow rate of water is needed to get a
# certain heat output of the radiator when the supply water temperature is
# fixed.
Vw_dot_rng, Qe_dot_rng = radiator.QV_characteristic(
    Tw_sup=Q_(80, 'degC'),
    Ti=Q_(20, 'degC'),
    Vw_dot_max=radiator.Vw_dot_nom
)
qv_chart = LineChart()
qv_chart.add_xy_data(
    label='QV-characteristic',
    x1_values=[Vw.to('L /min').m for Vw in Vw_dot_rng],
    y1_values=[Qe.to('W').m for Qe in Qe_dot_rng]
)
qv_chart.x1.add_title('volume flow rate, L/min')
qv_chart.y1.add_title('heat output, W')
qv_chart.show()

# Plot the QT-characteristic of the radiator for a given volume flow rate
# and room temperature. From the QT-characteristic of the radiator one can
# read what supply water temperature is needed to get a certain heat output
# of the radiator when the water volume flow rate is fixed.
Tw_sup_rng, Qe_dot_rng = radiator.QT_characteristic(
    Vw_dot=Q_(1.5, 'L / min'),
    Ti=Q_(20, 'degC'),
    Tw_sup_max=radiator.Tw_sup_nom
)
qt_chart = LineChart()
qt_chart.add_xy_data(
    label='QT-characteristic',
    x1_values=[Tw_sup.to('degC').m for Tw_sup in Tw_sup_rng],
    y1_values=[Qe.to('W').m for Qe in Qe_dot_rng]
)
qt_chart.x1.add_title('supply temperature, °C')
qt_chart.y1.add_title('heat output, W')
qt_chart.show()
