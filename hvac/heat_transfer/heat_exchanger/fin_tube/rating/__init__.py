from .circular_fin_tube_heat_exchanger import (
    CircularFinTubeCrossFlowHeatExchanger,
    CircularFinTubeCounterFlowHeatExchanger
)
from .plain_fin_tube_heat_exchanger import (
    PlainFinTubeCrossFlowHeatExchanger,
    PlainFinTubeCounterFlowHeatExchanger
)

CFT_CR_HEX = CircularFinTubeCrossFlowHeatExchanger
CFT_CO_HEX = CircularFinTubeCounterFlowHeatExchanger

PFT_CR_HEX = PlainFinTubeCrossFlowHeatExchanger
PFT_CO_HEX = PlainFinTubeCounterFlowHeatExchanger
