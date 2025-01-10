from typing import Optional
import warnings
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import SimpleCompressionCycle, StateContainer
import matplotlib.pyplot as plt
from .. import Quantity
from ..fluids import Fluid


Q_ = Quantity
warnings.filterwarnings('ignore', category=UserWarning)


class StandardVaporCompressionCycle:

    def __init__(
        self,
        Refrigerant: Fluid,
        evaporationTemperature: Quantity,
        condensationTemperature: Quantity,
        evaporatorSuperheat: Optional[Quantity] = None,
        subCooling: Optional[Quantity] = None,
        suctionLineSuperheat: Optional[Quantity] = None,
        isentropicEfficiency: Optional[Quantity] = None
    ):
        """
        Create a standard vapor compression cycle

        Parameters
        ----------
        Refrigerant: Fluid
            The refrigerant used in the vapor compression cycle.
        evaporationTemperature: Quantity
            Evaporation temperature of the cycle.
        condensationTemperature: Quantity
            Condensation temperature of the cycle.
        evaporatorSuperheat: Quantity, optional
            Degree of superheat at the evaporator's exit.
        subCooling: Quantity, optional
            Degree of subcooling at the condenser's exit.
        suctionLineSuperheat: Quantity, optional
            Amount of superheat taken up by the refrigerant vapor in the
            suction line.
        isentropicEfficiency: Quantity, optional
            Isentropic efficiency of compressor. If None, isentropic compression
            is assumed.
        """
        self.Refrigerant = Refrigerant
        self.evaporationTemperature = evaporationTemperature
        self.condensationTemperature = condensationTemperature
        self.evaporatorSuperheat = evaporatorSuperheat or Q_(0.001, 'delta_degC')
        self.subCooling = subCooling or Q_(0.001, 'delta_degC')
        self.suctionLineSuperheat = suctionLineSuperheat or Q_(0.001, 'delta_degC')
        self.isentropicEfficiency = isentropicEfficiency or Q_(1.0, 'frac')

        self.saturatedLiquid = Refrigerant(
            T=self.condensationTemperature,
            x=Q_(0, 'pct')
        )
        self.saturatedVapor = Refrigerant(
            T=self.evaporationTemperature,
            x=Q_(100, 'pct')
        )
        self.subCooledLiquid = Refrigerant(
            T=self.saturatedLiquid.T - self.subCooling,
            P=self.saturatedLiquid.P
        )
        try:
            self.liquidVaporMixture = Refrigerant(
                P=self.saturatedVapor.P,
                h=self.subCooledLiquid.h
            )
        except IndexError:
            self.liquidVaporMixture = Refrigerant(
                P=self.saturatedVapor.P,
                h=self.subCooledLiquid.h,
                x=Q_(0, 'pct')
            )
        self.superheatedVapor = Refrigerant(
            T=self.saturatedVapor.T + self.evaporatorSuperheat,
            P=self.saturatedVapor.P
        )
        self.suctionGas = Refrigerant(
            T=self.superheatedVapor.T + self.suctionLineSuperheat,
            P=self.saturatedVapor.P
        )
        try:
            self.dischargeGas_isentropic = Refrigerant(
                P=self.saturatedLiquid.P,
                s=self.suctionGas.s
            )
        except IndexError:
            self.dischargeGas_isentropic = Refrigerant(
                P=self.saturatedLiquid.P,
                s=self.suctionGas.s,
                T=self.condensationTemperature
            )
        self.isentropicCompressionWork = self.dischargeGas_isentropic.h - self.suctionGas.h
        self.actualCompressionWork = self.isentropicCompressionWork / self.isentropicEfficiency
        try:
            self.dischargeGas = Refrigerant(
                P=self.saturatedLiquid.P,
                h=self.suctionGas.h + self.actualCompressionWork
            )
        except IndexError:
            self.dischargeGas = Refrigerant(
                P=self.saturatedLiquid.P,
                h=self.suctionGas.h + self.actualCompressionWork,
                T=self.condensationTemperature
            )

    @property
    def netRefrigerationEffect(self) -> Quantity:
        h2 = self.superheatedVapor.h
        h1 = self.liquidVaporMixture.h
        return h2 - h1


class LogPhDiagram:

    def __init__(self, refrigerant: Fluid, **kwargs):
        """
        Create log(P)-h diagram for `refrigerant`.

        Keyword arguments
        -----------------
        size: Tuple[float, float], optional
            Size of the diagram.
        dpi: int, optional
            Dots per inch of diagram.
        """
        self.refrigerant = refrigerant

        size = kwargs.get('size')
        dpi = kwargs.get('dpi')

        fig, axis = plt.subplots(figsize=size, dpi=dpi, layout='constrained')

        self._propertyPlot = PropertyPlot(
            self.refrigerant.coolprop_abstract_state,
            'PH',
            unit_system='EUR'
        )

        # Close the figure made by `PropertyPlot` and pass it our own `Figure`
        # and `Axis` object.
        plt.close(self._propertyPlot.figure)
        self._propertyPlot.figure = fig
        self._propertyPlot.axis = axis
        self._propertyPlot.get_axis_limits()
        # noinspection PyProtectedMember
        self._propertyPlot._plot_default_annotations()
        self._propertyPlot.calc_isolines()

        self._scc_cycle: Optional[SimpleCompressionCycle] = None
        self._stateContainer: Optional[StateContainer] = None
        self.cycle: Optional[StandardVaporCompressionCycle] = None

    def setCycle(self, cycle: StandardVaporCompressionCycle):
        """Pass a `StandardVaporCompressionCycle` instance to draw this cycle on
        the log(P)-h diagram.
        """
        self.cycle = cycle
        self._scc_cycle = SimpleCompressionCycle(
            self.refrigerant.coolprop_abstract_state,
            'PH',
            unit_system='EUR'
        )
        self._scc_cycle.simple_solve(
            T0=cycle.suctionGas.T.to('K').magnitude,
            p0=cycle.suctionGas.P.to('Pa').magnitude,
            T2=cycle.subCooledLiquid.T.to('K').magnitude,
            p2=cycle.subCooledLiquid.P.to('Pa').magnitude,
            eta_com=cycle.isentropicEfficiency.to('frac').magnitude
        )
        self._scc_cycle.steps = 50
        self._stateContainer = self._scc_cycle.get_state_changes()
        plt.close(self._scc_cycle.figure)

    def show(self):
        """Show the diagram on screen."""
        if self._stateContainer is not None:
            self._propertyPlot.draw_process(self._stateContainer)
            self._annotate_points()
        self._propertyPlot.show()

    def _annotate_points(self):
        x_data = [
            self.cycle.liquidVaporMixture.h.to('kJ / kg').m,
            self.cycle.suctionGas.h.to('kJ / kg').m,
            self.cycle.dischargeGas.h.to('kJ / kg').m,
            self.cycle.subCooledLiquid.h.to('kJ / kg').m
        ]
        y_data = [
            self.cycle.liquidVaporMixture.P.to('bar').m,
            self.cycle.suctionGas.P.to('bar').m,
            self.cycle.dischargeGas.P.to('bar').m,
            self.cycle.subCooledLiquid.P.to('bar').m
        ]
        labels = [1, 2, 3, 4]
        xytext_list = [(0, -15), (0, -15), (0, 7), (0, 7)]
        for x, y, label, xytext in zip(x_data, y_data, labels, xytext_list):
            label = f"{label}"
            self._propertyPlot.axis.annotate(
                label,
                (x, y),
                textcoords='offset points',
                xytext=xytext,
                ha='center'
            )

    @property
    def figure(self):
        return self._propertyPlot.figure
