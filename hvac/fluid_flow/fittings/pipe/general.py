from typing import Optional
import math
from hvac import Quantity
from hvac.fluids import Fluid, FluidState
from hvac.fluid_flow.schedule import pipe_schedule_40
from hvac.fluid_flow.conduit import Pipe
from ..abstract_fitting import AbstractFitting


Q_ = Quantity
Water = Fluid('Water')


class FlowCoefficient:
    """The flow coefficient of a valve relates the volume flow rate through
    the valve with the pressure loss across this valve.

    In SI base units:
    ```
        V_dot = Av * math.sqrt(delta_P / rho)
    ```
    where:
        V_dot = volume flow rate in 'm**3 / s'
        delta_P = pressure drop across the valve in 'Pa'
        rho = (average) mass density of the fluid in 'kg / m**3'
        Av = the flow coefficient based on SI base units

    Another flow coefficient that is most frequently used, is defined like:
    ```
        V_dot = Kv * math.sqrt(delta_P / (rho / rho_w_15))
    ```
    where:
        V_dot = volume flow rate in 'm**3 / hr'
        delta_P = pressure drop across the valve in 'bar'
        rho = (average) mass density of the fluid in 'kg / m**3'
        rho_w_15 = mass density of water at 15 °C in 'kg / m**3'
        Kv = the flow coefficient based on the volume flow rate in m**3/hr and
             the pressure drop in bar.

    The ratio `rho / rho_w_15` is also called the specific gravity of the fluid.
    """
    water_15 = Water(T=Q_(15.0, 'degC'), P=Q_(101.325, 'kPa'))
    rho_15 = water_15.rho.to('kg / m ** 3').magnitude

    @classmethod
    def to_Kv(cls, Av: float) -> float:
        """Converts the flow coefficient Av, which is based on SI base units
        (volume flow rate in m³/s, pressure drop in Pa, and mass density in
        kg/m³), to the flow coefficient Kv, which is based on the volume flow
        rate being expressed in units m³/h, the pressure drop in units
        bar, and the mass density of water at 15 °C.
        """
        Kv = Av * 3.6e5 * math.sqrt(10) / math.sqrt(cls.rho_15)
        return Kv

    @classmethod
    def to_Av(cls, Kv: float) -> float:
        """Converts the flow coefficient Kv (float), which is based on volume
        flow rate in m³/hr and pressure drop in bar, to the flow coefficient Av,
        which is based on SI base units (volume flow rate in m³/s, pressure drop
        in Pa, and mass density in kg/m³).
        """
        Av = Kv * math.sqrt(cls.rho_15) / (3.6e5 * math.sqrt(10))
        return Av

    @classmethod
    def get_Kv(
        cls,
        volume_flow_rate: Quantity,
        pressure_drop: Quantity,
        fluid: Optional[FluidState] = None
    ) -> float:
        """Returns the flow coefficient Kv of the valve, which corresponds
        with the given volume flow rate through the valve and the given
        pressure drop across the valve.

        Parameters
        ----------
        volume_flow_rate:
            The volume flow rate through the valve.
        pressure_drop:
            The pressure drop across the valve.
        fluid:
            The state of the fluid flowing through the valve.

        Returns
        -------
        The corresponding flow coefficient Kv of the valve as a float.
        """
        V = volume_flow_rate.to('m ** 3 / hr').magnitude
        dp = pressure_drop.to('bar').magnitude
        specific_gravity = 1.0
        if isinstance(fluid, FluidState):
            rho_fluid = fluid.rho.to('kg / m ** 3').magnitude
            specific_gravity = rho_fluid / cls.rho_15
        Kv = V / math.sqrt(dp / specific_gravity)
        return Kv


class ResistanceCoefficient:
    """The resistance coefficient of a fitting, aka 'zeta-value', relates the
    pressure drop across the fitting with the velocity pressure at the inlet
    or outlet of the fitting.
    """
    @staticmethod
    def from_Av(Av: float, diameter: Quantity) -> float:
        """Converts the flow coefficient Av (float) of a valve to its equivalent
        resistance coefficient (float).

        Parameters
        ----------
        Av:
            The flow coefficient of the valve based on SI base units (volume
            flow rate in m³/s, pressure drop in Pa, and mass density in kg/m³).
        diameter:
            The inner diameter of the pipe to which the velocity pressure
            applies.
        """
        D = diameter.to('m').magnitude
        zeta = math.pi ** 2 * D ** 4 / (8.0 * Av ** 2)
        return zeta

    @staticmethod
    def from_Kv(Kv: float, diameter: Quantity) -> float:
        """Converts the flow coefficient Kv (float) of a valve to its equivalent
        resistance coefficient (float).

        Parameters
        ----------
        Kv:
            The flow coefficient of the valve based on volume flow rate in
            'm**3/hr' and pressure drop in 'bar'.
        diameter:
            The inner diameter of the pipe to which the velocity pressure
            applies.
        """
        Av = FlowCoefficient.to_Av(Kv)
        return ResistanceCoefficient.from_Av(Av, diameter)

    @staticmethod
    def from_ELR(ELR: float, diameter: Quantity) -> float:
        """Converts the 'Equivalent Length Ratio' (ELR) of a fitting to its
        equivalent resistance coefficient (see Crane's TP 410M).

        Parameters
        ----------
        ELR:
            The equivalent length ratio of the fitting (see Crane's Technical
            Paper No. 410M, metric version).
        diameter:
            The inner diameter of the pipe to which the velocity pressure
            applies.
        """
        D = diameter.to('mm').magnitude
        D_nom = pipe_schedule_40.get_closest_nominal_diameter(diameter)
        D_sch40 = pipe_schedule_40.get_internal_diameter(D_nom).to('mm').magnitude
        e = 0.046  # wall roughness of steel pipe schedule 40 in mm.
        f = 0.25 / (math.log10(e / (3.7 * D_sch40))) ** 4
        return f * ELR * (D / D_sch40) ** 4


class PipeFitting(AbstractFitting):
    """Represents a fitting in a pipe."""

    def __init__(
        self,
        pipe: Pipe,
        ID: str = '',
        Kv: float = float('nan'),
        zeta: float = float('nan'),
        zeta_inf: float = float('nan'),
        zeta_d: float = float('nan'),
        ELR: float = float('nan')
    ) -> None:
        """Creates a `PipeFitting` object.

        Parameters
        ----------
        pipe:
            The pipe to which the fittings needs to be added.
        ID:
            Identifier for the fitting in the pipe.
        Kv: optional
            The flow coefficient Kv of the valve.
        zeta: optional
            The resistance coefficient of the fitting.
        zeta_inf: optional
            A second resistance coefficient of the same fitting (see 3K-method).
        zeta_d: optional
            A third resistance coefficient of the same fitting (see 3K-method).
        ELR: optional
            The 'Equivalent Length Ratio' of the fitting.

        Notes
        -----
        1. Zeta-values are of type float. They are used to relate the pressure
        drop across the fitting to the volume flow rate through the fitting in
        the equation `dP = zeta * rho * v**2 / 2` where `dP` is the pressure drop
        expressed in 'Pa', `rho` is the mass density of the fluid expressed in
        'kg / m**3', and `v` is the flow velocity expressed in 'm / s'.

        2. If `zeta_inf` or `zeta_d` is specified, then `zeta`, `zeta_inf`, and
        `zeta_d` must all be specified, to calculate the pressure drop across
        the fitting using the 3K-method.

        3. Kv-values are of type float. They are used to relate the volume flow
        rate through a valve to the pressure drop across a valve in the equation
        `V_dot = Kv * math.sqrt(dP / (rho/rho_w_15))` where `V_dot` is the volume
        flow rate through the valve expressed in 'm**3 / hr', `dP` is the
        pressure drop across the valve expressed in `bar`, and `rho/rho_w_15` is
        the specific gravity of the fluid, being the ratio of the mass density
        `rho` of the fluid to the mass density `rho_w_15` of water at standard
        pressure (101_325 Pa) and a temperature of 15 °C.

        4. ELR-values are of type float. ELR relates the resistance coefficient
        `zeta` to the Darcy friction factor for fully turbulent flow:
        `zeta = ELR * fT` where fT is the Darcy friction factor for fully
        turbulent flow (see Crane's TP 410M).
        """
        super().__init__(ID)
        self.pipe = pipe
        self.Kv = Kv
        self._zeta = zeta
        self.zeta_inf = zeta_inf
        self.zeta_d = zeta_d
        self.ELR = ELR

    def _deltaP_from_Kv(self) -> Quantity:
        zeta = ResistanceCoefficient.from_Kv(
            self.Kv,
            self.pipe.cross_section.hydraulic_diameter
        )
        rho = self.pipe.fluid.rho
        v = self.pipe.velocity
        dp = zeta * rho * v ** 2 / 2
        return dp.to('Pa')

    def _deltaP_from_1K(self) -> Quantity:
        rho = self.pipe.fluid.rho
        v = self.pipe.velocity
        dp = self._zeta * rho * v ** 2 / 2
        return dp.to('Pa')

    def _deltaP_from_3K(self) -> Quantity:
        rho = self.pipe.fluid.rho
        mu = self.pipe.fluid.absolute_viscosity
        v = self.pipe.velocity
        pv = rho * v ** 2 / 2
        D = self.pipe.cross_section.hydraulic_diameter
        Re = v * D * rho / mu
        dp = (
            self._zeta / Re
            + self.zeta_inf * (1 + self.zeta_d / D ** 0.3)
        ) * pv
        return dp.to('Pa')

    def _deltaP_from_ELR(self) -> Quantity:
        rho = self.pipe.fluid.rho
        v = self.pipe.velocity
        zeta = ResistanceCoefficient.from_ELR(
            self.ELR,
            self.pipe.cross_section.hydraulic_diameter
        )
        dp = zeta * rho * v ** 2 / 2
        return dp.to('Pa')

    @property
    def pressure_drop(self) -> Quantity:
        """Returns the pressure drop across the fitting."""
        if not math.isnan(self.Kv):
            return self._deltaP_from_Kv()
        if not math.isnan(self.zeta_inf):
            return self._deltaP_from_3K()
        if not math.isnan(self.ELR):
            return self._deltaP_from_ELR()
        if not math.isnan(self._zeta):
            return self._deltaP_from_1K()

    @property
    def zeta(self) -> float:
        """Returns the resistance coefficient of the fitting."""
        if not math.isnan(self.Kv):
            return ResistanceCoefficient.from_Kv(
                self.Kv,
                self.pipe.cross_section.hydraulic_diameter
            )
        if not math.isnan(self.zeta_inf):
            dp = self._deltaP_from_3K()
            rho = self.pipe.fluid.rho
            v = self.pipe.velocity
            pv = rho * v ** 2 / 2
            return (dp / pv).magnitude
        if not math.isnan(self.ELR):
            return ResistanceCoefficient.from_ELR(
                self.ELR,
                self.pipe.cross_section.hydraulic_diameter
            )
        if not math.isnan(self._zeta):
            return self._zeta


class BalancingValve(AbstractFitting):
    """Represents a valve for static balancing."""

    def __init__(
        self,
        pipe: Pipe,
        ID: str = '',
        pressure_drop_full_open: Quantity = Q_(float('nan'), 'Pa')
    ) -> None:
        """Creates a `BalancingValve` object.

        Parameters
        ----------
        pipe:
            The pipe in which the balancing valve is to be installed.
        ID:
            Identifier for the balancing valve in the pipe network.
        pressure_drop_full_open:
            The target pressure drop across the valve when it is fully open.
            Normally it is recommended that the minimum pressure drop across
            the valve is at least 3 kPa.
        """
        super().__init__(ID)
        self.pipe = pipe
        self.fluid = self.pipe.fluid
        self.diameter = self.pipe.cross_section.hydraulic_diameter
        self.volume_flow_rate = self.pipe.volume_flow_rate
        self.pressure_drop_full_open = pressure_drop_full_open
        self._pressure_drop: Quantity = Q_(float('nan'), 'Pa')
        self.Kvs: float = float('nan')
        self.Kvr: float = float('nan')

    def calculate_preliminary_Kvs(self) -> float:
        """Calculates the required flow coefficient Kvs of the fully open
        balancing valve so that the target pressure drop is present across the
        fully open valve when the design volume flow rate in the pipe flows
        through this fully open valve.
        """
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        dp_100 = self.pressure_drop_full_open.to('Pa').magnitude
        Avs = V / math.sqrt(dp_100 / rho)
        Kvs = FlowCoefficient.to_Kv(Avs)
        return Kvs

    def set_Kvs(self, Kvs: float):
        """Based on the preliminary, minimal Kvs-value, a commercially available
        balancing valve with a given Kvs-value can be selected. With this method
        the actual Kvs-value of the selected balancing valve can be set.
        """
        self.Kvs = Kvs
        self._calculate_pressure_drop_full_open()

    def _calculate_pressure_drop_full_open(self):
        zeta = ResistanceCoefficient.from_Kv(self.Kvs, self.diameter)
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        D = self.diameter.to('m').magnitude
        A = math.pi * D ** 2 / 4
        v = V / A
        dp = zeta * rho * v ** 2 / 2
        self._pressure_drop = Q_(dp, 'Pa')

    def calculate_Kv_setting(self, extra_pressure_drop: Quantity) -> float:
        """Calculates and returns the Kv-setting of the balancing valve in the
        network that is required to compensate for the extra pressure drop that
        is needed to statically balance the pipe network.
        """
        dp_extra = extra_pressure_drop.to('Pa').magnitude
        dp_open = self._pressure_drop.to('Pa').magnitude
        dp_final = dp_open + dp_extra
        self._pressure_drop = Q_(dp_final, 'Pa')
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        Avr = V / math.sqrt(dp_final / rho)
        self.Kvr = FlowCoefficient.to_Kv(Avr)
        return self.Kvr

    @property
    def zeta(self) -> float:
        """Returns the corresponding resistance coefficient of the balancing
        valve.
        """
        if not math.isnan(self.Kvr):
            return ResistanceCoefficient.from_Kv(
                self.Kvr,
                self.pipe.cross_section.hydraulic_diameter
            )
        if not math.isnan(self.Kvs):
            return ResistanceCoefficient.from_Kv(
                self.Kvs,
                self.pipe.cross_section.hydraulic_diameter
            )
        return float('nan')

    @property
    def pressure_drop(self) -> Quantity:
        """Returns the pressure drop across the balancing valve."""
        return self._pressure_drop


class ControlValve(AbstractFitting):
    """Represents a control valve for regulating the flow."""

    def __init__(
        self,
        pipe: Pipe,
        ID: str = '',
        target_authority: float = 0.5,
        pressure_drop_crit_path: Quantity = Q_(float('nan'), 'Pa')
    ) -> None:
        """Creates a `ControlValve` object.

        Parameters
        ----------
        pipe:
            The pipe in which the balancing valve is to be installed.
        ID:
            Identifier for the balancing valve in the pipe network.
        target_authority:
            The initial valve authority targeted at. Normally it is recommended
            to target at a valve authority of around 0.5.
        pressure_drop_crit_path:
            The pressure loss along the critical path of the pipe network.
        """
        super().__init__(ID)
        self.pipe = pipe
        self.fluid = self.pipe.fluid
        self.diameter = self.pipe.cross_section.hydraulic_diameter
        self.volume_flow_rate = self.pipe.volume_flow_rate
        self.target_authority = target_authority
        self.pressure_drop_crit_path = pressure_drop_crit_path
        self._dp_full_open: Quantity = Q_(float('nan'), 'Pa')
        self.Kvs: float = float('nan')

    def calculate_preliminary_Kvs(self) -> float:
        """Calculates and returns the fully open flow coefficient Kvs of the
        control valve that is required to attain the valve authority targeted
        at.
        """
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        dp_crit_path = self.pressure_drop_crit_path.to('Pa').magnitude  # still wo control valve
        a = self.target_authority
        dp_valve_full_open = a * dp_crit_path / (1 - a)
        Avs = V / math.sqrt(dp_valve_full_open / rho)
        Kvs = FlowCoefficient.to_Kv(Avs)
        return Kvs

    def set_Kvs(self, Kvs: float):
        """Based on the preliminary Kvs-value, a commercially available control
        valve with a given Kvs-value can be selected. With this method
        the actual Kvs-value of the selected control valve can be set.
        """
        self.Kvs = Kvs
        self._calculate_pressure_drop_full_open()

    def _calculate_pressure_drop_full_open(self):
        """Calculates the pressure drop across the fully open control valve
        after the selected Kvs-value has been set.
        """
        zeta = ResistanceCoefficient.from_Kv(self.Kvs, self.diameter)
        rho = self.fluid.rho.to('kg / m ** 3').magnitude
        V = self.volume_flow_rate.to('m ** 3 / s').magnitude
        D = self.diameter.to('m').magnitude
        A = math.pi * D ** 2 / 4
        v = V / A
        dp = zeta * rho * v ** 2 / 2
        self._dp_full_open = Q_(dp, 'Pa')

    def get_valve_authority(self, pressure_drop_crit_path: Quantity) -> float:
        """Returns the actual valve authority of the control valve with the
        selected Kvs-value.
        """
        # to be called when the control valve has been added to the cross-over
        dp_valve_full_open = self._dp_full_open.to('Pa').magnitude
        dp_crit_path = pressure_drop_crit_path.to('Pa').magnitude
        return dp_valve_full_open / dp_crit_path

    @property
    def zeta(self) -> float:
        """Returns the resistance coefficient (float) of the fully open control
        valve.
        """
        return ResistanceCoefficient.from_Kv(
            self.Kvs,
            self.pipe.cross_section.hydraulic_diameter
        )

    @property
    def pressure_drop(self) -> Quantity:
        """Returns the pressure drop across the fully open control valve."""
        return self._dp_full_open

    def set_valve_opening(self, percent_open: int) -> float:
        """Returns the resistance coefficient (float) of the control valve that
        corresponds with a certain valve opening given as a percentage between
        fully closed (0 % open) and fully open (100 % open), assuming a linear
        relationship between the degree of valve opening and the flow
        coefficient of the valve.
        """
        Kv = (percent_open / 100) * self.Kvs
        zeta = ResistanceCoefficient.from_Kv(
            Kv,
            self.pipe.cross_section.hydraulic_diameter
        )
        return zeta
