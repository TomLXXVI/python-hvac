from typing import Optional
from dataclasses import dataclass
from .abc_process import Process, Equation, Variable
from hvac import Quantity
from hvac.fluids import HumidAir


Q_ = Quantity


class AirConditioningProcess(Process):
    
    def __init__(
        self,
        air_in: Optional[HumidAir] = None,
        T_ai: Optional[Quantity] = None,
        W_ai: Optional[Quantity] = None,
        h_ai: Optional[Quantity] = None,
        air_out: Optional[HumidAir] = None,
        T_ao: Optional[Quantity] = None,
        W_ao: Optional[Quantity] = None,
        h_ao: Optional[Quantity] = None,
        m_da: Optional[Quantity] = None,
        m_w: Optional[Quantity] = None,
        h_w: Optional[Quantity] = None,
        Q: Optional[Quantity] = None,
        Q_sen: Optional[Quantity] = None,
        Q_lat: Optional[Quantity] = None,
        SHR: Optional[Quantity] = None,
        ADP: Optional[HumidAir] = None,
        T_adp: Optional[Quantity] = None,
        W_adp: Optional[Quantity] = None,
        h_adp: Optional[Quantity] = None,
        beta: Optional[Quantity] = None
    ):
        """Creates an `AirConditioningProcess` object.

        Parameters
        ----------
        air_in: HumidAir
            state of entering air
        T_ai: Quantity
            dry-bulb temperature of entering air
        W_ai: Quantity
            humidity ratio of entering air
        h_ai: Quantity
            specific enthalpy of entering air
        air_out: HumidAir
            state of leaving air
        T_ao: Quantity
            dry-bulb temperature of leaving air
        W_ao: Quantity
            humidity ratio leaving air
        h_ao: Quantity
            specific enthalpy of leaving air
        m_da: Quantity
            mass flow rate of dry air
        m_w: Quantity
            mass flow rate of water entering or leaving the process
        h_w : Quantity
            specific enthalpy of water entering or leaving the process
        Q: Quantity
            total heat transfer in or out of the process
        Q_sen: Quantity
            sensible heat transfer in or out of the process
        Q_lat: Quantity
            latent transfer in or out of the process
        SHR: Quantity
            sensible heat ratio
        ADP: HumidAir
            apparatus dew point
        T_adp: Quantity
            dry-bulb temperature of ADP
        W_adp: Quantity
            humidity ratio of ADP
        h_adp: Quantity
            specific enthalpy of ADP
        beta : Quantity
            contact factor of air cooler or efficiency/effectiveness of humidifier/dehumidifier

        """
        # Set up the equations that describe air conditioning processes
        equations = {
            'heat_balance': Equation(
                variables=[
                    Variable(name='Q', unit='W'),
                    Variable(name='m_da', unit='kg / s'),
                    Variable(name='h_ao', unit='J / kg'),
                    Variable(name='h_ai', unit='J / kg'),
                    Variable(name='m_w', unit='kg / s'),
                    Variable(name='h_w', unit='J / kg')
                ],
                lhs="m_da * (h_ao - h_ai) - m_w * h_w - Q"
            ),
            'mass_balance': Equation(
                variables=[
                        Variable(name='m_da', unit='kg / s'),
                        Variable(name='W_ao', unit='kg / kg'),
                        Variable(name='W_ai', unit='kg / kg'),
                        Variable(name='m_w', unit='kg / s')
                ],
                lhs="m_da * (W_ao - W_ai) - m_w"
            ),
            'sensible_heat_balance': Equation(
                variables=[
                    Variable(name='m_da', unit='kg / s'),
                    Variable(name='T_ai', unit='K'),
                    Variable(name='T_ao', unit='K'),
                    Variable(name='Q_sen', unit='W'),
                    Variable(name='c_pa', value=1.02e3, unit='J / (kg * K)')
                    # mean specific heat of moist air between 0 and 60 °C
                ],
                lhs="m_da * c_pa * (T_ao - T_ai) - Q_sen"
            ),
            'latent_heat_balance': Equation(
                variables=[
                    Variable(name='m_da', unit='kg / s'),
                    Variable(name='W_ai', unit='kg / kg'),
                    Variable(name='W_ao', unit='kg / kg'),
                    Variable(name='h_w', unit='J / kg'),
                    Variable(name='Q_lat', unit='W'),
                    Variable(name='h_wg', value=2555e3, unit='J / kg')
                    # mean enthalpy of water vapor between 0 and 60 °C
                ],
                lhs="m_da * (W_ao - W_ai) * (h_wg - h_w) - Q_lat"
            ),
            'enthalpy-contactFactor-ADP': Equation(
                variables=[
                    Variable(name='h_ai', unit='J / kg'),
                    Variable(name='h_ao', unit='J / kg'),
                    Variable(name='h_adp', unit='J / kg'),
                    Variable(name='beta', unit='frac')  # contact factor
                ],
                lhs="(1 - beta) * h_ai + beta * h_adp - h_ao"
            ),
            'humidity_ratio-contactFactor-ADP': Equation(
                variables=[
                    Variable(name='W_ai', unit='kg / kg'),
                    Variable(name='W_ao', unit='kg / kg'),
                    Variable(name='W_adp', unit='kg / kg'),
                    Variable(name='beta', unit='frac')  # contact factor
                ],
                lhs="(1 - beta) * W_ai + beta * W_adp - W_ao"
            ),
            'dry_bulb_temperature-ContactFactor-ADP': Equation(
                variables=[
                    Variable(name='T_ai', unit='K'),
                    Variable(name='T_ao', unit='K'),
                    Variable(name='T_adp', unit='K'),
                    Variable(name='beta', unit='frac')  # contact factor
                ],
                lhs="(1 - beta) * T_ai + beta * T_adp - T_ao"
            ),
            'space_condition_line': Equation(
                variables=[
                    Variable(name='T_ai', unit='K'),
                    Variable(name='T_ao', unit='K'),
                    Variable(name='W_ai', unit='kg / kg'),
                    Variable(name='W_ao', unit='kg / kg'),
                    Variable(name='SHR', unit='frac'),
                    Variable(name='c_pa', value=1.02e3, unit='J / (kg * K)'),
                    # mean specific heat of moist air between 0 and 60 °C
                    Variable(name='h_wg', value=2555e3, unit='J / kg'),
                    # mean enthalpy of water vapor between 0 and 60 °C
                ],
                lhs="(W_ao - W_ai) / (T_ao - T_ai) - c_pa / h_wg * (1 / SHR - 1)"
            ),
            'total_sensible_latent_heat': Equation(
                variables=[
                    Variable(name='Q', unit='W'),
                    Variable(name='Q_sen', unit='W'),
                    Variable(name='Q_lat', unit='W')
                ],
                lhs="Q_sen + Q_lat - Q"
            ),
            'shr_sensible_latent_heat': Equation(
                variables=[
                    Variable(name='Q_lat', unit='W'),
                    Variable(name='Q_sen', unit='W'),
                    Variable(name='SHR', unit='frac')
                ],
                lhs="Q_lat - (1 - SHR) * Q_sen / SHR"
            ),
            'shr_sensible_total_heat': Equation(
                variables=[
                    Variable(name='Q', unit='W'),
                    Variable(name='SHR', unit='frac'),
                    Variable(name='Q_sen', unit='W')
                ],
                lhs="Q_sen - SHR * Q"
            )
        }
        super().__init__(equations)

        # Assign known input values (not None) of the air conditioning process
        # to the values of the corresponding variables in the equation.
        self._air_in = air_in
        if isinstance(self._air_in, HumidAir):
            self['T_ai'] = self._air_in.Tdb
            self['W_ai'] = self._air_in.W
            self['h_ai'] = self._air_in.h
        else:
            if T_ai is not None: self['T_ai'] = T_ai
            if W_ai is not None: self['W_ai'] = W_ai
            if h_ai is not None: self['h_ai'] = h_ai

        self._air_out = air_out
        if isinstance(self._air_out, HumidAir):
            self['T_ao'] = self._air_out.Tdb
            self['W_ao'] = self._air_out.W
            self['h_ao'] = self._air_out.h
        else:
            if T_ao is not None: self['T_ao'] = T_ao
            if W_ao is not None: self['W_ao'] = W_ao
            if h_ao is not None: self['h_ao'] = h_ao

        if m_da is not None: self['m_da'] = m_da

        if m_w is not None: self['m_w'] = m_w
        if h_w is not None: self['h_w'] = h_w

        if Q is not None: self['Q'] = Q
        if Q_sen is not None: self['Q_sen'] = Q_sen
        if Q_lat is not None: self['Q_lat'] = Q_lat

        self._ADP = ADP
        if isinstance(self._ADP, HumidAir):
            self['T_adp'] = self._ADP.Tdb
            self['W_adp'] = self._ADP.W
            self['h_adp'] = self._ADP.h
        else:
            if T_adp is not None:
                self['T_adp'] = T_adp
                self._ADP = HumidAir(Tdb=T_adp, RH=Q_(100, 'pct'))
                self['W_adp'] = self._ADP.W
                self['h_adp'] = self._ADP.h
            if W_adp is not None:
                self['W_adp'] = W_adp
                self._ADP = HumidAir(W=W_adp, RH=Q_(100, 'pct'))
                self['T_adp'] = self._ADP.Tdb
                self['h_adp'] = self._ADP.h
            if h_adp is not None:
                self['h_adp'] = h_adp
                self._ADP = HumidAir(h=h_adp, RH=Q_(100, 'pct'))
                self['T_adp'] = self._ADP.Tdb
                self['W_adp'] = self._ADP.W

        if beta is not None: self['beta'] = beta  # contact factor
        if SHR is not None: self['SHR'] = SHR

        # After assignment try to already solve for the inlet and outlet air state of the process and the ADP
        if self._air_in is None:
            try:
                _ = self.air_in
            except ValueError:
                pass

        if self._air_out is None:
            try:
                _ = self.air_out
            except ValueError:
                pass

        if self._ADP is None:
            try:
                _ = self.ADP
            except ValueError:
                pass

    @property
    def Q(self) -> Quantity:
        """Heat transfer rate to the air in the process."""
        return self['Q']

    @property
    def m_da(self) -> Quantity:
        """Mass flow rate of dry air through the process."""
        return self['m_da']

    @property
    def h_ai(self) -> Quantity:
        """Enthalpy of air at the inlet of the process."""
        return self['h_ai']

    @property
    def h_ao(self) -> Quantity:
        """Enthalpy of air at the outlet of the process"""
        return self['h_ao']

    @property
    def m_w(self) -> Quantity:
        """Mass flow rate of make-up water to or condensate from the process."""
        return self['m_w']

    @property
    def h_w(self) -> Quantity:
        """Enthalpy of make-up water to or condensate from the process."""
        return self['h_w']

    @property
    def Q_lat(self) -> Quantity:
        """Latent heat transfer rate to the air in the process."""
        return self['Q_lat']

    @property
    def W_ai(self) -> Quantity:
        """Humidity ratio of air at the inlet of the process."""
        return self['W_ai']

    @property
    def W_ao(self) -> Quantity:
        """Humidity ratio of air at the outlet of the process."""
        return self['W_ao']

    @property
    def Q_sen(self) -> Quantity:
        """Sensible heat transfer rate to the air in the process."""
        return self['Q_sen']

    @property
    def T_ai(self) -> Quantity:
        """Dry-bulb temperature of air at the inlet of the process."""
        return self['T_ai']

    @property
    def T_ao(self) -> Quantity:
        """Dry-bulb temperature of air at the outlet of the process."""
        return self['T_ao']

    @property
    def SHR(self) -> Quantity:
        """Sensible heat ratio of the process."""
        return self['SHR']

    @property
    def beta(self) -> Quantity:
        """Contact factor of air cooler or efficiency/effectiveness of
        humidifier/dehumidifier
        """
        return self['beta']

    @property
    def T_adp(self) -> Quantity:
        """Apparatus dew point temperature of the process."""
        return self['T_adp']

    @property
    def ADP(self) -> Optional[HumidAir]:
        """Air state at Apparatus Dew Point of the process."""
        MAX_ITER = 100
        if self._ADP is None:
            T_ai = self.T_ai.to('K').magnitude
            W_ai = self.W_ai.to('kg / kg').magnitude
            T_ao = self.T_ao.to('K').magnitude
            W_ao = self.W_ao.to('kg / kg').magnitude
            T_adp = self.air_out.Tdp.to('K').magnitude  # initial guess
            dT = T_ai - T_ao
            dW = W_ai - W_ao
            if dT == 0.0:
                return None
            if dW == 0.0:  # sensible process
                return HumidAir(Tdb=self.air_out.Tdp, W=self.air_out.W)
            i = 0
            while i < MAX_ITER:
                W_adp = W_ai - (dW / dT) * (T_ai - T_adp)
                ADP = HumidAir(Tdb=Q_(T_adp, 'K'), W=Q_(W_adp, 'kg / kg'))
                T_adp_new = ADP.Tdp.to('K').magnitude
                if abs(T_adp - T_adp_new) < 0.01:
                    self._ADP = HumidAir(Tdb=Q_(T_adp, 'K'), RH=Q_(100.0, 'pct'))
                    # update ADP variables in the register
                    self['T_adp'] = self._ADP.Tdb
                    self['W_adp'] = self._ADP.W
                    self['h_adp'] = self._ADP.h
                    break
                T_adp = T_adp_new
                i += 1
            else:
                raise ValueError('no solution found for ADP')
        return self._ADP

    @property
    def air_in(self) -> HumidAir:
        """Air state at the inlet of the process."""
        if self._air_in is None:
            try:
                self._air_in = HumidAir(Tdb=self.T_ai, W=self.W_ai)
                self['h_ai'] = self._air_in.h  # also update `h_ai` in the equations that contain this variable
            except ValueError:
                self._air_in = HumidAir(h=self.h_ai, W=self.W_ai)
                self['T_ai'] = self._air_in.Tdb  # also update `T_ai` in the equations that contain this variable
        return self._air_in

    @property
    def air_out(self) -> HumidAir:
        """Air state at the outlet of the process."""
        if self._air_out is None:
            try:
                self._air_out = HumidAir(Tdb=self.T_ao, W=self.W_ao)
                self['h_ao'] = self._air_out.h  # also update `h_ao` in the equations that contain this variable
            except ValueError:
                self._air_out = HumidAir(h=self.h_ao, W=self.W_ao)
                self['T_ao'] = self._air_out.Tdb  # also update `T_ao` in the equations that contain this variable
        return self._air_out


# noinspection PyUnresolvedReferences
@dataclass
class AirStream:
    """
    Class that groups the properties of an air stream. This class is used
    together with the `AdiabaticMixing` class.

    Parameters
    ----------
    state: HumidAir, optional
        thermodynamic state of the air in the stream
    m_da: Quantity, optional
        mass flow rate of dry air in the stream
    Tdb: Quantity, optional
        dry-bulb temperature of the air in the stream
    W: Quantity, optional
        humidity ratio of the air in the stream
    """
    state: Optional[HumidAir] = None
    m_da: Optional[Quantity] = None
    Tdb: Optional[Quantity] = None
    W: Optional[Quantity] = None
    h: Optional[Quantity] = None
    m_da_fractional: Optional[Quantity] = None
    
    def __post_init__(self):
        if isinstance(self.state, HumidAir):
            self.Tdb = self.state.Tdb
            self.W = self.state.W
            self.h = self.state.h


class AdiabaticMixing(Process):

    def __init__(
        self,
        in1: Optional[AirStream] = None,
        in2: Optional[AirStream] = None,
        out: Optional[AirStream] = None
    ) -> None:
        """Create `AdiabaticMixing` object.

        Parameters
        ----------
        in1: AirStream, optional
            first air stream entering the mixing process
        in2: AirStream, optional
            second air stream entering the mixing process
        out : AirStream, optional
            air stream coming out of the mixing process
        """
        # The equations that describe adiabatic air mixing
        equations = {
            'energy_balance': Equation(
                variables=[
                    Variable(name='h_ai1', unit='J / kg'),
                    Variable(name='h_ai2', unit='J / kg'),
                    Variable(name='h_ao', unit='J / kg'),
                    Variable(name='m_da_i1', unit='kg / s'),
                    Variable(name='m_da_i2', unit='kg / s'),
                    Variable(name='m_da_o', unit='kg / s')
                ],
                lhs="m_da_i1 * h_ai1 + m_da_i2 * h_ai2 - m_da_o * h_ao"
            ),
            'dry_air_mass_balance': Equation(
                variables=[
                    Variable(name='m_da_i1', unit='kg / s'),
                    Variable(name='m_da_i2', unit='kg / s'),
                    Variable(name='m_da_o', unit='kg / s')
                ],
                lhs="m_da_i1 + m_da_i2 - m_da_o"
            ),
            'vapor_mass_balance': Equation(
                variables=[
                    Variable(name='W_ai1', unit='kg / kg'),
                    Variable(name='W_ai2', unit='kg / kg'),
                    Variable(name='W_ao', unit='kg / kg'),
                    Variable(name='m_da_i1', unit='kg / s'),
                    Variable(name='m_da_i2', unit='kg / s'),
                    Variable(name='m_da_o', unit='kg / s')
                ],
                lhs="m_da_i1 * W_ai1 + m_da_i2 * W_ai2 - m_da_o * W_ao"
            ),
            'mass_flow_rate_entering_stream2': Equation(
                variables=[
                    Variable(name='W_ai1', unit='kg / kg'),
                    Variable(name='W_ai2', unit='kg / kg'),
                    Variable(name='W_ao', unit='kg / kg'),
                    Variable(name='m_da_i2', unit='kg / s'),
                    Variable(name='m_da_o', unit='kg / s')
                ],
                lhs="m_da_i2 * ((W_ao - W_ai2) / (W_ai1 - W_ao) + 1) - m_da_o"
            ),
            'mass_flow_rate_entering_stream1': Equation(
                variables=[
                    Variable(name='W_ai1', unit='kg / kg'),
                    Variable(name='W_ai2', unit='kg / kg'),
                    Variable(name='W_ao', unit='kg / kg'),
                    Variable(name='m_da_i1', unit='kg / s'),
                    Variable(name='m_da_o', unit='kg / s')
                ],
                lhs="m_da_i1 * (1 + (W_ai1 - W_ao) / (W_ao - W_ai2)) - m_da_o"
            ),
            'air_states': Equation(
                variables=[
                    Variable(name='W_ai1', unit='kg / kg'),
                    Variable(name='W_ai2', unit='kg / kg'),
                    Variable(name='W_ao', unit='kg / kg'),
                    Variable(name='T_ai1', unit='K'),
                    Variable(name='T_ai2', unit='K'),
                    Variable(name='T_ao', unit='K')
                ],
                lhs="W_ai1 + (W_ai2 - W_ai1) / (T_ai2 - T_ai1) * (T_ao - T_ai1) - W_ao"
            ),
            'fractional_energy_balance': Equation(
                variables=[
                    Variable(name='h_ai1', unit='J / kg'),
                    Variable(name='h_ai2', unit='J / kg'),
                    Variable(name='h_ao', unit='J / kg'),
                    Variable(name='x', unit='frac')
                ],
                lhs="x * h_ai1 + (1 - x) * h_ai2 - h_ao"
            ),
            'fractional_vapor_mass_balance': Equation(
                variables=[
                    Variable(name='W_ai1', unit='kg / kg'),
                    Variable(name='W_ai2', unit='kg / kg'),
                    Variable(name='W_ao', unit='kg / kg'),
                    Variable(name='x', unit='frac')
                ],
                lhs="x * W_ai1 + (1 - x) * W_ai2 - W_ao"
            )
        }
        super().__init__(equations)

        # Assign the input values of the adiabatic mixing process to the values of the variables in the equations
        self._in1 = in1 if in1 is not None else AirStream()
        self._in2 = in2 if in2 is not None else AirStream()
        self._out = out if out is not None else AirStream()

        if self._in1.m_da is not None: self['m_da_i1'] = self._in1.m_da
        if self._in2.m_da is not None: self['m_da_i2'] = self._in2.m_da
        if self._out.m_da is not None: self['m_da_o'] = self._out.m_da
        if self._in1.m_da_fractional is not None: self['x'] = self._in1.m_da_fractional
        if self._in2.m_da_fractional is not None: self['x'] = 1.0 - self._in2.m_da_fractional

        if isinstance(self._in1.state, HumidAir):
            self['T_ai1'] = self._in1.state.Tdb
            self['W_ai1'] = self._in1.state.W
            self['h_ai1'] = self._in1.state.h
        else:
            if self._in1.Tdb is not None: self['T_ai1'] = self._in1.Tdb
            if self._in1.W is not None: self['W_ai1'] = self._in1.W
            if self._in1.h is not None: self['h_ai1'] = self._in1.h

        if isinstance(self._in2.state, HumidAir):
            self['T_ai2'] = self._in2.state.Tdb
            self['W_ai2'] = self._in2.state.W
            self['h_ai2'] = self._in2.state.h
        else:
            if self._in2.Tdb is not None: self['T_ai2'] = self._in2.Tdb
            if self._in2.W is not None: self['W_ai2'] = self._in2.W
            if self._in2.h is not None: self['h_ai2'] = self._in2.h

        if isinstance(self._out.state, HumidAir):
            self['T_ao'] = self._out.state.Tdb
            self['W_ao'] = self._out.state.W
            self['h_ao'] = self._out.state.h
        else:
            if self._out.Tdb is not None: self['T_ao'] = self._out.Tdb
            if self._out.W is not None: self['W_ao'] = self._out.W
            if self._out.h is not None: self['h_ao'] = self._out.h

        # Try to solve for the air streams
        if self._in1.state is None or self._in1.m_da is None:
            try:
                _ = self.stream_in1
            except ValueError:
                pass

        if self._in2.state is None or self._in2.m_da is None:
            try:
                _ = self.stream_in2
            except ValueError:
                pass

        if self._out.state is None or self._out.m_da is None:
            try:
                _ = self.stream_out
            except ValueError:
                pass

    @property
    def m_da_o(self) -> Quantity:
        """Mass flow rate of the leaving air stream."""
        self._out.m_da = self['m_da_o']
        return self['m_da_o']

    @property
    def m_da_i1(self) -> Quantity:
        """Mass flow rate of the first entering air stream."""
        self._in1.m_da = self['m_da_i1']
        return self['m_da_i1']

    @property
    def m_da_i2(self) -> Quantity:
        """Mass flow rate of the second entering air stream."""
        self._in2.m_da = self['m_da_i2']
        return self['m_da_i2']

    @property
    def T_ai1(self) -> Quantity:
        """Dry-bulb temperature of the first entering air stream."""
        self._in1.Tdb = self['T_ai1']
        return self['T_ai1']

    @property
    def W_ai1(self) -> Quantity:
        """Humidity ratio of the first entering air stream."""
        self._in1.W = self['W_ai1']
        return self['W_ai1']

    @property
    def h_ai1(self) -> Quantity:
        """Enthalpy of the first entering air stream."""
        self._in1.h = self['h_ai1']
        return self['h_ai1']

    @property
    def T_ai2(self) -> Quantity:
        """Dry-bulb temperature of the second entering air stream."""
        self._in2.Tdb = self['T_ai2']
        return self['T_ai2']

    @property
    def W_ai2(self) -> Quantity:
        """Humidity ratio of the second entering air stream."""
        self._in2.W = self['W_ai2']
        return self['W_ai2']

    @property
    def h_ai2(self) -> Quantity:
        """Enthalpy of the second entering air stream."""
        self._in2.h = self['h_ai2']
        return self['h_ai2']

    @property
    def T_ao(self) -> Quantity:
        """Dry-bulb temperature of the leaving air stream."""
        self._out.Tdb = self['T_ao']
        return self['T_ao']

    @property
    def W_ao(self) -> Quantity:
        """Humidity ratio of the leaving air stream."""
        self._out.W = self['W_ao']
        return self['W_ao']

    @property
    def h_ao(self) -> Quantity:
        """Enthalpy of the leaving air stream."""
        self._out.h = self['h_ao']
        return self['h_ao']

    @property
    def stream_in1(self) -> AirStream:
        """Returns the state and mass flow rate of the first entering air
        stream inside an `AirStream` object."""
        try:
            self._in1.state = HumidAir(Tdb=self.T_ai1, W=self.W_ai1)
            self['h_ai1'] = self._in1.state.h
        except ValueError:
            self._in1.state = HumidAir(h=self.h_ai1, W=self.W_ai1)
            self['T_ai1'] = self._in1.state.Tdb
        self._in1.Tdb = self._in1.state.Tdb
        self._in1.W = self._in1.state.W
        self._in1.h = self._in1.state.h
        self._in1.m_da = self.m_da_i1
        return self._in1

    @property
    def stream_in2(self) -> AirStream:
        """Returns the state and mass flow rate of the second entering air
        stream inside an `AirStream` object."""
        try:
            self._in2.state = HumidAir(Tdb=self.T_ai2, W=self.W_ai2)
            self['h_ai2'] = self._in2.state.h
        except ValueError:
            self._in2.state = HumidAir(h=self.h_ai2, W=self.W_ai2)
            self['T_ai2'] = self._in2.state.Tdb
        self._in2.Tdb = self._in2.state.Tdb
        self._in2.W = self._in2.state.W
        self._in2.h = self._in2.state.h
        self._in2.m_da = self.m_da_i2
        return self._in2

    @property
    def stream_out(self) -> AirStream:
        """Returns the state and mass flow rate of the leaving air
        stream inside an `AirStream` object."""
        try:
            self._out.state = HumidAir(Tdb=self.T_ao, W=self.W_ao)
            self['h_ao'] = self._out.state.h
        except ValueError:
            # try:
            self._out.state = HumidAir(h=self.h_ao, W=self.W_ao)
            self['T_ao'] = self._out.state.Tdb
            # except ValueError:
            #     self._out.state = HumidAir(h=self.h_ao, RH=Q_(100.0, 'pct'))
            #     self['T_ao'] = self._out.state.Tdb
            #     self['h_ao'] = self._out.state.h
            #     self['W_ao'] = self._out.state.W
        self._out.Tdb = self._out.state.Tdb
        self._out.W = self._out.state.W
        self._out.h = self._out.state.h
        self._out.m_da = self.m_da_o
        return self._out

    @property
    def m_da_i1_fractional(self) -> Quantity:
        """Returns the mass flow rate of the first entering air stream as a
        fraction of the mass flow rate of the leaving air stream."""
        return self['x']

    @property
    def m_da_i2_fractional(self) -> Quantity:
        """Returns the mass flow rate of the second entering air stream as a
        fraction of the mass flow rate of the leaving air stream."""
        x = self.m_da_i1_fractional
        return 1 - x


class Fan(Process):

    def __init__(
        self,
        air_in: Optional[HumidAir] = None,
        air_out: Optional[HumidAir] = None,
        eta_fan: Optional[Quantity] = None,
        eta_motor: Optional[Quantity] = None,
        dP_fan: Optional[Quantity] = None,
        m_da: Optional[Quantity] = None
    ) -> None:
        """Creates `Fan` instance

        Parameters
        ----------
        air_in: optional
            The air state at the fan inlet.
        air_out: optional
            The air state at the fan outlet.
        eta_fan: optional
            Fan efficiency.
        eta_motor: optional
            Motor efficiency.
        dP_fan:
            The pressure difference across the fan at the known volume flow rate
            of air the fan delivers.
        m_da:
            The mass flow rate of air through the fan.
        """
        self.m_da = m_da

        # Set up the equation that describes heating of air due to fan inefficiency
        equations = {
            'fan_heating': Equation(
                variables=[
                    Variable(name='T_ai', unit='K'),
                    Variable(name='T_ao', unit='K'),
                    Variable(name='eta_f', unit='frac'),
                    Variable(name='eta_m', unit='frac'),
                    Variable(name='c_pa', unit='J / (kg * K)'),
                    Variable(name='rho_a', unit='kg / m ** 3'),
                    Variable(name='dp', unit='Pa'),
                ],
                lhs="T_ao - T_ai - dp / (eta_f * eta_m * c_pa * rho_a)"
            )
        }
        super().__init__(equations)

        # Assign known input values (not None) of the fan heating process to the
        # values of the corresponding variables in the equation.
        self._air_in = air_in
        if isinstance(self._air_in, HumidAir):
            self['T_ai'] = self._air_in.Tdb
        else:
            self['T_ai'] = None

        self._air_out = air_out
        if isinstance(self._air_out, HumidAir):
            self['T_ao'] = self._air_out.Tdb
        else:
            self['T_ao'] = None

        if eta_fan is not None: self['eta_f'] = eta_fan
        self['eta_m'] = Q_(1.0, 'frac') if eta_motor is None else eta_motor

        if isinstance(self._air_in, HumidAir):
            self['c_pa'] = self._air_in.cp
            self['rho_a'] = self._air_in.rho
        else:
            self['c_pa'] = self._air_out.cp
            self['rho_a'] = self._air_out.rho

        if dP_fan is not None: self['dp'] = dP_fan

    @property
    def air_out(self) -> HumidAir:
        """The air state at the fan outlet."""
        if self._air_out is None:
            T_ao = self['T_ao']
            W_ao = self._air_in.W
            # fan heating is a sensible process: humidity ratio = cst.
            self._air_out = HumidAir(Tdb=T_ao, W=W_ao)
        return self._air_out

    @property
    def air_in(self) -> HumidAir:
        """The air state at the fan outlet."""
        if self._air_in is None:
            T_ai = self['T_ai']
            W_ai = self._air_out.W
            self._air_in = HumidAir(Tdb=T_ai, W=W_ai)
        return self._air_in

    @property
    def Q(self) -> Quantity:
        """The heat given off by the fan to the air."""
        T_ai = self.air_in.Tdb.to('K')
        T_ao = self.air_out.Tdb.to('K')
        dT = T_ao - T_ai
        cp_ai = self.air_in.cp
        cp_ao = self.air_out.cp
        cp = (cp_ai + cp_ao) / 2
        Q = self.m_da * cp * dT
        return Q

    @property
    def W_input(self) -> Quantity:
        """The input power taken up by the fan."""
        dP = self['dp'].to('Pa')
        rho_ai = self.air_in.rho.to('kg / m ** 3')
        rho_ao = self.air_out.rho.to('kg / m ** 3')
        rho = (rho_ai + rho_ao) / 2
        W_fluid = self.m_da.to('kg / s') * dP / rho
        eta = self['eta_m'].to('frac') * self['eta_f'].to('frac')
        W_input = W_fluid / eta
        return W_input
