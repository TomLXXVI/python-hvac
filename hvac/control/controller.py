import numpy as np
from hvac import Quantity

Q_ = Quantity


class PIDController:

    def __init__(
        self,
        SP: Quantity,
        Kc: Quantity,
        ti: Quantity,
        td: Quantity,
        bias: Quantity,
        PV_range: tuple[Quantity, Quantity],
        dt: Quantity
    ) -> None:
        """Creates a `PIDController` object.

        Parameters
        ----------
        SP:
            Setpoint of the controller.
        Kc:
            Controller's proportional gain factor in units of percent span.
        ti:
            Integration time.
        td:
            Derivative time.
        bias:
            Bias (controller output when deviation is zero) in percent
            (between 0 and 100 %).
        PV_range:
            Measuring range of the controller.
        dt:
            Sampling time of the controller.
        """
        self.SP = SP
        self.Kc = Kc
        self.ti = ti
        self.td = td
        self.bias = bias
        self.PV_min = PV_range[0]
        self.PV_max = PV_range[1]
        self.dt = dt

        self._SP = self.SP.to_base_units().magnitude
        self._Kc = self.Kc.to_base_units().magnitude
        self._ti = self.ti.to_base_units().magnitude
        self._td = self.td.to_base_units().magnitude
        self._bias = self.bias.to_base_units().magnitude
        self._PV_min = self.PV_min.to_base_units().magnitude
        self._PV_max = self.PV_max.to_base_units().magnitude
        self._dt = self.dt.to_base_units().magnitude

        self._pv = None
        self._e = 0.0
        self._I = 0.0
        self._D = 0.0
        self._reg_I = []
        self._reg_D = []

        # Set point as fraction of the controller's measuring range:
        self._sp = (self._SP - self._PV_min) / (self._PV_max - self._PV_min)

    def setpoint(self, SP: Quantity) -> None:
        """Changes the controller's setpoint."""
        self.SP = SP
        self._SP = self.SP.to_base_units().magnitude
        self._sp = (self._SP - self._PV_min) / (self._PV_max - self._PV_min)

    def _input(self, t: Quantity | float, PV: Quantity) -> None:
        """Passes the process value PV at time moment t to the controller."""
        t = t.to_base_units().magnitude if isinstance(t, Quantity) else t
        # If t is an integer multiple of the controller's sampling time, look
        # for new PV:
        if t % self._dt == 0:
            _PV = PV.to_base_units().magnitude
            if self._PV_min <= _PV <= self._PV_max:
                self._pv = (_PV - self._PV_min) / (self._PV_max - self._PV_min)
            elif _PV < self._PV_min:
                self._pv = 0.0
            else:
                self._pv = 1.0
        # Calculate the fractional error e
        self._e = self._pv - self._sp

    def _P_action(self) -> float:
        """Proportional control action."""
        return self._Kc * self._e

    def _I_action(self) -> float:
        """Integrating control action."""
        if self._ti != 0.0:
            if len(self._reg_I) < 3:
                # Add current deviation e to the I-register:
                self._reg_I.append(self._e)
            if len(self._reg_I) == 3:
                # When register contains three e-values, calculate integral
                # with Simpson's rule:
                _I = (
                    (self._reg_I[0] + 4 * self._reg_I[1] + self._reg_I[2])
                    * (self._dt / 3.0)
                )
                # Reset I-register:
                self._reg_I = [self._e]
                # Accumulated value of I-action:
                self._I += (self._Kc / self._ti) * _I
            return self._I
        else:
            return 0.0

    def _D_action(self):
        """Derivative control action."""
        if self._td != 0.0:
            if len(self._reg_D) < 3:
                # Add current deviation e to the D-register:
                self._reg_D.append(self._e)
            if len(self._reg_D) == 3:
                # When register contains three values, calculate derivative
                # de/dt using a backward finite difference approximation of
                # order O(h^2):
                _der_e = (
                    (self._reg_D[0] - 4 * self._reg_D[1] + 3 * self._reg_D[2])
                    / (2 * self._dt)
                )
                # Move the last two elements to the left in the D-register:
                self._reg_D = self._reg_D[-2:]
                # Value of D-action:
                self._D = self._Kc * self._td * _der_e
            return self._D
        else:
            return 0.0

    def _output(self) -> float:
        """Returns the control output as fraction of the control range."""
        p = self._P_action()
        i = self._I_action()
        d = self._D_action()
        out = p + i + d + self._bias
        if out < 0.0:
            out = 0.0
        elif out > 1.0:
            out = 1.0
        return out

    def __call__(self, t: Quantity | float, PV: Quantity) -> Quantity:
        """Passes the process value PV at time moment t to the controller and
        and returns the control output in percent of the control range.

        Parameters
        ----------
        t:
            Time passed since reference time (t_ref = 0).
        PV:
            The process value at time t.
        """
        self._input(t, PV)
        out = self._output()
        return Q_(out, 'frac').to('pct')

    def P_characteristic(self) -> tuple[Quantity, Quantity]:
        """Returns the percentage P-control characteristic. The abscissa is
        the control output in percent. The ordinate is the measured process
        value in percent.
        """
        pv = np.linspace(0.0, 1.0, endpoint=True)
        out = self.Kc * (pv - self._sp) + self._bias
        return Q_(out, 'frac').to('pct'), Q_(pv, 'frac').to('pct')

    def P_operating_point(self, Kp: Quantity, Z: Quantity) -> tuple[Quantity, Quantity]:
        """Returns the static operating point of the P-controller. This is the
        intersection of the static control characteristic and the static process
        characteristic for a constant value of the disturbance Z.

        Parameters
        ----------
        Kp:
            Percentage static gain of the process.
        Z:
            Constant value of the disturbance (with the same measuring unit as
            the process value PV)

        Returns
        -------
        The control output and process value in percent at the intersection.
        """
        _Z = Z.to_base_units.magnitude
        _Kp = Kp.to_base_units.magnitude
        _z = (_Z - self._PV_min) / (self._PV_max - self._PV_min)
        out = (self._Kc * (_z - self._sp) + self._bias) / (1 - self._Kc * _Kp)
        pv = (_z - _Kp * self._Kc * self._sp + _Kp * self._bias) / (1 - self._Kc * _Kp)
        return Q_(out, 'frac').to('pct'), Q_(pv, 'frac').to('pct')


class OnOffController:

    def __init__(
        self,
        SP: Quantity,
        HL_offset: Quantity,
        LL_offset: Quantity,
        PV_range: tuple[Quantity, Quantity],
        dt: Quantity,
        ctrl_dir: str = 'inverse'
    ) -> None:
        """Creates an `OnOffController` object.

        Parameters
        ----------
        SP:
            Setpoint.
        HL_offset:
            High limit offset of dead band with respect to SP.
        LL_offset:
            Low limit offset of dead band with respect to SP.
        PV_range:
            Measuring range of the controller.
        dt:
            Sampling time of the controller.
        ctrl_dir: {'inverse' (default), 'direct'}
            Control direction.
        """
        self.SP = SP
        self.HL_offset = HL_offset
        self.LL_offset = LL_offset
        self.PV_min = PV_range[0]
        self.PV_max = PV_range[1]
        self.dt = dt
        self.ctrl_dir = ctrl_dir

        self._SP = self.SP.to_base_units().magnitude
        self._HL = self._SP + self.HL_offset.to_base_units().m
        self._LL = self._SP + self.LL_offset.to_base_units().m
        self._PV_min = self.PV_min.to_base_units().m
        self._PV_max = self.PV_max.to_base_units().m
        self._dt = self.dt.to_base_units().m
        if self.ctrl_dir.lower() == 'inverse':
            self._ctrl_dir = -1
        else:
            self._ctrl_dir = 1

        self._sp = (self._SP - self._PV_min) / (self._PV_max - self._PV_min)
        self._hl = (self._HL - self._PV_min) / (self._PV_max - self._PV_min)
        self._ll = (self._LL - self._PV_min) / (self._PV_max - self._PV_min)
        self._e_hl = self._hl - self._sp
        self._e_ll = self._ll - self._sp

        self._pv = None
        self._e = None
        self._out = None

    def setpoint(self, SP: Quantity) -> None:
        """Changes the controller's setpoint."""
        self.SP = SP
        self._SP = self.SP.to_base_units().magnitude
        self._sp = (self._SP - self._PV_min) / (self._PV_max - self._PV_min)
        self._HL = self._SP + self.HL_offset.to_base_units().m
        self._LL = self._SP + self.LL_offset.to_base_units().m
        self._hl = (self._HL - self._PV_min) / (self._PV_max - self._PV_min)
        self._ll = (self._LL - self._PV_min) / (self._PV_max - self._PV_min)
        self._e_hl = self._hl - self._sp
        self._e_ll = self._ll - self._sp

    def _input(self, t: Quantity | float, PV: Quantity) -> None:
        """Passes the process value PV at time moment t to the controller."""
        t = t.to_base_units().magnitude if isinstance(t, Quantity) else t
        if t % self._dt == 0:
            _PV = PV.to_base_units().magnitude
            self._pv = (_PV - self._PV_min) / (self._PV_max - self._PV_min)
            self._e = self._pv - self._sp

    def _output(self):
        """Returns the control output as fraction of the control range."""
        if self._e <= self._e_ll:
            self._out = 1.0 if self._ctrl_dir == -1 else 0.0
        elif self._e >= self._e_hl:
            self._out = 0.0 if self._ctrl_dir == -1 else 0.0
        # if self._e is between self._e_ll and self._e_hl -> controller output
        # doesn't change.
        return self._out

    def __call__(self, t: Quantity | float, PV: Quantity) -> Quantity:
        """Passes the process value PV at time moment t to the controller and
        and returns the control output in percent of the control range.

        Parameters
        ----------
        t:
            Time passed since reference time (t_ref = 0).
        PV:
            The process value at time t.
        """
        self._input(t, PV)
        out = self._output()
        return Q_(out, 'frac').to('pct')


class PWMController:
    
    def __init__(
        self,
        SP: Quantity,
        Kc: Quantity,
        ti: Quantity,
        td: Quantity,
        bias: Quantity,
        n_cycle: int,
        PV_range: tuple[Quantity, Quantity],
        dt: Quantity
    ) -> None:
        """Creates a `PWMController` object.

        Parameters
        ----------
        SP:
            Setpoint of the controller.
        Kc:
            Controller's proportional gain factor in units of percent span.
        ti:
            Integration time.
        td:
            Derivative time.
        bias:
            Bias (controller output when deviation is zero) in percent
            (between 0 and 100 %).
        n_cycle:
            Number of sampling time steps that constitutes one PWM cycle period.
        PV_range:
            Measuring range of the controller.
        dt:
            Sampling time of the controller.
        """
        self._pid_controller = PIDController(SP, Kc, ti, td, bias, PV_range, dt)
        self._dt = dt.to_base_units().magnitude
        self.n_cycle = n_cycle
        self._T_cycle = self._dt * n_cycle
        self._cycle_output = []

    def setpoint(self, SP: Quantity):
        self._pid_controller.setpoint(SP)

    def __call__(self, t: Quantity | float, PV: Quantity) -> None:
        """Passes the process value PV at time moment t to the controller."""
        t = t.to_base_units().magnitude if isinstance(t, Quantity) else t
        if t % self._T_cycle == 0:
            # If time t is an integer multiple of the PWM cycle period, a PWM
            # output cycle has been completed:
            self._cycle_output = [0.0] * self.n_cycle
            out = self._pid_controller(t, PV)
            T_on = out.to('frac').m * self._T_cycle  # cycle ON time
            n_on = int(T_on // self._dt)             # number of ON time steps
            self._cycle_output[:n_on] = [1.0] * n_on
        if t % self._dt == 0.0:
            # If time t is an integer multiple of the controller's sampling time,
            # return and remove the first element from the PWM output cycle:
            return self._cycle_output.pop(0)
        else:
            # Only return (but don't remove) the first element:
            return self._cycle_output[0]
