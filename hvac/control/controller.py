import numpy as np
from hvac import Quantity

Q_ = Quantity


class PIDController:

    def __init__(
        self,
        setpoint: Quantity,
        K_ctrl: float,
        t_integral: Quantity,
        t_deriv: Quantity,
        bias: Quantity,
        measuring_range: tuple[Quantity, Quantity],
        dt: Quantity
    ) -> None:
        """Creates a `PIDController` object.

        Parameters
        ----------
        setpoint:
            Controller setpoint.
        K_ctrl:
            Proportional gain.
        t_integral:
            Integration time.
        t_deriv:
            Derivative time.
        bias:
            Bias (i.e. the controller's output when the deviation is zero) in
            percent (between 0 and 100 %).
        measuring_range:
            Measuring range of the controller.
        dt:
            Sampling time of the controller.
        """
        self.setpoint = setpoint
        self.K_ctrl = K_ctrl
        self.t_integral = t_integral
        self.t_deriv = t_deriv
        self.bias = bias
        self.PV_min = measuring_range[0]
        self.PV_max = measuring_range[1]
        self.dt = dt

        self._SP = self.setpoint.to_base_units().magnitude
        self._ti = self.t_integral.to_base_units().magnitude
        self._td = self.t_deriv.to_base_units().magnitude
        self._B = self.bias.to_base_units().magnitude
        self._PV_min = self.PV_min.to_base_units().magnitude
        self._PV_max = self.PV_max.to_base_units().magnitude
        self._dt = self.dt.to_base_units().magnitude

        self._pv = None
        self._e = 0.0
        self._I = 0.0
        self._D = 0.0
        self._I_register = []
        self._D_register = []

        # Setpoint as a fraction of the controller's measuring range:
        self._sp = (self._SP - self._PV_min) / (self._PV_max - self._PV_min)

    def setpoint(self, value: Quantity) -> None:
        """Changes the setpoint."""
        self.setpoint = value
        self._SP = self.setpoint.to_base_units().magnitude
        self._sp = (self._SP - self._PV_min) / (self._PV_max - self._PV_min)

    def __input__(self, t: Quantity | float, PV: Quantity) -> None:
        """Sets the process value `PV` at time moment `t`."""
        t = t.to_base_units().magnitude if isinstance(t, Quantity) else t
        # If `t` is an integer multiple of the controller's sampling time,
        # accept the new PV:
        if t % self._dt == 0:
            _PV = PV.to_base_units().magnitude
            if self._PV_min <= _PV <= self._PV_max:
                self._pv = (_PV - self._PV_min) / (self._PV_max - self._PV_min)
            elif _PV < self._PV_min:
                self._pv = 0.0
            else:
                self._pv = 1.0
        # Calculate the relative error `e`:
        self._e = self._pv - self._sp

    def __P_action__(self) -> float:
        return self.K_ctrl * self._e

    def __I_action__(self) -> float:
        if self._ti != 0.0:
            if len(self._I_register) < 3:
                # Add current deviation `e` to the I-register:
                self._I_register.append(self._e)
            if len(self._I_register) == 3:
                # When the register contains three values, calculate the integral
                # using Simpson's rule:
                e0 = self._I_register[0]
                e1 = self._I_register[1]
                e2 = self._I_register[2]
                _I = (e0 + 4 * e1 + e2) * (self._dt / 3.0)
                # Reset the register:
                self._I_register = [self._e]
                # Accumulated value:
                self._I += (self.K_ctrl / self._ti) * _I
            return self._I
        else:
            return 0.0

    def __D_action__(self) -> float:
        if self._td != 0.0:
            if len(self._D_register) < 3:
                # Add current deviation `e` to the D-register:
                self._D_register.append(self._e)
            if len(self._D_register) == 3:
                # When the register contains three values, calculate the
                # derivative using a backward finite difference approximation
                # of order O(h²):
                e0 = self._D_register[0]
                e1 = self._D_register[1]
                e2 = self._D_register[2]
                _der_e = (e0 - 4 * e1 + 3 * e2) / (2 * self._dt)
                # Keep the last two elements in the D-register:
                self._D_register = self._D_register[-2:]
                # Value of D-action:
                self._D = self.K_ctrl * self._td * _der_e
            return self._D
        else:
            return 0.0

    def __output__(self) -> float:
        """Returns the controller's output as a percentage of the control
        range.
        """
        p = self.__P_action__()
        i = self.__I_action__()
        d = self.__D_action__()
        out = p + i + d + self._B
        if out < 0.0:
            out = 0.0
        elif out > 1.0:
            out = 1.0
        return out

    def __call__(self, t: Quantity | float, PV: Quantity) -> Quantity:
        """Sets the process value `PV` at time `t` and returns the controller's
        output as a percentage of the control range.

        Parameters
        ----------
        t:
            Time measured from t = 0.
        PV:
            Process value at time t.
        """
        self.__input__(t, PV)
        out = self.__output__()
        return Q_(out, 'frac').to('pct')

    def P_characteristic(self) -> tuple[Quantity, Quantity]:
        """Returns the percent P-control characteristic. The abscissa is
        the controller's output in percent. The ordinate is the process value in
        percent.
        """
        pv = np.linspace(0.0, 1.0, endpoint=True)
        out = self.K_ctrl * (pv - self._sp) + self._B
        return Q_(out, 'frac').to('pct'), Q_(pv, 'frac').to('pct')

    def P_operating_point(self, K_proc: float, Z: Quantity) -> tuple[Quantity, Quantity]:
        """Returns the static operating point of the P-controller, being the
        intersection of the static control characteristic and the static process
        characteristic for a constant value of the disturbance `Z`.

        Parameters
        ----------
        K_proc:
            Static gain of the process as a percentage.
        Z:
            Constant value of the disturbance (having the same measuring unit as
            the process value `PV`).

        Returns
        -------
        The controller's output and the process value as a percentage.
        """
        _Z = Z.to_base_units().magnitude
        _z = (_Z - self._PV_min) / (self._PV_max - self._PV_min)
        out = (self.K_ctrl * (_z - self._sp) + self._B) / (1 - self.K_ctrl * K_proc)
        pv = (_z - K_proc * self.K_ctrl * self._sp + K_proc * self._B) / (1 - self.K_ctrl * K_proc)
        return Q_(out, 'frac').to('pct'), Q_(pv, 'frac').to('pct')


class OnOffController:

    def __init__(
        self,
        setpoint: Quantity,
        upper_limit: Quantity,
        lower_limit: Quantity,
        measuring_range: tuple[Quantity, Quantity],
        dt: Quantity,
        control_action: str = 'reverse'
    ) -> None:
        """Creates an `OnOffController` object.

        Parameters
        ----------
        setpoint:
        upper_limit:
            Dead band upper limit with respect to SP.
        lower_limit:
            Lower limit of dead band with respect to SP.
        measuring_range:
            Measuring range of the controller.
        dt:
            Sampling time of the controller.
        control_action: {'reverse' (default), 'direct'}
            Control action.
        """
        self.setpoint = setpoint
        self.upper_limit = upper_limit
        self.lower_limit = lower_limit
        self.PV_min = measuring_range[0]
        self.PV_max = measuring_range[1]
        self.dt = dt
        self.ctrl_action = control_action

        self._SP = self.setpoint.to_base_units().magnitude
        self._H_LIM = self._SP + self.upper_limit.to_base_units().m
        self._L_LIM = self._SP + self.lower_limit.to_base_units().m
        self._PV_min = self.PV_min.to_base_units().m
        self._PV_max = self.PV_max.to_base_units().m
        self._dt = self.dt.to_base_units().m
        if self.ctrl_action.lower() == 'reverse':
            self._sign = -1
        else:
            self._sign = 1

        self._sp = (self._SP - self._PV_min) / (self._PV_max - self._PV_min)
        self._h_lim = (self._H_LIM - self._PV_min) / (self._PV_max - self._PV_min)
        self._l_lim = (self._L_LIM - self._PV_min) / (self._PV_max - self._PV_min)
        self._e_hl = self._h_lim - self._sp
        self._e_ll = self._l_lim - self._sp

        self._pv = None
        self._e = None
        self._out = None

    def setpoint(self, SP: Quantity) -> None:
        """Changes the setpoint."""
        self.setpoint = SP
        self._SP = self.setpoint.to_base_units().magnitude
        self._sp = (self._SP - self._PV_min) / (self._PV_max - self._PV_min)
        self._H_LIM = self._SP + self.upper_limit.to_base_units().m
        self._L_LIM = self._SP + self.lower_limit.to_base_units().m
        self._h_lim = (self._H_LIM - self._PV_min) / (self._PV_max - self._PV_min)
        self._l_lim = (self._L_LIM - self._PV_min) / (self._PV_max - self._PV_min)
        self._e_hl = self._h_lim - self._sp
        self._e_ll = self._l_lim - self._sp

    def __input__(self, t: Quantity | float, PV: Quantity) -> None:
        """Gets the process value `PV` at time moment `t`."""
        t = t.to_base_units().magnitude if isinstance(t, Quantity) else t
        if t % self._dt == 0:
            _PV = PV.to_base_units().magnitude
            self._pv = (_PV - self._PV_min) / (self._PV_max - self._PV_min)
            self._e = self._pv - self._sp

    def __output__(self) -> float:
        """Returns the controller's output as a fraction of the control range."""
        if self._e <= self._e_ll:
            self._out = 1.0 if self._sign == -1 else 0.0
        elif self._e >= self._e_hl:
            self._out = 0.0 if self._sign == -1 else 1.0
        # elif self._e_ll < self._e < self._e_hl --> no change.
        return self._out

    def __call__(self, t: Quantity | float, PV: Quantity) -> Quantity:
        """Sets the process value `PV` at time moment `t` and returns the
        controller's output as a percentage of the control range. In case
        of an on/off-controller this is either 0% (off), or 100% (on).

        Parameters
        ----------
        t:
            Time measured from t = 0.
        PV:
            The process value at time t.
        """
        self.__input__(t, PV)
        out = self.__output__()
        return Q_(out, 'frac').to('pct')


class PWMController:
    
    def __init__(
        self,
        setpoint: Quantity,
        K_ctrl: float,
        t_integral: Quantity,
        t_deriv: Quantity,
        bias: Quantity,
        measuring_range: tuple[Quantity, Quantity],
        dt: Quantity,
        n_cycle: int
    ) -> None:
        """Creates a `PWMController` object.

        Parameters
        ----------
        setpoint:
            Controller's setpoint.
        K_ctrl:
            Controller's proportional gain factor in units of percent span.
        t_integral:
            Integration time.
        t_deriv:
            Derivative time.
        bias:
            Bias (i.e. the controller's output when the deviation is zero) as a
            percentage (between 0 and 100 %).
        measuring_range:
            Controller's measuring range.
        dt:
            Sampling time of the controller.
        n_cycle:
            Number of sampling time steps that constitutes the period of the
            PWM-cycle.
        """
        self._pid_controller = PIDController(
            setpoint, K_ctrl,
            t_integral, t_deriv,
            bias, measuring_range, dt
        )
        self._dt = dt.to_base_units().magnitude
        self.n_cycle = n_cycle
        self._T_cycle = self._dt * n_cycle
        self._cycle_output = []

    def setpoint(self, SP: Quantity):
        self._pid_controller.setpoint(SP)

    def __call__(self, t: Quantity | float, PV: Quantity) -> Quantity:
        """Sets the process value `PV` at time moment `t` and returns the
        controller's output as a percentage of the control range.
        """
        t = t.to_base_units().magnitude if isinstance(t, Quantity) else t
        if t % self._T_cycle == 0:
            # If time `t` is an integer multiple of the PWM-cycle period, a PWM
            # output cycle has been completed:
            self._cycle_output = [0.0] * self.n_cycle
            out = self._pid_controller(t, PV)
            T_on = out.to('frac').m * self._T_cycle  # cycle ON-time
            n_on = int(T_on // self._dt)             # number of ON-time steps
            self._cycle_output[:n_on] = [1.0] * n_on
        if t % self._dt == 0.0:
            # If time `t` is an integer multiple of the controller's sampling time,
            # return and remove the first element from the PWM-output cycle:
            return Q_(self._cycle_output.pop(0), 'frac').to('pct')
        else:
            # Only return (but don't remove) the first element:
            return Q_(self._cycle_output[0], 'frac').to('pct')
