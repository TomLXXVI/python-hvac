from typing import Dict, Tuple
import pandas as pd
from scipy import interpolate
from scipy import optimize
from hvac import Quantity
from .load import Load

Q_ = Quantity


class HeatPump:
    default_units = {'T': 'degC', 'W': 'kW', 'Q': 'kW', 'COP': 'kW / kW'}

    def __init__(self, file_path: str, units: Dict[str, str] | None = None, Cd: float = 0.25):
        """
        Create `HeatPump` object.

        Parameters
        ----------
        file_path: str
            The file path to the csv-file with the heat pump performance data
            at varying outdoor air temperature. The csv-file has 3 columns:
            1.  First column must be outdoor air temperature.
            2.  Second column must be the corresponding heating capacity.
            3.  Third column must be the corresponding electric input power.
        units: Dict[str, str], optional (default None)
            Dictionary that specifies the units in which the quantities are
            expressed in the heat pump performance datafile. If not specified or
            None the default units of the class are assumed:
            - Temperature: key 'T' -> value 'degC'
            - Electric power: key 'W' -> value 'kW'
            - Thermal power: key 'Q' -> value 'kW'
            - COP: key 'COP' -> value 'kW / kW'
        Cd: float, default 0.25
            Degradation coefficient to account for degradation of steady-state
            efficiency (COP) due to cycling at part-load [see "Heating and Cooling
            of Buildings" (3rd Ed.) by T. Agami Reddy et al. §14.7 (p.407)]
        """
        self.file_path = file_path
        self.units = units if units is not None else self.default_units
        self.Cd = Cd
        self._read_file()
        self._To: Quantity | None = None
        self._load: Load | None = None
        self._num_hours: Quantity = Q_(0, 'hr')

    def _read_file(self):
        df = pd.read_csv(self.file_path)
        self.To_range = df.iloc[:, 0]
        self.Q_range = df.iloc[:, 1]
        self.W_range = df.iloc[:, 2]
        self.COP_range = self.Q_range / self.W_range
        self._Q_fun = interpolate.interp1d(self.To_range, self.Q_range)
        self._W_fun = interpolate.interp1d(self.To_range, self.W_range)
        self._COP_fun = interpolate.interp1d(self.To_range, self.COP_range)

    @property
    def To(self) -> Quantity:
        return self._To

    @To.setter
    def To(self, q: Quantity) -> None:
        """Set outdoor air temperature."""
        self._To = q

    @property
    def load(self) -> Load:
        return self._load

    @load.setter
    def load(self, o: Load) -> None:
        """Set `Load` object. This is the heating load the heat pump must
        compensate for at the set outdoor air temperature."""
        self._load = o

    @property
    def num_hours(self) -> Quantity:
        return self._num_hours

    @num_hours.setter
    def num_hours(self, v: int) -> None:
        """Set the number of hours the heat pump must compensate for the set
        load."""
        self._num_hours = Q_(v, 'hr')

    @property
    def Q(self) -> Quantity:
        """Get the steady-state available heating capacity of the heat pump at
        the set outdoor air temperature.
        """
        To = self._To.to(self.units['T']).m
        Q = self._Q_fun(To)
        return Q_(Q, self.units['Q'])

    @property
    def W(self) -> Quantity:
        """Get the steady-state input power the heat pump consumes when delivering
        the steady-state heating capacity at the set outdoor air temperature."""
        To = self._To.to(self.units['T']).m
        W = self._W_fun(To)
        return Q_(W, self.units['W'])

    @property
    def COP(self) -> Quantity:
        """Get the steady-state COP of the heat pump at the set outdoor air
        temperature."""
        To = self._To.to(self.units['T']).m
        COP = self._COP_fun(To)
        return Q_(COP, self.units['COP'])

    @property
    def balance_point(self) -> Tuple[Quantity, ...]:
        """Get the balance point of the heat pump at the set load. Returns the
        outdoor air temperature and the corresponding heating capacity of the
        heat pump that matches with the heating load at this temperature.
        """
        def eq(To: float) -> float:
            To = Quantity(To, self.units['T'])
            self._load.To = To
            self.To = To
            Q_load = self._load.Q
            Q_hp = self.Q
            out = (Q_hp - Q_load).to(self.units['Q']).m
            return out

        To_min = self._load.To_min.to(self.units['T']).m
        T_bal_load = self._load.T_bal.to(self.units['T']).m
        T_bal_hp = optimize.brentq(eq, To_min, T_bal_load)
        T_bal_hp = Q_(T_bal_hp, self.units['T'])
        self.To = T_bal_hp
        Q_bal_hp = self.Q
        return T_bal_hp, Q_bal_hp

    @property
    def PLR(self) -> float:
        """Get the part-load ratio of the heat pump, i.e. the ratio of the
        heating load to the available, steady-state heating capacity."""
        PLR = min(1, (self._load.Q / self.Q).to('frac').m)
        return PLR

    @property
    def PLF(self) -> float:
        """Get the calculated part-load factor of the heat pump, i.e. the ratio
        of the expected COP at part-load to the COP at steady-state operation
        of the heat pump. It depends on the current part-load ratio of the
        heat pump.
        """
        PLF = 1 - self.Cd * (1 - self.PLR)
        return PLF

    @property
    def COP_pl(self) -> Quantity:
        """Get the COP of the heat pump at the current part-load condition
        (i.e. current part-load ratio)"""
        COP_pl = self.PLF * self.COP
        return COP_pl

    @property
    def Q_pl(self) -> Quantity:
        """Get the heat output of the heat pump at the current part-load
        condition."""
        Q_pl = min(self._load.Q, self.Q)
        return Q_pl

    @property
    def W_pl(self) -> Quantity:
        """Get the power input taken up by the heat pump at the current part-load
        condition."""
        W_pl = self.Q_pl / self.COP_pl
        return W_pl

    @property
    def Q_aux(self) -> Quantity:
        """Get the required heat output of an auxiliary heater. Returns zero
        if the available heat pump capacity is sufficient to cover the set load
        (i.e. when the set outdoor air temperature is higher than the balance
        point temperature of the heat pump)."""
        Q_aux = self._load.Q - self.Q_pl
        if Q_aux > 0:
            return Q_aux
        else:
            return Q_(0, self.units['Q'])

    @property
    def E(self) -> Quantity:
        """Get the amount of electric energy taken up by the heat pump at the set
        outdoor air temperature and the set number of hours that the given load
        is present.
        """
        E = self._num_hours * self.W_pl
        return E

    @property
    def E_aux(self) -> Quantity:
        """Get the amount of thermal energy that the auxiliary heater must deliver
        at the set outdoor air temperature and the set number of hours that the
        set load is present.
        """
        E_aux = self._num_hours * self.Q_aux
        return E_aux
