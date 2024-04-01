from copy import deepcopy
import pandas as pd
from scipy import interpolate
from scipy import optimize
from hvac import Quantity
from .load import HeatingLoad

Q_ = Quantity


class HeatPump:
    default_units = {'T': 'degC', 'W_dot': 'kW', 'Q_dot': 'kW', 'COP': 'kW / kW'}

    def __init__(
        self,
        file_path: str,
        units: dict[str, str] | None = None,
        Cd: float = 0.25
    ) -> None:
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
            - Electric power: key 'W_dot' -> value 'kW'
            - Thermal power: key 'Q_dot' -> value 'kW'
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
        self._T_ext: Quantity | None = None
        self._load: HeatingLoad | None = None
        self._num_hours: Quantity = Q_(0, 'hr')

    def _read_file(self):
        df = pd.read_csv(self.file_path)
        self.T_ext_range = df.iloc[:, 0]
        self.Q_dot_range = df.iloc[:, 1]
        self.W_dot_range = df.iloc[:, 2]
        self.COP_range = self.Q_dot_range / self.W_dot_range
        self._Q_dot_fun = interpolate.interp1d(self.T_ext_range, self.Q_dot_range)
        self._W_dot_fun = interpolate.interp1d(self.T_ext_range, self.W_dot_range)
        self._COP_fun = interpolate.interp1d(self.T_ext_range, self.COP_range)

    @property
    def T_ext(self) -> Quantity:
        return self._T_ext

    @T_ext.setter
    def T_ext(self, q: Quantity) -> None:
        """Set outdoor air temperature."""
        self._T_ext = q.to('K')

    @property
    def load(self) -> HeatingLoad:
        return self._load

    @load.setter
    def load(self, o: HeatingLoad) -> None:
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
    def Q_dot(self) -> Quantity:
        """Get the steady-state available heating capacity of the heat pump at
        the set outdoor air temperature.
        """
        T_ext = self._T_ext.to(self.units['T']).m
        Q = self._Q_dot_fun(T_ext)
        return Q_(Q, self.units['Q_dot']).to('W')

    @property
    def W_dot(self) -> Quantity:
        """Get the steady-state input power the heat pump consumes when delivering
        the steady-state heating capacity at the set outdoor air temperature."""
        To = self._T_ext.to(self.units['T']).m
        W = self._W_dot_fun(To)
        return Q_(W, self.units['W_dot']).to('W')

    @property
    def COP(self) -> Quantity:
        """Get the steady-state COP of the heat pump at the set outdoor air
        temperature."""
        T_ext = self._T_ext.to(self.units['T']).m
        COP = self._COP_fun(T_ext)
        return Q_(COP, self.units['COP']).to('W / W')

    @property
    def balance_point(self) -> tuple[Quantity, Quantity]:
        """Get the balance point of the heat pump at the set load. Returns the
        outdoor air temperature and the corresponding heating capacity of the
        heat pump that matches with the heating load at this temperature.
        """
        _load = deepcopy(self._load)

        def eq(T_ext: float) -> float:
            T_ext = Q_(T_ext, 'K')
            _load.T_ext = T_ext
            self.T_ext = T_ext
            Q_load = _load.Q_dot_in
            Q_hp = self.Q_dot
            out = (Q_hp - Q_load).m
            return out

        T_ext_min = self._load.T_ext_min.m
        T_bal_load = self._load.T_bal.m
        T_bal_hp = optimize.brentq(eq, T_ext_min, T_bal_load)
        T_bal_hp = Q_(T_bal_hp, 'K')
        Q_bal_hp = self._Q_dot_fun(T_bal_hp.to(self.units['T']).m)
        Q_bal_hp = Q_(Q_bal_hp, self.units['Q_dot']).to('W')
        return T_bal_hp, Q_bal_hp

    @property
    def PLR(self) -> float:
        """Get the part-load ratio of the heat pump, i.e. the ratio of the
        heating load to the available, steady-state heating capacity."""
        PLR = min(
            1,
            (self._load.Q_dot_in.to('kW') / self.Q_dot.to('kW')).to('frac').m
        )
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
    def Q_dot_pl(self) -> Quantity:
        """Get the heat output of the heat pump at the current part-load
        condition."""
        Q_pl = min(self._load.Q_dot_in, self.Q_dot)
        return Q_pl

    @property
    def W_dot_pl(self) -> Quantity:
        """Get the power input taken up by the heat pump at the current part-load
        condition."""
        W_pl = self.Q_dot_pl / self.COP_pl
        return W_pl

    @property
    def Q_dot_aux(self) -> Quantity:
        """Get the required heat output of an auxiliary heater. Returns zero
        if the available heat pump capacity is sufficient to cover the set load
        (i.e. when the set outdoor air temperature is higher than the balance
        point temperature of the heat pump)."""
        Q_aux = self._load.Q_dot_in - self.Q_dot_pl
        if Q_aux > 0:
            return Q_aux
        else:
            return Q_(0, self.units['Q_dot'])

    @property
    def W(self) -> Quantity:
        """Get the amount of electric energy taken up by the heat pump at the set
        outdoor air temperature and the set number of hours that the given load
        is present.
        """
        W = self._num_hours * self.W_dot_pl
        return W

    @property
    def Q_aux(self) -> Quantity:
        """Get the amount of thermal energy that the auxiliary heater must deliver
        at the set outdoor air temperature and the set number of hours that the
        set load is present.
        """
        Q_aux = self._num_hours * self.Q_dot_aux
        return Q_aux
