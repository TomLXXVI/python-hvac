from copy import deepcopy
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy import optimize
from hvac import Quantity
from .load import HeatingLoad

Q_ = Quantity


class HeatPump:
    """Represents an air-to-air or air-to-water heat pump in heating mode
    operation that may operate under varying heat source temperatures, while the
    heat sink temperature is kept constant.
    """
    units = {
        'T': 'degC',
        'W_dot': 'kW',
        'Q_dot': 'kW'
    }

    def __init__(
        self,
        file_path: str,
        units: dict[str, str] | None = None,
        Cd: float = 0.25
    ) -> None:
        """
        Creates a `HeatPump` object.

        Parameters
        ----------
        file_path: str
           File path to the csv-file with the heat pump performance table.
           This csv-file must have 3 columns:
            1.  First column must be the heat source temperature, i.e., the
                outdoor air temperature.
            2.  Second column must be the corresponding heating capacity of the
                heat pump for a constant heat sink temperature, i.e. a constant
                water or air leaving temperature.
            3.  Third column must be the corresponding electric input power to
                the heat pump.
        units: dict[str, str], optional (default None)
            Dictionary that specifies the units in which quantities are
            expressed in the heat pump performance file. If `None`, following
            units are assumed:
            - temperature: key 'T' -> value 'degC'
            - electric power: key 'W_dot' -> value 'kW'
            - thermal power: key 'Q_dot' -> value 'kW'
        Cd: float, default 0.25
            Degradation coefficient to account for degradation of the
            steady-state efficiency (COP) due to cycling at part-load.
            [See "Heating and Cooling of Buildings" (3rd Ed.) by T. Agami Reddy
            et al., chapter 14.7 (p. 407)].
        """
        self.file_path = file_path
        self.Cd = Cd
        if units is not None: self.units.update(units)
        self._read_file()
        self._T_ext: Quantity | None = None
        self._load: HeatingLoad | None = None
        self._num_hours: Quantity = Q_(0, 'hr')

    def _read_file(self):
        """Reads the heat pump performance file and creates interpolation
        functions that return heat capacity, input power and COP as function of
        heat source temperature.
        """
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
        """Returns the set outdoor air temperature."""
        return self._T_ext

    @T_ext.setter
    def T_ext(self, q: Quantity) -> None:
        """Sets the outdoor air temperature.
        If a heating load is attached to the heat pump, the outdoor air
        temperature is also set on this load.
        """
        self._T_ext = q.to('K')
        if self._load is not None:
            self._load.T_ext = self._T_ext

    @property
    def load(self) -> HeatingLoad:
        """Returns the heating load served by the heat pump."""
        return self._load

    @load.setter
    def load(self, o: HeatingLoad) -> None:
        """Attaches a `HeatingLoad` object to the `HeatPump` object. This is the
        heating load the heat pump must serve at the set outdoor air temperature.
        """
        self._load = o

    @property
    def num_hours(self) -> Quantity:
        """Returns the number of hours the heating pump is running with the set
        heating load. This is used to estimate energy consumption with the bin
        table method in class `EnergyEstimator` in module `energy_estimation.py`.
        """
        return self._num_hours

    @num_hours.setter
    def num_hours(self, v: int) -> None:
        """Sets the number of hours the heat pump is running with the set
        heating load. This is used to estimate energy consumption with the bin
        table method in class `EnergyEstimator` in module `energy_estimation.py`.
        If a heating load is attached to the heat pump, this number of hours is
        also set on this load.
        """
        self._num_hours = Q_(v, 'hr')
        if self._load is not None:
            self._load.num_hours = v

    @property
    def Q_dot(self) -> Quantity:
        """Returns the steady-state heating capacity of the heat pump at the set
        outdoor air temperature.
        """
        T_ext = self._T_ext.to(self.units['T']).m
        Q = self._Q_dot_fun(T_ext)
        return Q_(Q, self.units['Q_dot']).to('W')

    @property
    def W_dot(self) -> Quantity:
        """Returns the input power consumed by the heat pump in providing the
        heating capacity at the set outdoor air temperature under steady-state
        operation.
        """
        To = self._T_ext.to(self.units['T']).m
        W = self._W_dot_fun(To)
        return Q_(W, self.units['W_dot']).to('W')

    @property
    def COP(self) -> Quantity:
        """Returns the steady-state COP of the heat pump at the set outdoor air
        temperature.
        """
        T_ext = self._T_ext.to(self.units['T']).m
        COP = self._COP_fun(T_ext)
        unit_COP = self.units['Q_dot'] + " / " + self.units['W_dot']
        return Q_(COP, unit_COP).to('W / W')

    @property
    def balance_point(self) -> tuple[Quantity, Quantity]:
        """Returns the balance point of the heat pump with the set heating load.
        The coordinates of the balance point are the outdoor air temperature and
        the corresponding heating capacity that matches the heating load at this
        outdoor air temperature.
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
        """Returns the part-load ratio of the heat pump, i.e. the ratio of the
        set heating load to the available steady-state heating capacity of the
        heat pump.
        """
        PLR = min(
            1,
            (self._load.Q_dot_in.to('kW') / self.Q_dot.to('kW')).to('frac').m
        )
        return PLR

    @property
    def PLF(self) -> float:
        """Returns the calculated part-load factor of the heat pump, i.e. the
        ratio of the expected COP at part-load to the COP at steady-state
        operation of the heat pump. The PLF depends on the current part-load
        ratio of the heat pump and the degradation coefficient `Cd`.
        """
        PLF = 1 - self.Cd * (1 - self.PLR)
        return PLF

    @property
    def COP_pl(self) -> Quantity:
        """Returns the COP of the heat pump under the current part-load
        condition (i.e. current part-load ratio).
        """
        COP_pl = self.PLF * self.COP
        return COP_pl

    @property
    def Q_dot_pl(self) -> Quantity:
        """Returns the heat output of the heat pump at the current part-load
        condition.
        """
        Q_pl = min(self._load.Q_dot_in, self.Q_dot)
        return Q_pl

    @property
    def W_dot_pl(self) -> Quantity:
        """Returns the power input drawn by the heat pump at the current
        part-load condition.
        """
        W_pl = self.Q_dot_pl / self.COP_pl
        return W_pl

    @property
    def Q_dot_aux(self) -> Quantity:
        """Returns the required heat rate an auxiliary heater must provide to
        supplement the heat pump's heating capacity deficit. Returns zero when
        the available heat pump capacity is able to offset the set heating load
        (i.e. when the set outdoor air temperature is higher than the balance
        point temperature of the heat pump).
        """
        Q_aux = self._load.Q_dot_in - self.Q_dot_pl
        if Q_aux > 0:
            return Q_aux
        else:
            return Q_(0.0, 'W')

    @property
    def W(self) -> Quantity:
        """Returns the amount of electric energy drawn by the heat pump at the
        set outdoor air temperature during the number of hours the set heating
        load is present.
        """
        W = self._num_hours * self.W_dot_pl
        return W

    @property
    def Q_aux(self) -> Quantity:
        """Returns the amount of heating energy that the auxiliary heater must
        deliver at the set outdoor air temperature during the number of hours
        the set heating load is present.
        """
        Q_aux = self._num_hours * self.Q_dot_aux
        return Q_aux

    def get_heating_capacity_characteristic(self) -> tuple[Quantity, Quantity]:
        """Returns a range of outdoor temperatures and the corresponding heating
        capacity of the heat pump at these temperatures.
        """
        if self._load is not None:
            T_ext_min = self._load.T_ext_min.to('degC').m
            T_ext_max = self._load.T_int.to('degC').m
            T_ext_rng = Q_(np.linspace(T_ext_min, T_ext_max), 'degC')
        else:
            T_ext_rng = Q_(self.T_ext_range.to_numpy(), self.units['T']).to('degC')
        Q_dot_rng = []
        for T_ext in T_ext_rng:
            self.T_ext = T_ext
            Q_dot_rng.append(self.Q_dot)
        Q_dot_rng = Quantity.from_list(Q_dot_rng)
        return T_ext_rng, Q_dot_rng
