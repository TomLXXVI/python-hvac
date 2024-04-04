import warnings
import numpy as np
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from .bin_table import TimeSegment

warnings.filterwarnings('ignore', category=CoolPropWarning)

Q_ = Quantity
Air = Fluid('Air')


class HeatingLoad:
    """Represents the heating load of a building (or a space inside a building).
    The heating load is considered to be solely a function of outdoor and indoor
    temperature.
    """

    def __init__(
        self,
        T_int: Quantity,
        T_ext_min: Quantity,
        H_trm: Quantity,
        V_dot_ven: Quantity = Q_(0.0, 'm ** 3 / s'),
        Q_dot_ihg: Quantity = Q_(0.0, 'W'),
        eta_sys: Quantity = Q_(100, 'pct'),
        Q_dot_dhw: Quantity = Q_(0.0, 'W'),
        time_segment: TimeSegment | None = None
    ) -> None:
        """Creates a `HeatingLoad` object.

        Parameters
        ----------
        T_int: Quantity
            Indoor air temperature used in the heating load calculation of the
            building.
        T_ext_min: Quantity
            Minimal outdoor air temperature used to calculate the heating load
            under design conditions.
        H_trm: Quantity
            Building's transmission heat loss coefficient derived from the
            building's heating load calculation.
        V_dot_ven: Quantity
            Ventilation/infiltration air volume flow rate.
        Q_dot_ihg: Quantity
            The building's internal heat gain.
        eta_sys: Quantity
            Heating system efficiency, i.e., the ratio of the system's heat
            output to the system's heat input. This takes any heat energy losses
            into account between the point where heat energy is generated and
            the point(s) where heat energy is delivered to the building.
        Q_dot_dhw: Quantity, optional
            (Average) heating power for domestic hot water preparation.
        time_segment: TimeSegment, default None
            The time period of the day the load is present. (This parameter is
            used in class `EnergyEstimator` in module `energy_estimation.py`.)
        """
        self.time_segment = time_segment
        self.T_int = T_int.to('K')
        self.H_trm = H_trm.to('W / K')
        self.V_dot_ven = V_dot_ven.to('m**3 / s')
        self.Q_dot_ihg = Q_dot_ihg.to('W')
        self.Q_dot_dhw = Q_dot_dhw.to('W')
        self.eta_sys = eta_sys.to('frac')
        self.T_ext_min = T_ext_min.to('K')
        self._outdoor_air = Air(T=T_ext_min, P=Q_(101_325, 'Pa'))
        self._num_hours: Quantity = Q_(0, 'hr')

    @property
    def T_ext(self) -> Quantity:
        """Returns the current set outdoor temperature."""
        return self._outdoor_air.T.to('K')

    @T_ext.setter
    def T_ext(self, q: Quantity) -> None:
        """Sets the current outdoor air temperature."""
        self._outdoor_air = Air(T=q, P=Q_(101_325, 'Pa'))

    @property
    def num_hours(self) -> Quantity:
        """Returns the number of hours the heating load is present.
        This is used to estimate energy consumption with the bin table method
        in class `EnergyEstimator` in module `energy_estimation.py`.
        """
        return self._num_hours

    @num_hours.setter
    def num_hours(self, v: int) -> None:
        """Sets the number of hours the load is present.
        This is used to estimate energy consumption with the bin table method
        in class `EnergyEstimator` in module `energy_estimation.py`.
        """
        self._num_hours = Q_(v, 'hr')

    @property
    def H_ven(self) -> Quantity:
        """Returns the heat loss coefficient due to air ventilation or
        infiltration.
        """
        H_ven = self._outdoor_air.rho * self._outdoor_air.cp * self.V_dot_ven
        return H_ven.to('W / K')

    @property
    def H_tot(self) -> Quantity:
        """Returns the total heat loss coefficient of the building including
        heat loss due to transmission through the building envelope and heat
        loss due to air ventilation/infiltration.
        It is the sum of transmission heat loss coefficient and the ventilation/
        infiltration heat loss coefficient.
        """
        H_tot = self.H_trm + self.H_ven
        return H_tot.to('W / K')

    @property
    def T_bal(self) -> Quantity:
        """Returns the balance temperature of the building.
        This is the outdoor temperature at which the internal heat gain in
        the building balances the total heat loss of the building. It depends
        on the size of the internal heat gain and the set indoor temperature.
        """
        T_bal = self.T_int - self.Q_dot_ihg / self.H_tot
        return T_bal.to('K')

    @property
    def Q_dot_out(self) -> Quantity:
        """Returns the required heating power that must be emitted in the
        building to maintain the set indoor air temperature when the set
        outdoor air temperature is present. If there is a domestic hot water
        load, it is also included.
        """
        Q_dot_out = self.H_tot * (self.T_bal - self.T_ext) + self.Q_dot_dhw
        Q_dot_out = max(Q_(0.0, 'W'), Q_dot_out.to('W'))
        return Q_dot_out

    @property
    def Q_dot_trm(self) -> Quantity:
        """Returns the heat loss due to transmission through the building's
        envelope at the set outdoor and indoor temperature.
        """
        Q_dot_trm = self.H_trm * (self.T_int - self.T_ext)
        Q_dot_trm = max(Q_(0.0, 'W'), Q_dot_trm.to('W'))
        return Q_dot_trm

    @property
    def Q_dot_ven(self) -> Quantity:
        """Returns the heat loss due to air ventilation and/or air infiltration
        at the set outdoor and indoor temperature.
        """
        Q_dot_ven = self.H_ven * (self.T_int - self.T_ext)
        Q_dot_ven = max(Q_(0.0, 'W'), Q_dot_ven.to('W'))
        return Q_dot_ven

    @property
    def Q_dot_in(self) -> Quantity:
        """Returns the thermal power that the heating system must deliver to the
        building to maintain the set indoor temperature at the set outdoor
        temperature. This also takes the efficiency of the heating system into
        account.
        """
        Q_dot_in = self.Q_dot_out / self.eta_sys
        return Q_dot_in.to('W')

    @property
    def Q_in(self) -> Quantity:
        """Returns the amount of energy the heating system must deliver to
        maintain the set indoor temperature at the set outdoor air temperature
        depending on the set number of hours that the load is present.
        This is used to estimate energy consumption with the bin table method
        in class `EnergyEstimator` in module `energy_estimation.py`.
        """
        Q_in = self._num_hours * self.Q_dot_in
        return Q_in.to('Wh')

    def get_characteristic(self) -> tuple[Quantity, Quantity]:
        """Returns a range of outdoor temperatures between the set minimum and
        the set indoor temperature and the corresponding heat rates demanded
        from the heating system.
        """
        T_ext_min = self.T_ext_min.to('degC').m
        T_ext_max = self.T_int.to('degC').m
        T_ext_rng = Q_(np.linspace(T_ext_min, T_ext_max), 'degC')
        Q_dot_in_rng = []
        for T_ext in T_ext_rng:
            self.T_ext = T_ext
            Q_dot_in_rng.append(self.Q_dot_in)
        Q_dot_in_rng = Quantity.from_list(Q_dot_in_rng)
        return T_ext_rng, Q_dot_in_rng
