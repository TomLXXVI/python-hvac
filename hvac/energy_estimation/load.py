import warnings
from hvac import Quantity
from hvac.fluids import Fluid, CoolPropWarning
from .bin_table import TimeSegment

warnings.filterwarnings('ignore', category=CoolPropWarning)

Q_ = Quantity
Air = Fluid('Air')


class Load:

    def __init__(
        self,
        T_int: Quantity,
        T_ext_min: Quantity,
        H_trm: Quantity,
        V_dot_ven: Quantity = Q_(0.0, 'm ** 3 / s'),
        Q_dot_ihg: Quantity = Q_(0.0, 'W'),
        eta_sys: Quantity = Q_(100, 'pct'),
        time_segment: TimeSegment | None = None
    ) -> None:
        """Creates a `Load` object.

        Parameters
        ----------
        T_int: Quantity
            Indoor air temperature.
        T_ext_min: Quantity
            Minimal outdoor air temperature (design value).
        H_trm: Quantity
            Building's transmission heat loss coefficient derived from the
            building's heat load calculation.
        V_dot_ven: Quantity
            Ventilation/infiltration air volume flow rate.
        Q_dot_ihg: Quantity
            The building's internal heat gain.
        eta_sys: Quantity
            System efficiency, i.e., the ratio of the system's heat output to
            the system's heat input.
        time_segment: TimeSegment, default None
            The time period of the day the load is present.
        """
        self.time_segment = time_segment
        self.T_int = T_int
        self.H_trm = H_trm
        self.V_dot_ven = V_dot_ven
        self.Q_dot_ihg = Q_dot_ihg
        self.eta_sys = eta_sys
        self.T_ext_min = T_ext_min
        self._outdoor_air = Air(T=T_ext_min, P=Q_(101_325, 'Pa'))
        self._num_hours: Quantity = Q_(0, 'hr')

    @property
    def T_ext(self) -> Quantity:
        return self._outdoor_air.T

    @T_ext.setter
    def T_ext(self, q: Quantity) -> None:
        """Sets the outdoor air temperature."""
        self._outdoor_air = Air(T=q, P=Q_(101_325, 'Pa'))

    @property
    def num_hours(self) -> Quantity:
        return self._num_hours

    @num_hours.setter
    def num_hours(self, v: int) -> None:
        """Sets the number of hours the load is present."""
        self._num_hours = Q_(v, 'hr')

    @property
    def H_ven(self) -> Quantity:
        """Returns the heat loss coefficient due to air ventilation or
        infiltration.
        """
        H_ven = self._outdoor_air.rho * self._outdoor_air.cp * self.V_dot_ven
        return H_ven

    @property
    def H_tot(self) -> Quantity:
        """Returns the total heat loss coefficient of the building including
        heat loss due to transmission through the building envelope and heat
        loss due to air ventilation/infiltration.
        """
        H_tot = self.H_trm + self.H_ven
        return H_tot

    @property
    def T_bal(self) -> Quantity:
        """Returns the balance temperature of the building."""
        T_bal = self.T_int - self.Q_dot_ihg / self.H_tot
        return T_bal

    @property
    def Q_dot_out(self) -> Quantity:
        """Returns the required thermal power that must be emitted in the
        building to maintain the desired indoor air temperature at the set
        outdoor air temperature.
        """
        Qe = self.H_tot * (self.T_bal - self.T_ext)
        return Qe

    @property
    def Q_dot_trm(self) -> Quantity:
        """Returns the heat loss due to transmission through the building
        envelope at the set outdoor air temperature.
        """
        Q_dot_trm = self.H_trm * (self.T_int - self.T_ext)
        return Q_dot_trm

    @property
    def Q_dot_ven(self) -> Quantity:
        """Returns the heat loss due to air ventilation and/or air infiltration
        at the set outdoor air temperature.
        """
        Q_dot_ven = self.H_ven * (self.T_int - self.T_ext)
        return Q_dot_ven

    @property
    def Q_dot_in(self) -> Quantity:
        """Returns the required thermal power that the heating system must
        deliver to maintain the desired indoor air temperature at the set outdoor
        air temperature.
        """
        Q_dot_in = self.Q_dot_out / self.eta_sys
        return Q_dot_in

    @property
    def Q_in(self) -> Quantity:
        """Returns the required amount of thermal energy that the heating system
        must deliver to maintain the desired indoor air temperature at the set
        outdoor air temperature and the set number of hours that the load is
        present.
        """
        Q_in = self._num_hours * self.Q_dot_in
        return Q_in
