from hvac import Quantity
from hvac.fluids import HumidAir
from .bin_table import TimeSegment

Q_ = Quantity


class Load:

    def __init__(
        self,
        Ti: Quantity,
        K_tr: Quantity,
        V_vi: Quantity = Q_(0.0, 'm ** 3 / s'),
        Qig: Quantity = Q_(0.0, 'W'),
        eta_inst: Quantity = Q_(100, 'pct'),
        outdoor_air: HumidAir = HumidAir(Tdb=Q_(0, 'degC'), RH=Q_(0, 'pct')),
        To_min: Quantity | None = None,
        time_segment: TimeSegment | None = None
    ):
        """
        Create a `Load` object.

        Parameters
        ----------
        Ti: Quantity
            Indoor air temperature.
        K_tr: Quantity
            Building heat loss coefficient by transmission that results
            from the heat load calculation of the building.
        V_vi: Quantity
            Ventilation volume flow rate and/or infiltration air volume flow rate.
        Qig: Quantity
            Internal heat gains.
        eta_inst: Quantity
            Installation efficiency.
        outdoor_air: HumidAir
            Outdoor air state.
        To_min: Quantity
            Minimal outdoor air temperature (design value).
        time_segment: TimeSegment, default None
            The time period of the day the load is present.
        """
        self.time_segment = time_segment
        self.Ti = Ti
        self.K_tr = K_tr
        self.V_vi = V_vi
        self.Qig = Qig
        self.eta_inst = eta_inst
        self.To_min = To_min
        self._outdoor_air = outdoor_air
        self._num_hours: Quantity = Q_(0, 'hr')

    @property
    def To(self) -> Quantity:
        return self._outdoor_air.Tdb

    @To.setter
    def To(self, q: Quantity) -> None:
        """Set outdoor air temperature."""
        self._outdoor_air = HumidAir(Tdb=q, RH=self._outdoor_air.RH)

    @property
    def outdoor_air(self) -> HumidAir:
        return self._outdoor_air

    @outdoor_air.setter
    def outdoor_air(self, oa: HumidAir) -> None:
        """Set outdoor air state."""
        self._outdoor_air = oa

    @property
    def num_hours(self) -> Quantity:
        return self._num_hours

    @num_hours.setter
    def num_hours(self, v: int) -> None:
        """Set number of hours the load is present."""
        self._num_hours = Q_(v, 'hr')

    @property
    def K_vi(self) -> Quantity:
        """Get the heat loss coefficient due to air ventilation or infiltration."""
        K_vi = self._outdoor_air.rho * self._outdoor_air.cp * self.V_vi
        return K_vi

    @property
    def K(self) -> Quantity:
        """Get the total heat loss coefficient of the building including heat loss
        due to transmission through the building envelope and heat loss due
        to air ventilation/infiltration."""
        K = self.K_tr + self.K_vi
        return K

    @property
    def T_bal(self) -> Quantity:
        """Get the balance temperature of the building."""
        T_bal = self.Ti - self.Qig / self.K
        return T_bal

    @property
    def Qe(self) -> Quantity:
        """Get the required thermal power that must be emitted in the building
        to maintain the desired indoor air temperature at the set outdoor
        air temperature."""
        Qe = self.K * (self.T_bal - self.To)
        return Qe

    @property
    def Q_loss_tr(self) -> Quantity:
        """Get the heat loss due to transmission through the building
        envelope at the set outdoor air temperature."""
        Ql_tr = self.K_tr * (self.Ti - self.To)
        return Ql_tr

    @property
    def Q_loss_vi(self) -> Quantity:
        """Get the heat loss due to air ventilation and/or air infiltration at
        the set outdoor air temperature."""
        Ql_vi = self.K_vi * (self.Ti - self.To)
        return Ql_vi

    @property
    def Q(self) -> Quantity:
        """Get the required thermal power that the heating system must deliver
        to maintain the desired indoor air temperature at the set outdoor
        air temperature."""
        Q = self.Qe / self.eta_inst
        return Q

    @property
    def E(self) -> Quantity:
        """Get the required amount of thermal energy that the heating system must
        deliver to maintain the desired indoor air temperature at the set outdoor
        air temperature and the set number of hours that the load is present."""
        E = self._num_hours * self.Q
        return E
