from dataclasses import dataclass
import pandas as pd
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.energy_estimation import TimeSegment
from . import vrf

Q_ = Quantity


@dataclass
class DesignValues:
    """Data class containing the design heat load together with the associated
    outdoor and indoor temperature design value.

    Parameters
    ----------
    Ql: Quantity
        Design heating/cooling load
    Tia: Quantity
        Design indoor air temperature
    Toa: Quantity
        Design outdoor air temperature
    """
    Ql: Quantity
    Tia: Quantity
    Toa: Quantity


class Load:
    """Class that models the building heating load associated with the indoor
    air setpoint temperature active at the given time segment."""

    def __init__(
        self,
        time_segment: TimeSegment | None,
        design_values: DesignValues,
        Tia_avg: Quantity,
        Qhu: Quantity = Q_(0.0, 'kW')
    ) -> None:
        """Creates a `Load` instance.

        Parameters
        ----------
        time_segment: optional, default None
            The time segment that the given setpoint indoor air temperature
            is active.
        design_values:
            Instance of `DesignValues` holding the design heating load and
            the associated design indoor and outdoor temperatures. From this
            the heat transmission coefficient K of the building is determined.
        Tia_avg:
            The (average) indoor air setpoint temperature of the building that
            is supposed to remain constant during the given time segment.
        Qhu: default 0.0 kW
            The extra heating-up power supplied to the building during the
            given time segment.

        Notes
        -----
        To use the `Load` instance, you also need to set its attribute `Toa`,
        which is the outdoor air temperature, after its instantiation. Then you
        can use property `Q` of the `Load` instance to get the heating load of
        the building associated with `self.Tia_avg` and the current outdoor air
        temperature `Toa`.
        """
        self.time_segment = time_segment
        self.K = design_values.Ql / (design_values.Tia - design_values.Toa)
        self.Tia_avg = Tia_avg
        self.Qhu = Qhu
        self.Toa: Quantity | None = None

    @property
    def Q(self) -> Quantity:
        """Get the heating load of the building associated with `self.Tia_avg`
        and the current outdoor air temperature."""
        Q = self.K * (self.Tia_avg - self.Toa) + self.Qhu
        return max(Q_(0, 'W'), Q)

    def __call__(self, Toa: Quantity) -> Quantity:
        """Get the heating load of the building associated with `self.Tia_avg`
        and the given outdoor air temperature `Toa`."""
        self.Toa = Toa
        return self.Q


class _VRFSystem:
    """Underlying helper-class used inside class `EnergyConsumption` that
    encapsulates the 'real' `VRFSystem` instance and also the building load."""

    def __init__(self, vrf_system: vrf.VRFSystem):
        """Creates a `_VRFSystem` instance.

        Parameters
        ----------
        vrf_system:
            Instance of `VRFSystem` class, representing the VRF-system of a
            building.

        Notes
        -----
        The `_VRFSystem` instance has an attribute `load` through which a
        `Load`-instance is passed after instantiation of the `_VRFSystem`
        object.
        """
        self._vrf_system = vrf_system
        self.load: Load | None = None
        self._Toa: Quantity | None = None
        self._num_hours: Quantity | None = None

    @property
    def Toa(self) -> Quantity:
        """Get outdoor air temperature. In fact, this property (getter) is only
        implemented to have the setter."""
        return self._Toa

    @Toa.setter
    def Toa(self, v: Quantity):
        """Set the outdoor air temperature for the VRF-system and for the building
        load."""
        self._Toa = v
        self.load.Toa = v

    @property
    def num_hours(self) -> Quantity:
        """Get the number of hours in an outdoor air temperature bin. In fact,
        this property (getter) is only implemented to have the setter."""
        return self._num_hours

    @num_hours.setter
    def num_hours(self, v: int):
        """Set the number of hours in the current outdoor air temperature bin
        and time segment."""
        self._num_hours = Q_(v, 'hr')

    @property
    def Q_vrf(self) -> Quantity:
        """Get the available capacity of the VRF-system at the currently set
        outdoor air temperature."""
        Q_vrf = self._vrf_system.get_available_capacity(self.Toa)
        return Q_vrf

    @property
    def Q_aux(self) -> Quantity:
        """Get the additional capacity that is needed from an auxiliary heating
        source to compensate the building load at the given indoor
        air temperature and the currently set outdoor air temperature."""
        Q_aux = self.load.Q - self.Q_vrf
        if Q_aux > 0:
            return Q_aux
        else:
            return Q_(0.0, 'W')

    @property
    def W_in(self) -> Quantity:
        """Get the input power to the VRF-system at the given indoor
        air temperature and the currently set outdoor air temperature."""
        Tia_avg = self.load.Tia_avg
        Q_ou = self._vrf_system.get_full_load_outdoor_unit_capacity(Tia_avg, self.Toa)
        PLR = min(Q_(1.0, 'frac'), self.load.Q / Q_ou)
        W_in = self._vrf_system.get_input_power(Tia_avg, self.Toa, PLR)
        return W_in

    @property
    def E_in(self) -> Quantity:
        """Get the energy consumption of the VRF-system in the current outdoor
        air temperature bin and time segment."""
        return self.W_in * self.num_hours

    @property
    def E_aux(self) -> Quantity:
        """Get the energy consumption of the auxiliary heating source in the
        current outdoor air temperature bin and time segment."""
        return self.Q_aux * self.num_hours


class EnergyConsumption:
    """Class for calculating the energy consumption of a VRF-system using the
    bin-method."""

    def __init__(
        self,
        bin_table: pd.DataFrame,
        loads: list[Load],
        vrf_system: vrf.VRFSystem
    ) -> None:
        """Creates an `EnergyConsumption` instance.

        Parameters
        ----------
        bin_table:
            A bin table created with an instance of class `BinTableCreator`
            (see hvac.energy_consumption.bin_table.py).
        loads:
            A list of `Load` objects associated with the time segments of the
            bin table. Within each time segment of the day a different building
            load can be set depending on the indoor air setpoint temperature or
            the need for additional heating-up power.
        vrf_system:
            Instance of the `VRFSystem` class representing the VRF-system of
            the building.
        """
        self.bin_table = bin_table
        self.loads = loads
        self.vrf_system = _VRFSystem(vrf_system)
        self._E_consumption_data: list[list[float]] = []

    def estimate(self, T_unit: str = 'degC', E_unit: str = 'kWh') -> pd.DataFrame:
        """Estimates the energy consumption of the VRF-system and the auxiliary
        heating source within the given bin table and associated building
        loads.

        Parameters
        ----------
        T_unit: default 'degC'
            Temperature unit used in the bin table.
        E_unit: default 'kWh'
            The unit to be used to express energy consumption.

        Returns
        -------
        A Pandas `DataFrame` object with the energy consumption in each bin and
        time segment of the bin table.
        """
        for T_bin in self.bin_table.index:
            Toa = T_bin.mid
            r = []
            for j, time_segment in enumerate(self.bin_table.columns):
                load = self.loads[j]
                self.vrf_system.load = load
                self.vrf_system.Toa = Q_(Toa, T_unit)
                self.vrf_system.num_hours = self.bin_table.loc[T_bin, time_segment]
                E_vrf = self.vrf_system.E_in.to(E_unit).m
                E_aux = self.vrf_system.E_aux.to(E_unit).m
                E_tot = E_vrf + E_aux
                r.extend([E_vrf, E_aux, E_tot])
            self._E_consumption_data.append(r)
        day_periods = self.bin_table.columns
        consumption = ['VRF', 'Aux.', 'Tot.']
        columns = pd.MultiIndex.from_product(
            [day_periods, consumption],
            names=['Day Period', 'Consumption']
        )
        df = pd.DataFrame(
            data=self._E_consumption_data,
            index=self.bin_table.index,
            columns=columns
        )
        df.loc['TOTAL'] = df.sum()
        df_tmp = df.groupby(level=1, axis=1, sort=False).sum()
        df_tmp.columns = pd.MultiIndex.from_product([['TOTAL'], df_tmp.columns])
        df = df.join(df_tmp)
        return df


def find_balance_point(
    load: Load,
    vrf_system: vrf.VRFSystem,
    Toa_limits: tuple[Quantity, Quantity]
) -> tuple[Quantity, Quantity]:
    """Finds the balance point where building load line and VRF-system capacity
    curve intersect.

    Parameters
    ----------
    load:
        `Load` object representing the building load.
    vrf_system:
        `VRFSystem` instance representing the VRF-system.
    Toa_limits:
        Tuple with the lower and upper outdoor air temperature limits

    Returns
    -------
    Tuple with the 1st element the balancing outdoor air temperature and the
    2nd element the building load/VRF-system capacity in the balance point.
    """
    def eq(unknown: float) -> float:
        Toa = Q_(unknown, 'degC')
        Qc = vrf_system.get_available_capacity(Toa).to('kW').m
        Ql = load(Toa).to('kW').m
        out = Qc - Ql
        return out

    sol = root_scalar(
        eq,
        method='brentq',
        bracket=[Toa_limits[0].to('degC').m, Toa_limits[1].to('degC').m]
    )
    Toa_bal = Q_(sol.root, 'degC')
    Q_bal = load(Toa_bal)
    return Toa_bal, Q_bal
