from dataclasses import dataclass
import pandas as pd
from hvac import Quantity
from hvac.energy_estimation import TimeSegment
from . import vrf

Q_ = Quantity


@dataclass
class DesignValues:
    Ql: Quantity
    Tia: Quantity
    Toa: Quantity


class Load:

    def __init__(
        self,
        time_segment: TimeSegment,
        design_values: DesignValues,
        Tia_avg: Quantity,
        Qhu: Quantity = Q_(0.0, 'kW')
    ):
        self.time_segment = time_segment
        self.K = design_values.Ql / (design_values.Tia - design_values.Toa)
        self.Tia_avg = Tia_avg
        self.Qhu = Qhu
        self.Toa: Quantity | None = None

    @property
    def Q(self) -> Quantity:
        Q = self.K * (self.Tia_avg - self.Toa) + self.Qhu
        return max(Q_(0, 'W'), Q)


class _VRFSystem:

    def __init__(self, vrf_system: vrf.VRFSystem):
        self._vrf_system = vrf_system
        self.load: Load | None = None
        self._Toa: Quantity | None = None
        self._num_hours: Quantity | None = None

    @property
    def Toa(self) -> Quantity:
        return self._Toa

    @Toa.setter
    def Toa(self, v: Quantity):
        self._Toa = v
        self.load.Toa = v

    @property
    def num_hours(self) -> Quantity:
        return self._num_hours

    @num_hours.setter
    def num_hours(self, v: int):
        self._num_hours = Q_(v, 'hr')

    @property
    def Q_vrf(self) -> Quantity:
        Tia_avg = self.load.Tia_avg
        Q_vrf = self._vrf_system.get_available_capacity(Tia_avg, self.Toa)
        return Q_vrf

    @property
    def Q_aux(self) -> Quantity:
        Q_aux = self.load.Q - self.Q_vrf
        if Q_aux > 0:
            return Q_aux
        else:
            return Q_(0.0, 'W')

    @property
    def W_in(self) -> Quantity:
        Tia_avg = self.load.Tia_avg
        Q_ou = self._vrf_system.get_full_load_outdoor_unit_capacity(Tia_avg, self.Toa)
        PLR = self.load.Q / Q_ou
        W_in = self._vrf_system.get_input_power(Tia_avg, self.Toa, PLR)
        return W_in

    @property
    def E_in(self) -> Quantity:
        return self.W_in * self.num_hours

    @property
    def E_aux(self) -> Quantity:
        return self.Q_aux * self.num_hours


class EnergyConsumption:

    def __init__(
        self,
        bin_table: pd.DataFrame,
        loads: list[Load],
        vrf_system: vrf.VRFSystem
    ):
        self.bin_table = bin_table
        self.loads = loads
        self.vrf_system = _VRFSystem(vrf_system)
        self._E_consumption_data: list[list[float]] = []

    def estimate(self, T_unit: str = 'degC', E_unit: str = 'kWh') -> pd.DataFrame:
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
