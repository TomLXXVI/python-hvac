import pandas as pd
from hvac import Quantity
from .load import HeatingLoad
from .heat_pump import HeatPump

Q_ = Quantity


class EnergyEstimator:

    def __init__(
        self,
        bin_table: pd.DataFrame,
        loads: list[HeatingLoad],
        heat_pump: HeatPump | None = None
    ) -> None:
        """
        Create and initialize `EnergyEstimator` object.

        Parameters
        ----------
        bin_table: Pandas-DataFrame
            Monthly bin table retrieved from a `BinTableCreator` object.
        loads: List[Load]
            A list of `Load` objects in the same order of the time segments
            used in the bin table. They represent the heating load of the
            building for each period of the day.
        heat_pump: HeatPump, optional
            Heat pump of which the energy consumption must be estimated. If
            not specified or None, the required thermal energy to maintain
            indoor air temperature inside the building is estimated.
        """
        self.bin_table = bin_table
        self.loads = loads
        self.heat_pump = heat_pump
        self._E_load_table: list[list[HeatingLoad]] = []

    def _estimate(
        self,
        T_unit: str = 'degC',
        E_unit: str = 'kWh'
    ) -> pd.DataFrame:
        for T_bin in self.bin_table.index:
            T_ext = T_bin.mid
            r = []
            for j, time_segment in enumerate(self.bin_table.columns):
                load = self.loads[j]
                load.T_ext = Q_(T_ext, T_unit)
                load.num_hours = self.bin_table.loc[T_bin, time_segment]
                r.append(load.Q_in.to(E_unit).m)
            self._E_load_table.append(r)
        df = pd.DataFrame(
            data=self._E_load_table,
            index=self.bin_table.index,
            columns=self.bin_table.columns
        )
        df.loc['TOTAL'] = df.sum()
        df.loc[:, 'TOTAL'] = df.sum(axis=1)
        return df

    def _estimate_with_heat_pump(
        self,
        T_unit: str = 'degC',
        E_unit: str = 'kWh'
    ) -> pd.DataFrame:
        for T_bin in self.bin_table.index:
            T_ext = T_bin.mid
            r = []
            for j, time_segment in enumerate(self.bin_table.columns):
                load = self.loads[j]
                load.T_ext = Q_(T_ext, T_unit)
                self.heat_pump.load = load
                self.heat_pump.T_ext = Q_(T_ext, T_unit)
                self.heat_pump.num_hours = self.bin_table.loc[T_bin, time_segment]
                W_hp = self.heat_pump.W.to(E_unit).m
                Q_aux = self.heat_pump.Q_aux.to(E_unit).m
                E_tot = W_hp + Q_aux
                r.extend([W_hp, Q_aux, E_tot])
            self._E_load_table.append(r)
        day_periods = self.bin_table.columns
        consumption = ['HP', 'Aux.', 'Total']
        columns = pd.MultiIndex.from_product(
            [day_periods, consumption],
            names=['Day Period', 'Consumption']
        )
        df = pd.DataFrame(
            data=self._E_load_table,
            index=self.bin_table.index,
            columns=columns
        )
        # take the sum of all the columns and add them in a new row TOTAL at the
        # bottom of the dataframe
        df.loc['TOTAL'] = df.sum()
        # take the sum of all the rows
        df_tmp = df.T.groupby(level=1, sort=False).sum()
        df_tmp = df_tmp.T
        # give the temporary dataframe the same multi-column index as the
        # original dataframe
        df_tmp.columns = pd.MultiIndex.from_product([['TOTAL'], df_tmp.columns])
        # join the temporary dataframe with the original one
        df = df.join(df_tmp)
        return df

    def estimate(self, T_unit: str = 'degC', E_unit: str = 'kWh') -> pd.DataFrame:
        """Returns a Pandas-DataFrame that is based on the given bin table with
        the estimated required thermal energy for the building or with the
        consumed electric energy of the heat pump filled in.
        """
        if self.heat_pump is not None:
            return self._estimate_with_heat_pump(T_unit, E_unit)
        else:
            return self._estimate(T_unit, E_unit)
