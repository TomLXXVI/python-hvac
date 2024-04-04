from copy import deepcopy
import numpy as np
import pandas as pd
from hvac import Quantity
from hvac.charts import BarChart
from .load import HeatingLoad
from .heat_pump import HeatPump

Q_ = Quantity


class EnergyEstimator:
    """
    Estimation of the energy need of a building or estimation of the electrical
    energy consumption of an air-to-water/air-to-air heat pump system using the
    bin table method.
    """
    def __init__(
        self,
        bin_table: pd.DataFrame,
        loads: list[HeatingLoad] | None = None,
        heat_pumps: list[HeatPump] | None = None
    ) -> None:
        """Creates and initializes the `EnergyEstimator` object.

        Parameters
        ----------
        bin_table: Pandas-DataFrame
            Monthly bin table retrieved from a `BinTableCreator` object.
        loads: list[HeatingLoad], optional
            A list of `HeatingLoad` objects having the same order as the time
            segments used in the bin table. They represent the different heating
            loads of the building during each time segment of the day.
        heat_pumps: list[HeatPump], optional
            A list of `HeatPump` objects with an attached `HeatingLoad` object
            having the same order as the time segments used in the bin table.
            They represent the different heat pump system operating states
            during each time segment of the day.

        If parameter `loads` is specified, the energy need of the building is
        estimated. If parameter `heat_pumps` is specified, the electrical energy
        consumption of the heat pump system is estimated.
        """
        if loads is None and heat_pumps is None:
            raise ValueError(
                "Either parameter `loads` or parameter `heat_pumps` must be set."
            )
        self.loads = loads
        self.heat_pumps = heat_pumps
        self.bin_table = bin_table
        self._load_table: list[list[HeatingLoad]] = []
        self._output: pd.DataFrame | None = None
        self._T_unit: str = ''
        self._E_unit: str = ''

    def _estimate(self) -> pd.DataFrame:
        """Estimates the energy need of a building."""
        for T_bin in self.bin_table.index:
            T_ext = T_bin.mid
            r = []
            for j, time_segment in enumerate(self.bin_table.columns):
                load = self.loads[j]
                load.T_ext = Q_(T_ext, self._T_unit)
                load.num_hours = self.bin_table.loc[T_bin, time_segment]
                r.append(load.Q_in.to(self._E_unit).m)
            self._load_table.append(r)
        df = pd.DataFrame(
            data=self._load_table,
            index=self.bin_table.index,
            columns=self.bin_table.columns
        )
        df.loc['TOTAL'] = df.sum()
        df.loc[:, 'TOTAL'] = df.sum(axis=1)
        return df

    def _estimate_with_heat_pump(self) -> pd.DataFrame:
        """Estimates the energy consumption of the heat pump system."""
        for T_bin in self.bin_table.index:
            T_ext = T_bin.mid
            r = []
            for j, time_segment in enumerate(self.bin_table.columns):
                heat_pump = self.heat_pumps[j]
                heat_pump.T_ext = Q_(T_ext, self._T_unit)
                heat_pump.num_hours = self.bin_table.loc[T_bin, time_segment]
                W_hp = heat_pump.W.to(self._E_unit).m
                Q_aux = heat_pump.Q_aux.to(self._E_unit).m
                E_tot = W_hp + Q_aux
                r.extend([W_hp, Q_aux, E_tot])
            self._load_table.append(r)
        day_periods = self.bin_table.columns
        consumption = ['HP', 'Aux.', 'Total']
        columns = pd.MultiIndex.from_product(
            [day_periods, consumption],
            names=['Day Period', 'Consumption']
        )
        df = pd.DataFrame(
            data=self._load_table,
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

    def estimate(
        self,
        T_unit: str = 'degC',
        E_unit: str = 'kWh'
    ) -> pd.DataFrame:
        """Returns a Pandas-DataFrame that is based on the given bin table with
        the estimated required thermal energy for the building or with the
        consumed electric energy of the heat pump filled in.

        Parameters
        ----------
        T_unit: str, default 'degC'
            The unit in which temperature values are expressed.
        E_unit: str, default 'kWh'
            The unit in which energy values are expressed.
        """
        self._T_unit, self._E_unit = T_unit, E_unit
        if self.heat_pumps is not None:
            self._output = self._estimate_with_heat_pump()
            return self._output
        elif self.loads is not None:
            self._output = self._estimate()
            return self._output

    def _draw_bar_chart(self) -> BarChart:
        if self._output is not None:
            T_ext_rng = [bin_.mid for bin_ in self._output.index[:-1]]
            chart = BarChart()
            bottom = np.zeros(len(T_ext_rng))
            for column in self._output.columns[:-1]:
                chart.add_xy_data(
                    label=column,
                    x1_values=T_ext_rng,
                    y1_values=self._output[column][:-1].values,
                    style_props={'width': 1.0, 'bottom': deepcopy(bottom)}
                )
                bottom += self._output[column][:-1].values
            chart.add_legend(columns=4)
            T_unit = f"{Q_(0, self._T_unit).units:~P.0f}"
            E_unit = f"{Q_(0, self._E_unit).units:~P.0f}"
            chart.x1.add_title(f'outdoor temperature, {T_unit}')
            chart.y1.add_title(f'energy, {E_unit}')
            return chart

    def _draw_bar_chart_with_heat_pump(self) -> BarChart:
        if self._output is not None:
            T_ext_rng = [bin_.mid for bin_ in self._output.index[:-1]]
            df = self._output['TOTAL']
            E_hp_rng = df['HP'][:-1].values
            E_aux_rng = df['Aux.'][:-1].values
            chart = BarChart()
            chart.add_xy_data(
                label='HP',
                x1_values=T_ext_rng,
                y1_values=E_hp_rng,
                style_props={'width': 1.0}
            )
            chart.add_xy_data(
                label='Aux.',
                x1_values=T_ext_rng,
                y1_values=E_aux_rng,
                style_props={'width': 1.0, 'bottom': E_hp_rng}
            )
            chart.add_legend()
            T_unit = f"{Q_(0, self._T_unit).units:~P.0f}"
            E_unit = f"{Q_(0, self._E_unit).units:~P.0f}"
            chart.x1.add_title(f'outdoor temperature, {T_unit}')
            chart.y1.add_title(f'energy, {E_unit}')
            return chart

    def draw_bar_chart(self) -> BarChart:
        if self.loads is not None:
            return self._draw_bar_chart()
        if self.heat_pumps is not None:
            return self._draw_bar_chart_with_heat_pump()
