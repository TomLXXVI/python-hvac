"""Analyzing a single-stage vapor compression machine at different compressor
speeds and at different mass flow rates of outdoor air through the condenser.
Other operating conditions (entering temperature of air at condenser = outdoor
temperature, entering air temperature and mass flow rate at evaporator, set
degree of superheat) remain fixed.

Part 2: Read the simulation results from disk and analyze (visualize) the
steady-state performance of the machine.

Continues from `07_multi_analysis_3a.py`. The file with simulation results is
saved at `analysis_results/07_multi_analysis_3.ods`.
"""
import pandas as pd
from hvac.charts import LineChart


def read_dataframe() -> pd.DataFrame:
    """Reads ods-file into dataframe object with multi-index (n_cmp,
    cnd_m_dot_air). Rows that only contain NaN-values are dropped.
    """
    # noinspection PyTypeChecker
    df = pd.read_excel('./analysis_results/07_multi_analysis_3.ods', usecols='B:U')
    df.set_index(['n_cmp [rpm]', 'cnd_m_dot_air [kg/h]'], inplace=True)
    df.dropna(inplace=True, axis=0, how='all')
    return df


def _get_values_at_n_cmp(
    df: pd.DataFrame,
    n_cmp: int,
    column: str
) -> tuple[int, list, list]:
    """Gets all the values of `column` at compressor speed `n_cmp`. Returns a
    tuple with the compressor speed `n_cmp`, the range of condenser air mass
    flow rates that correspond with the values in `column`, and the values
    in `column`.
    """
    selection = df.loc[(n_cmp,), column]
    cnd_m_dot_air_rng = selection.index.values
    column_values = selection.values
    return n_cmp, cnd_m_dot_air_rng, column_values


def _get_values_at_cnd_m_dot_air(
    df: pd.DataFrame,
    cnd_m_dot_air: float,
    column: str
) -> tuple[float, list, list]:
    """Gets all the values of `column` at condenser air mass flow rate
    `cnd_m_dot_air`. Returns a tuple with the condenser air mass flow rate
    `cnd_m_dot_air`, the range of compressor speeds that correspond with the
    values in `column`, and the values in `column`.
    """
    idx = pd.IndexSlice
    selection = df.loc[idx[:, cnd_m_dot_air], column]
    n_cmp_rng = [t[0] for t in selection.index.values]
    column_values = selection.values
    return cnd_m_dot_air, n_cmp_rng, column_values


def _create_chart(
    name: str,
    data: list[tuple[int | float, list, list]],
    y_unit: str,
    y_scale: tuple[float | int, ...],
    x_ax_label: str,
    p_unit: str
) -> LineChart:
    """Creates a line chart of the data contained in `data`.

    Parameters
    ----------
    name:
        label to identify the kind of y-values
    data:
        list of one or more tuples:
        (`n_cmp`, `cnd_m_dot_air_rng` -> x-values, `column_values` -> y-values)
        or (`cnd_m_dot_air`, `n_cmp_rng` -> x-values, `column_values` -> y-values)
    y_unit:
        the measuring unit of the y-values
    y_scale:
        tuple with lower limit, upper limit and step size of the y-axis
    x_ax_label:
        label for the x-axis (either condenser air mass flow rate, or compressor
        speed)
    p_unit:
        the measuring unit of the fixed parameter (either compressor speed, or
        condenser air mass flow rate)
    """
    chart = LineChart((8, 6))
    for t in data:
        chart.add_xy_data(
            label=f'{name} at {t[0]} {p_unit}',
            x1_values=t[1],
            y1_values=t[2]
        )
    chart.add_legend(columns=3)
    chart.x1.add_title(x_ax_label)
    chart.y1.add_title(f'{name}, {y_unit}')
    chart.y1.scale(y_scale[0], y_scale[1], y_scale[2])
    return chart


def _plot_n_cmp_fixed(
    df,
    n_cmp_rng: list[int],
    column: str,
    name: str,
    unit: str,
    y_scale: tuple
) -> None:
    """Plots line chart on screen of values in `column` as function of condenser
     air mass flow rate at given compressor speeds in `n_cmp_rng`.

    Parameters
    ----------
    df:
        dataframe with the simulation results.
    n_cmp_rng:
        list of compressor speeds for which the `column`-values (y-values) will
        be plotted in the line chart
    column:
        label of the column in the dataframe of which the values at the given
        range of compressor speeds will be plotted as function of outdoor air
        temperature
    name:
        label to identify the kind of column values
    unit:
        the measuring unit of the y-values
    y_scale:
        tuple with lower limit, upper limit and step size of the y-axis
    """
    data = [_get_values_at_n_cmp(df, n_cmp, column) for n_cmp in n_cmp_rng]
    chart = _create_chart(
        name, data, unit, y_scale,
        'condenser air mass flow rate, kg/h',
        'rpm'
    )
    chart.show()


def _plot_cnd_m_dot_air_fixed(
    df,
    cnd_m_dot_air_rng: list[float],
    column: str,
    name: str,
    unit: str,
    y_scale: tuple
) -> None:
    """Plots line chart on screen of values in `column` as function of compressor
     speed at given condenser air mass flow rates in `cnd_m_dot_air_rng`.

    Parameters
    ----------
    df:
        dataframe with the simulation results.
    cnd_m_dot_air_rng:
        list of condenser air mass flow rates for which the `column`-values
        (y-values) will be plotted in the line chart
    column:
        label of the column in the dataframe of which the values at the given
        range of outdoor temperatures will be plotted as function of compressor
        speed
    name:
        label to identify the kind of column values
    unit:
        the measuring unit of the y-values
    y_scale:
        tuple with lower limit, upper limit and step size of the y-axis
    """
    data = [
        _get_values_at_cnd_m_dot_air(df, cnd_m_dot_air, column)
        for cnd_m_dot_air in cnd_m_dot_air_rng
    ]
    chart = _create_chart(
        name, data, unit, y_scale,
        'compressor speed, rpm', 'kg/h'
    )
    chart.show()


def plot_COP(
    df,
    n_cmp_rng: list[int] | None = None,
    cnd_m_dot_air_rng: list[float] | None = None
) -> None:
    """If `n_cmp_rng` is not None: plots a line chart of COP against condenser
     air mass flow rate for each compressor speed in `n_cmp_rng`.
    If `cnd_m_dot_air_rng` is not None: plots a line chart of COP against
    compressor speed for each condenser air mass flow rate in `cnd_m_dot_air_rng`.
    """
    args = 'COP [-]', 'COP', '-', (2.0, 5.5, 0.5)
    if n_cmp_rng is not None:
        _plot_n_cmp_fixed(df, n_cmp_rng, *args)
    if cnd_m_dot_air_rng is not None:
        _plot_cnd_m_dot_air_fixed(df, cnd_m_dot_air_rng, *args)


def plot_Q_evp(
    df,
    n_cmp_rng: list[int] | None = None,
    cnd_m_dot_air_rng: list[float] | None = None
) -> None:
    """If `n_cmp_rng` is not None: plots a line chart of refrigeration capacity
    against condenser air mass flow rate for each compressor speed in `n_cmp_rng`.
    If `cnd_m_dot_air_rng` is not None: plots a line chart of refrigeration
    capacity against compressor speed for each condenser air mass flow rate in
    `cnd_m_dot_air_rng`.
    """
    args = 'Q_evp [kW]', 'Q_evp', 'kW', (0, 12, 1)
    if n_cmp_rng is not None:
        _plot_n_cmp_fixed(df, n_cmp_rng, *args)
    if cnd_m_dot_air_rng is not None:
        _plot_cnd_m_dot_air_fixed(df, cnd_m_dot_air_rng, *args)


def plot_Q_cnd(
    df,
    n_cmp_rng: list[int] | None = None,
    cnd_m_dot_air_rng: list[float] | None = None
) -> None:
    """If `n_cmp_rng` is not None: plots a line chart of heat rejection rate
    against condenser air mass flow rate for each compressor speed in `n_cmp_rng`.
    If `cnd_m_dot_air_rng` is not None: plots a line chart of heat rejection
    rate against compressor speed for each condenser air mass flow rate in
    `cnd_m_dot_air_rng`.
    """
    args = 'Q_cnd [kW]', 'Q_cnd', 'kW', (0, 12, 1)
    if n_cmp_rng is not None:
        _plot_n_cmp_fixed(df, n_cmp_rng, *args)
    if cnd_m_dot_air_rng is not None:
        _plot_cnd_m_dot_air_fixed(df, cnd_m_dot_air_rng, *args)


def plot_W_cmp(
    df,
    n_cmp_rng: list[int] | None = None,
    cnd_m_dot_air_rng: list[float] | None = None
) -> None:
    """If `n_cmp_rng` is not None: plots a line chart of compressor power
    against condenser air mass flow rate for each compressor speed in `n_cmp_rng`.
    If `cnd_m_dot_air_rng` is not None: plots a line chart of compressor power
    against compressor speed for each condenser air mass flow rate in
    `cnd_m_dot_air_rng`.
    """
    args = 'W_cmp [kW]', 'W_cmp', 'kW', (0, 12, 1)
    if n_cmp_rng is not None:
        _plot_n_cmp_fixed(df, n_cmp_rng, *args)
    if cnd_m_dot_air_rng is not None:
        _plot_cnd_m_dot_air_fixed(df, cnd_m_dot_air_rng, *args)


def plot_WQQ(
    df,
    n_cmp: int | None = None,
    cnd_m_dot_air: float | None = None
) -> None:
    """Plot `W_cmp`, `Q_evp` and `Q_cnd` in a single line chart. If `n_cmp` is
    not None, the values are plotted against condenser air mass flow rate at
    compressor speed `n_cmp`. If `cnd_m_dot_air` is not None (and `n_cmp` is
    None), the values are plotted against compressor speed at condenser air
    mass flow rate `cnd_m_dot_air`.
    """
    W_cmp_data: tuple = ()
    Q_evp_data: tuple = ()
    Q_cnd_data: tuple = ()
    y_labels = None
    x_ax_label = None
    if n_cmp is not None:
        W_cmp_data = _get_values_at_n_cmp(df, n_cmp, 'W_cmp [kW]')
        Q_evp_data = _get_values_at_n_cmp(df, n_cmp, 'Q_evp [kW]')
        Q_cnd_data = _get_values_at_n_cmp(df, n_cmp, 'Q_cnd [kW]')
        y_labels = [
            f'W_cmp @ {n_cmp} rpm',
            f'Q_evp @ {n_cmp} rpm',
            f'Q_cnd @ {n_cmp} rpm'
        ]
        x_ax_label = 'condenser air mass flow rate, kg/h'
    elif cnd_m_dot_air is not None:
        W_cmp_data = _get_values_at_cnd_m_dot_air(df, cnd_m_dot_air, 'W_cmp [kW]')
        Q_evp_data = _get_values_at_cnd_m_dot_air(df, cnd_m_dot_air, 'Q_evp [kW]')
        Q_cnd_data = _get_values_at_cnd_m_dot_air(df, cnd_m_dot_air, 'Q_cnd [kW]')
        y_labels = [
            f'W_cmp @ {cnd_m_dot_air} kg/h',
            f'Q_evp @ {cnd_m_dot_air} kg/h',
            f'Q_cnd @ {cnd_m_dot_air} kg/h'
        ]
        x_ax_label = 'compressor speed, rpm'
    x_values = W_cmp_data[1]
    y_values = [W_cmp_data[2], Q_evp_data[2], Q_cnd_data[2]]

    chart = LineChart((8, 6))
    for i in range(len(y_values)):
        chart.add_xy_data(
            label=y_labels[i],
            x1_values=x_values,
            y1_values=y_values[i]
        )
    chart.x1.add_title(x_ax_label)
    chart.y1.add_title('power, kW')
    chart.y1.scale(0, 12, 1)
    chart.add_legend(columns=3)
    chart.show()


def main(param='n_cmp'):
    df = read_dataframe()

    if param == 'n_cmp':
        # plot COP as function of condenser air mass flow rate at different
        # compressor speeds
        plot_COP(df, n_cmp_rng=[3000, 3300, 3600, 3900, 4200])
        # plot refrigeration capacity as function of condenser air mass flow
        # rate at different compressor speeds
        plot_Q_evp(df, n_cmp_rng=[3000, 3300, 3600, 3900, 4200])
        # plot compressor power as function of condenser air mass flow rate at
        # different compressor speeds
        plot_W_cmp(df, n_cmp_rng=[3000, 3300, 3600, 3900, 4200])
        # plot compressor power, refrigeration capacity and heat rejection rate
        # as function of condenser air mass flow rate at a given compressor speed.
        plot_WQQ(df, n_cmp=4200)

    if param == 'cnd_m_dot_air':
        # plot COP as function of compressor speed at different condenser air
        # mass flow rates
        plot_COP(df, cnd_m_dot_air_rng=[2000, 2500, 3000, 3500, 4000])
        # plot refrigeration capacity as function of compressor speed at different
        # condenser air mass flow rates
        plot_Q_evp(df, cnd_m_dot_air_rng=[2000, 2500, 3000, 3500, 4000])
        # plot compressor power as function of compressor speed at different
        # condenser air mass flow rates
        plot_W_cmp(df, cnd_m_dot_air_rng=[2000, 2500, 3000, 3500, 4000])
        # plot compressor power, refrigeration capacity and heat rejection rate
        # as function of compressor speed at a given condenser air mass flow rate.
        plot_WQQ(df, cnd_m_dot_air=3000)


if __name__ == '__main__':
    main(param='n_cmp')
