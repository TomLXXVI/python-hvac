"""Analyzing a single-stage vapor compression machine at different compressor
speeds and at different outdoor air temperatures entering the condenser.

Part 2: Read the simulation results from disk and analyze (visualize) the
steady-state performance of the machine.

Continues from `07_multi_analysis_2a.py`. The file with simulation results is
saved at `analysis_results/07_multi_analysis_2.ods`.

This module has a number of functions that can be used in a Jupyter Notebook to
do the analysis.
"""
import pandas as pd
import numpy as np
from scipy import optimize
from hvac.charts import LineChart
from hvac.logging import ModuleLogger

logger = ModuleLogger.get_logger(__name__)


def read_dataframe() -> pd.DataFrame:
    """Reads ods-file into dataframe object with multi-index (n_cmp, T_outdoor).
    Rows that only contain NaN-values are dropped.
    """
    # noinspection PyTypeChecker
    df = pd.read_excel('./analysis_results/07_multi_analysis_2.ods', usecols='B:U')
    df.set_index(['n_cmp [rpm]', 'T_outdoor [degC]'], inplace=True)
    df.dropna(inplace=True, axis=0, how='all')
    return df


def _get_values_at_n_cmp(
    df: pd.DataFrame,
    n_cmp: int,
    column: str
) -> tuple[int, list, list]:
    """Gets all the values of `column` at compressor speed `n_cmp`. Returns a
    tuple with the compressor speed `n_cmp`, the range of outdoor air
    temperatures that correspond with the values in `column`, and the values
    in `column`.
    """
    selection = df.loc[(n_cmp,), column]
    T_outdoor_rng = selection.index.values
    column_values = selection.values
    return n_cmp, T_outdoor_rng, column_values


def _get_values_at_T_outdoor(
    df: pd.DataFrame,
    T_outdoor: float,
    column: str
) -> tuple[float, list, list]:
    """Gets all the values of `column` at outdoor temperature `T_outdoor`.
    Returns a tuple with the outdoor temperature `T_outdoor`, the range of
    compressor speeds that correspond with the values in `column`, and the values
    in `column`.
    """
    idx = pd.IndexSlice
    selection = df.loc[idx[:, T_outdoor], column]
    n_cmp_rng = [t[0] for t in selection.index.values]
    column_values = selection.values
    return T_outdoor, n_cmp_rng, column_values


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
        list of one or more tuples
        (`n_cmp`, `T_outdoor_rng` --> x-values, `column_values` --> y-values)
    y_unit:
        the measuring unit of the y-values
    y_scale:
        tuple with lower limit, upper limit and step size of the y-axis
    x_ax_label:
        label for the x-axis
    p_unit:
        the measuring unit of the fixed parameter
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
    df: pd.DataFrame,
    n_cmp_rng: list[int],
    column: str,
    name: str,
    unit: str,
    y_scale: tuple
) -> None:
    """Plots a line chart on screen of values in `column` as function of outdoor
    temperature at given compressor speeds in `n_cmp_rng`.

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
    chart = _create_chart(name, data, unit, y_scale, 'T_outdoor, °C', 'rpm')
    chart.show()


def _plot_T_outdoor_fixed(
    df: pd.DataFrame,
    T_out_rng: list[float],
    column: str,
    name: str,
    unit: str,
    y_scale: tuple
) -> None:
    """Plots line chart on screen of values in `column` as function of compressor
     speed at given outdoor temperatures in `T_out_rng`.

    Parameters
    ----------
    df:
        dataframe with the simulation results.
    T_out_rng:
        list of outdoor temperatures for which the `column`-values (y-values) will
        be plotted in the line chart
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
    data = [_get_values_at_T_outdoor(df, T_out, column) for T_out in T_out_rng]
    chart = _create_chart(name, data, unit, y_scale, 'n_cmp, rpm', '°C')
    chart.show()


def plot_COP(
    df: pd.DataFrame,
    n_cmp_rng: list[int] | None = None,
    T_out_rng: list[float] | None = None
) -> None:
    """If `n_cmp_rng` is not None: plots a line chart of COP against outdoor air
    temperature for each compressor speed in `n_cmp_rng`.
    If `T_out_rng` is not None: plots a line chart of COP against compressor
    speed for each outdoor temperature in `T_out_rng`.
    """
    args = 'COP [-]', 'COP', '-', (2.0, 5.5, 0.5)
    if n_cmp_rng is not None:
        _plot_n_cmp_fixed(df, n_cmp_rng, *args)
    if T_out_rng is not None:
        _plot_T_outdoor_fixed(df, T_out_rng, *args)


def plot_Q_evp(
    df: pd.DataFrame,
    n_cmp_rng: list[int] | None = None,
    T_out_rng: list[float] | None = None
) -> None:
    """If `n_cmp_rng` is not None: plots a line chart of refrigeration capacity
    against outdoor air temperature for each compressor speed in `n_cmp_rng`.
    If `T_out_rng` is not None: plots a line chart of refrigeration capacity
    against compressor speed for each outdoor temperature in `T_out_rng`.
    """
    args = 'Q_evp [kW]', 'Q_evp', 'kW', (0, 12, 1)
    if n_cmp_rng is not None:
        _plot_n_cmp_fixed(df, n_cmp_rng, *args)
    if T_out_rng is not None:
        _plot_T_outdoor_fixed(df, T_out_rng, *args)


def plot_Q_cnd(
    df: pd.DataFrame,
    n_cmp_rng: list[int] | None = None,
    T_out_rng: list[float] | None = None
) -> None:
    """If `n_cmp_rng` is not None: plots a line chart of heat rejection rate
    against outdoor air temperature for each compressor speed in `n_cmp_rng`.
    If `T_out_rng` is not None: plots a line chart of heat rejection rate
    against compressor speed for each outdoor temperature in `T_out_rng`.
    """
    args = 'Q_cnd [kW]', 'Q_cnd', 'kW', (0, 12, 1)
    if n_cmp_rng is not None:
        _plot_n_cmp_fixed(df, n_cmp_rng, *args)
    if T_out_rng is not None:
        _plot_T_outdoor_fixed(df, T_out_rng, *args)


def plot_W_cmp(
    df: pd.DataFrame,
    n_cmp_rng: list[int] | None = None,
    T_out_rng: list[float] | None = None
) -> None:
    """If `n_cmp_rng` is not None: plots a line chart of compressor power
    against outdoor air temperature for each compressor speed in `n_cmp_rng`.
    If `T_out_rng` is not None: plots a line chart of compressor power
    against compressor speed for each outdoor temperature in `T_out_rng`.
    """
    args = 'W_cmp [kW]', 'W_cmp', 'kW', (0, 12, 1)
    if n_cmp_rng is not None:
        _plot_n_cmp_fixed(df, n_cmp_rng, *args)
    if T_out_rng is not None:
        _plot_T_outdoor_fixed(df, T_out_rng, *args)


def plot_WQQ(
    df: pd.DataFrame,
    n_cmp: int | None = None,
    T_outdoor: float | None = None
) -> None:
    """Plot `W_cmp`, `Q_evp` and `Q_cnd` in a single line chart. If `n_cmp` is
    not None, the values are plotted against outdoor temperature at compressor
    speed `n_cmp`. If `T_outdoor` is not None (and `n_cmp` is None), the
    values are plotted against compressor speed at outdoor temperature `T_out`.
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
        x_ax_label = 'outdoor temperature, °C'
    elif T_outdoor is not None:
        W_cmp_data = _get_values_at_T_outdoor(df, T_outdoor, 'W_cmp [kW]')
        Q_evp_data = _get_values_at_T_outdoor(df, T_outdoor, 'Q_evp [kW]')
        Q_cnd_data = _get_values_at_T_outdoor(df, T_outdoor, 'Q_cnd [kW]')
        y_labels = [
            f'W_cmp @ {T_outdoor} °C',
            f'Q_evp @ {T_outdoor} °C',
            f'Q_cnd @ {T_outdoor} °C'
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


def plot_T_cnd_vs_T_evp(
    df: pd.DataFrame,
    n_cmp_rng: list[float]
) -> None:
    """Plots condensation temperature against evaporation temperature at
    compressor speeds in `n_cmp_rng`.
    """
    T_evp_rng_list, T_cnd_rng_list = [], []
    for n_cmp in n_cmp_rng:
        *_, T_evp_rng = _get_values_at_n_cmp(df, n_cmp, 'T_evp [degC]')
        *_, T_cnd_rng = _get_values_at_n_cmp(df, n_cmp, 'T_cnd [degC]')
        T_evp_rng_list.append(T_evp_rng)
        T_cnd_rng_list.append(T_cnd_rng)
    chart = LineChart((8, 6))
    for n_cmp, T_evp_rng, T_cnd_rng in zip(n_cmp_rng, T_evp_rng_list, T_cnd_rng_list):
        chart.add_xy_data(
            label=f'T_cnd vs. T_evp @ {n_cmp} rpm',
            x1_values=T_evp_rng,
            y1_values=T_cnd_rng,
            style_props={'marker': 'o', 'linestyle': 'none'}
        )
    chart.x1.add_title('T_evp, °C')
    chart.y1.add_title('T_cnd, °C')
    chart.x1.scale(0, 11, 1)
    chart.y1.scale(0, 80, 10)
    chart.add_legend()
    chart.show()


def get_ctrl_n_cmp(df: pd.DataFrame, target_Q_evp: float) -> tuple[list, list]:
    """Returns from `df` the range of outdoor temperatures and the corresponding
    range of compressor speeds at which the refrigeration capacity would match
    with `target_Q_evp`.
    """
    n_cmp_rng, T_out_rng = [], []
    # Group the dataframe by outdoor temperature:
    grouped_obj = df.groupby(level=1)
    # In each group (i.e. for each outdoor temperature) determine by linear
    # interpolation the compressor speed that would result in `target_Q_evp`:
    for T_out, grp_df in grouped_obj:
        try:
            # Get dataframe with the `Q_evp` values smaller than `target_Q_evp`:
            df_smaller = grp_df[grp_df['Q_evp [kW]'] < target_Q_evp]
            # Get dataframe with the `Q_evp` values greater than `target_Q_evp`:
            df_greater = grp_df[grp_df['Q_evp [kW]'] > target_Q_evp]
            # Get the smaller `Q_evp` closest to `target_Q_evp`:
            Q_evp_1 = df_smaller['Q_evp [kW]'].max()
            # Get the greater `Q_evp` closest to `target_Q_evp`:
            Q_evp_2 = df_greater['Q_evp [kW]'].min()
            # Get compressor speed that corresponds with `Q_evp_1`:
            n_cmp_1 = df_smaller['Q_evp [kW]'].idxmax()[0]
            # Get compressor speed that corresponds with `Q_evp_2`:
            n_cmp_2 = df_greater['Q_evp [kW]'].idxmin()[0]
        except ValueError:
            n_cmp_rng.append(float('nan'))
            T_out_rng.append(T_out)
        else:
            # Determine compressor speed that corresponds with `target_Q_evp`
            # using linear interpolation:
            k = (Q_evp_2 - Q_evp_1) / (n_cmp_2 - n_cmp_1)
            n_cmp = n_cmp_1 + (1 / k) * (target_Q_evp - Q_evp_1)
            # Put `n_cmp` and `T_out` in a list:
            logger.debug(
                f"at {T_out} °C: "
                f"to get refrigeration capacity {target_Q_evp} kW, "
                f"use speed {n_cmp} rpm"
            )
            n_cmp_rng.append(n_cmp)
            T_out_rng.append(T_out)
    return T_out_rng, n_cmp_rng


def curve_fit_ctrl_n_cmp(
    T_out_rng: list[float],
    n_cmp_rng: list[float]
) -> tuple[np.ndarray, np.ndarray, tuple]:
    """Returns a range of outdoor temperatures and the corresponding range of
    compressor speeds at which the refrigeration capacity match with the
    `target_Q_evp` determined with a straight-line curve-fit of the data-points
    (`T_out_rng`, `n_cmp_rng`) returned from the function `get_ctrl_n_cmp`.
    Also returns a tuple with the parameters `a` and `b`of the straight-line
    equation `a * x + b`.
    """
    # filter out any nan-values
    arr = np.array([T_out_rng, n_cmp_rng])
    arr = arr[:, ~np.any(np.isnan(arr), axis=0)]
    if arr.shape[1] > 1:
        T_out_rng, n_cmp_rng = arr[0], arr[1]

        def f(x, a, b):
            return a * x + b

        popt, *_ = optimize.curve_fit(f, T_out_rng, n_cmp_rng)

        T_out_rng = np.linspace(T_out_rng[0], T_out_rng[-1], endpoint=True)
        n_cmp_rng = np.array([f(T_out, *popt) for T_out in T_out_rng])
        return T_out_rng, n_cmp_rng, popt
    else:
        raise ValueError(
            "cannot perform curve-fitting due to invalid data points"
        )


def plot_ctrl_n_cmp(df: pd.DataFrame, target_Q_evp: float) -> tuple:
    """Plots a line chart of the compressor speed against outdoor temperature
    for which the refrigeration capacity would match the `target_Q_evp`.
    """
    try:
        T_out_rng_data, n_cmp_rng_data = get_ctrl_n_cmp(df, target_Q_evp)
        T_out_rng_crvf, n_cmp_rng_crvf, popt = curve_fit_ctrl_n_cmp(T_out_rng_data, n_cmp_rng_data)
    except ValueError as err:
        raise err
    else:
        # Create line chart:
        fig = LineChart((8, 6))
        fig.add_xy_data(
            label=f'Q_evp = {target_Q_evp} kW (data)',
            x1_values=T_out_rng_data,
            y1_values=n_cmp_rng_data,
            style_props={'marker': 'o', 'linestyle': 'none'}
        )
        fig.add_xy_data(
            label=f'Q_evp = {target_Q_evp} kW (curve-fit)',
            x1_values=T_out_rng_crvf,
            y1_values=n_cmp_rng_crvf,
        )
        fig.x1.add_title('outdoor temperature, °C')
        fig.y1.add_title('compressor speed, rpm')
        fig.add_legend(anchor='upper left', position=(0.01, 0.99))
        # Show line chart:
        fig.show()
        # Return the curve-fit parameters:
        return popt


def test_01(param='n_cmp'):
    df = read_dataframe()

    if param == 'n_cmp':
        # plot COP as function of outdoor temperature at different compressor
        # speeds
        plot_COP(df, n_cmp_rng=[3000, 3300, 3600, 3900, 4200])
        # plot refrigeration capacity as function of outdoor temperature at
        # different compressor speeds
        plot_Q_evp(df, n_cmp_rng=[3000, 3300, 3600, 3900, 4200])
        # plot compressor power as function of outdoor temperature at different
        # compressor speeds
        plot_W_cmp(df, n_cmp_rng=[3000, 3300, 3600, 3900, 4200])
        # plot compressor power, refrigeration capacity and heat rejection rate
        # as function of outdoor temperature at a given compressor speed.
        plot_WQQ(df, n_cmp=4200)

    if param == 'T_outdoor':
        # plot COP as function of compressor speed at different outdoor
        # temperatures
        plot_COP(df, T_out_rng=[26, 30, 36])
        # plot refrigeration capacity as function of compressor speed at different
        # outdoor temperatures
        plot_Q_evp(df, T_out_rng=[26, 30, 36])
        # plot compressor power as function of compressor speed at different
        # outdoor temperatures
        plot_W_cmp(df, T_out_rng=[26, 30, 36])
        # plot compressor power, refrigeration capacity and heat rejection rate
        # as function of compressor speed at a given outdoor temperature.
        plot_WQQ(df, T_outdoor=38)


def test_02(target_Q_evp: float):
    df = read_dataframe()
    plot_ctrl_n_cmp(df, target_Q_evp)


def test_03(n_cmp_rng: list[float]):
    df = read_dataframe()
    plot_T_cnd_vs_T_evp(df, n_cmp_rng)


if __name__ == '__main__':
    test_03([2700, 3000, 3300, 3600, 3900, 4200, 4500, 4800, 5100, 5400])
