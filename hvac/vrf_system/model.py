from typing import Callable, Optional, Dict
import numpy as np
from scipy import interpolate
import dill as pickle
import pandas as pd
from hvac.logging import ModuleLogger

logger = ModuleLogger.get_logger(__name__)


class VRFModel:
    default_units = {
        'temperature': 'degC',
        'length': 'm',
        'height': 'm'
    }

    def __init__(self, units: Dict[str, str] | None = None):
        """
        Create a VRF model. This object contains all the functions that a
        `VRFSystem` instance will use to correct available capacity and input
        power depending on temperature and part-load conditions, and also on
        installation conditions (combination ratio, equivalent piping length,
        piping height).

        A separate model has to be created for cooling and for heating using the
        manufacturer's performance data.

        Parameters
        ----------
        units: Dict[str, str], optional
            Most input data to the model is dimensionless. Only temperature,
            piping length, and piping height are quantities.
            By default, it is assumed that units of temperature, length, and
            height are degrees Celsius, meter, and meter respectively. If this
            would not be the case, parameter `units` must be given a dictionary
            with keys 'temperature', 'length', and 'height' together with their
            appropriate measuring unit.
        """
        self.units = self.default_units if units is None else units
        self.CAPFT_fun: Optional[Callable] = None
        self.EIRFT_fun: Optional[Callable] = None
        self.EIRFPLR_fun: Optional[Callable] = None
        self.Leq_pipe_corr_fun: Optional[Callable] = None
        self.CR_corr_fun: Optional[Callable] = None
        self.HPRTF_fun: Optional[Callable] = None
        self.defrost_corr_fun: Optional[Callable] = None

    def set_CAPFT_fun(self, CAP_data: pd.Series | pd.DataFrame) -> None:
        """
        Set the function that returns the capacity ratio of the outdoor unit for a
        given indoor air temperature and a given outdoor air temperature (i.e. the
        ratio of available full-load capacity at the given temperatures to the
        rated full-load capacity).

        Parameters
        ----------
        CAP_data: pd.Series | pd.DataFrame
            - A Pandas Series object in case the capacity ratio should be independent
            of indoor air temperatures. In that case the index is composed of the
            outdoor air temperatures and the values are the corresponding capacity
            ratios.
            - A Pandas DataFrame object if the capacity ratio is also dependent of
            indoor air temperatures. In that case the index is composed of the indoor
            air temperatures and the columns are the outdoor air temperatures.

        Notes
        -----
        If a Pandas DataFrame object was passed in, CAPFT_fun is a function object
        with calling signature CAPFT_fun(Tia, Toa), wherein Tia is the indoor air
        temperature and Toa is the outdoor air temperature, that returns the
        capacity ratio that corresponds with the two given temperatures.

        If a Pandas Series object was passed in, CAPFT_fun is a function object
        with calling signature CAPFT_fun(Toa), wherein Toa is the outdoor air
        temperature, that returns the capacity ratio that corresponds with this
        outdoor air temperature.
        """
        if isinstance(CAP_data, pd.DataFrame):
            CAP_data = CAP_data.transpose()
            CAPFT_fun = interpolate.interp2d(
                x=CAP_data.columns,
                y=CAP_data.index,
                z=CAP_data.to_numpy()
            )
            Tia_min = min(CAP_data.columns[0], CAP_data.columns[-1])
            Tia_max = max(CAP_data.columns[0], CAP_data.columns[-1])
            Toa_min = min(CAP_data.index[0], CAP_data.index[-1])
            Toa_max = max(CAP_data.index[0], CAP_data.index[-1])
            
            def f(Tia: float, Toa: float) -> float:
                if not (Tia_min <= Tia <= Tia_max):
                    logger.warning(f"Indoor air temperature {Tia} outside interpolation domain")
                if not (Toa_min <= Toa <= Toa_max):
                    logger.warning(f"Outdoor air temperature {Toa} outside interpolation domain")
                return CAPFT_fun(Tia, Toa)[0]

            self.CAPFT_fun = f

        else:
            CAPFT_fun = interpolate.interp1d(
                x=CAP_data.index,
                y=CAP_data.values,
                bounds_error=False,
                fill_value='extrapolate'
            )
            Toa_min = min(CAP_data.index[0], CAP_data.index[-1])
            Toa_max = max(CAP_data.index[0], CAP_data.index[-1])
            
            def f(Toa: float) -> float:
                if not (Toa_min <= Toa <= Toa_max):
                    logger.warning(f"Outdoor air temperature {Toa} outside interpolation domain")
                return CAPFT_fun(Toa)

            self.CAPFT_fun = f

    def set_EIRFT_fun(self, EIR_data: pd.DataFrame) -> None:
        """
        Set the function that returns the ratio of the energy input ratio (EIR) for
        a given indoor air temperature and a given outdoor air temperature to the
        rated energy input ratio (EIR_rated) (i.e. the ratio of input power at
        rated full-load capacity to this rated full-load capacity).

        Parameters
        ----------
        EIR_data: Pandas DataFrame object
            A Pandas DataFrame with columns the outdoor air temperatures, with
            index the indoor air temperatures, and with data the EIR ratios
            derived from the manufacturer's data (EIR-ratio = input-power-ratio /
            capacity-ratio)

        Notes
        -----
        EIRFT_fun is a function object with calling signature EIRFT_fun(Tia, Toa),
        wherein Tia is the indoor air temperature and Toa is the outdoor air
        temperature, that returns the ratio EIR to EIR_rated that corresponds
        with the given temperatures.
        """
        EIR_data = EIR_data.transpose()
        EIRFT_fun = interpolate.interp2d(
            x=EIR_data.columns,
            y=EIR_data.index,
            z=EIR_data.to_numpy(),
        )
        Tia_min = min(EIR_data.columns[0], EIR_data.columns[-1])
        Tia_max = max(EIR_data.columns[0], EIR_data.columns[-1])
        Toa_min = min(EIR_data.index[0], EIR_data.index[-1])
        Toa_max = max(EIR_data.index[0], EIR_data.index[-1])
        
        def f(Tia: float, Toa: float) -> float:
            if not (Tia_min <= Tia <= Tia_max):
                logger.warning(f"Indoor air temperature {Tia} outside interpolation domain")
            if not (Toa_min <= Toa <= Toa_max):
                logger.warning(f"Outdoor air temperature {Toa} outside interpolation domain")
            return EIRFT_fun(Tia, Toa)[0]

        self.EIRFT_fun = f

    def set_EIRFPLR_fun(self, EIRFPLR_data: pd.Series) -> None:
        """
        Set the function that returns the ratio of actual input power to rated input
        power (this ratio is then defined as EIRFPLR) for a given part-load ratio
        (PLR) when PLR >= PLR_min. The EIRFPLR factor accounts for compressor speed
        changes above PLR_min (i.e. the minimum compressor part-load ratio).
        Below PLR_min the compressor will cycle on and off (see method
        `set_HPRTF_fun` for this case)

        Parameters
        ----------
        EIRFPLR_data: pd.Series
            A Pandas Series of which the index consists of PLR values and the data
            (values) of the corresponding input power ratios.

        Notes
        -----
        1. EIRFPLR_fun is a function object with calling signature EIRFPLR_fun(PLR),
        wherein PLR is the part-load ratio, that returns the ratio of actual
        input power to rated input power.

        2. The term EIRFPLR was badly chosen, as it has nothing to do with EIR
        (energy input ratio). It returns an input power ratio, and not
        a ratio of input power to capacity, which is the true definition of EIR .

        3. The EIRFPLR_data can be deduced from manufacturer's data that gives the
        ratio of power input in function of total capacity of indoor units. The
        capacity of the indoor units will finally balance with the load on the
        indoor units, and a such they can be considered the same. Consequently, the
        PLR can also be considered as a ratio of the actual total capacity of indoor
        units to the rated total capacity of indoor units.
        """
        PLR_min = min(EIRFPLR_data.index[0], EIRFPLR_data.index[-1])
        PLR_max = max(EIRFPLR_data.index[0], EIRFPLR_data.index[-1])

        EIRFPLR_fun = interpolate.interp1d(
            x=EIRFPLR_data.index,
            y=EIRFPLR_data.values,
            bounds_error=False,
            fill_value='extrapolate'
        )

        def f(PLR: float) -> float:
            if not (PLR_min <= PLR <= PLR_max):
                logger.warning(f"PLR-value {PLR} outside interpolation domain")
            if PLR >= PLR_min:
                return EIRFPLR_fun(PLR)
            else:
                return EIRFPLR_fun(PLR_min)

        self.EIRFPLR_fun = f

    def set_CR_corr_fun(self, CR_corr_data: pd.Series) -> None:
        """
        Set the combination ratio correction function to correct for outdoor unit
        capacity.

        Parameters
        ----------
        CR_corr_data: pd.Series
            A Pandas Series object of which the index consists of combination ratio
            (CR) values and the data are the corresponding capacity ratio's.

        Notes
        -----
        1. CR_corr_fun is a function object with calling signature CR_corr_fun(CR),
        wherein CR is the fixed combination ratio (or model size ratio) of the
        system, that returns the capacity correction factor for the given
        combination ratio.

        2. The CR_corr_data can be deduced from manufacturer's data that gives the
        ratio of capacity of the outdoor unit in function of total capacity of
        indoor units.
        """
        CR_min = min(CR_corr_data.index[0], CR_corr_data.index[-1])
        CR_max = max(CR_corr_data.index[0], CR_corr_data.index[-1])

        CR_corr_fun = interpolate.interp1d(
            x=CR_corr_data.index,
            y=CR_corr_data.values,
            bounds_error=False,
            fill_value=(CR_min, CR_max)
        )

        def f(CR: float) -> float:
            if not (CR_min <= CR <= CR_max):
                logger.warning(f"Combination ratio {CR} outside interpolation domain")
            return CR_corr_fun(CR)

        self.CR_corr_fun = f

    def set_defrost_corr_fun(self, defrost_data: pd.Series):
        """
        Set the defrost correction function to correct for outdoor unit capacity
        (only used in heating mode).

        Parameters
        ----------
        defrost_data:
            A Pandas Series object of which the index consists of outdoor air
            temperatures and the values are corresponding correction factors.
        """
        defrost_corr_fun = interpolate.interp1d(
            x=defrost_data.index,
            y=defrost_data.values,
            bounds_error=False,
            fill_value=(defrost_data.values[-1], 1.0)
        )

        def f(defrost_cf: float) -> float:
            return defrost_corr_fun(defrost_cf)

        self.defrost_corr_fun = f

    def set_Leq_pipe_corr_fun(
        self,
        L_eq_corr_data: pd.DataFrame | pd.Series,
        cf_height: float = 0.0
    ) -> None:
        """
        Set the equivalent piping length and height correction function for
        the correction of available capacity.

        Parameters
        ----------
        L_eq_corr_data: pd.DataFrame | pd.Series
            - A Pandas DataFrame object having a list of equivalent piping lengths as
            index, the columns are combination ratios and the data are the corresponding
            correction factors.
            - A Pandas Series objects having a list of equivalent piping lengths as
            index, the data are the corresponding correction factors.
        cf_height: float, default 0.0
            Correction factor for the height difference between outdoor unit and
            highest or lowest indoor unit.

        Notes
        -----
        If a Pandas DataFrame object is passed in (cooling mode) and cf_height is
        not zero, Leq_pipe_corr_fun is a function object with calling signature
        Leq_pipe_corr_fun(Leq_pipe, CR, h), wherein Leq_pipe is the equivalent
        pipe length between the outdoor unit and the farthest indoor unit, CR is
        the combination ratio (model size ratio) of the VRF system, and
        h is the height between the outdoor unit and the highest or lowest indoor unit
        or the height between the lowest and highest indoor unit in case indoor units
        are above and below the outdoor unit. The function returns the correction factor
        that must be applied to the available capacity of the outdoor unit.
        If cf_height is zero, no correction for piping height is incorporated and
        the calling signature is f(Leq_pipe, CR).

        If a Pandas Series is passed in (heating mode) and cf_height is
        not zero, the calling signature becomes Leq_pipe_corr_fun(Leq_pipe, h).
        If cf_height is zero, the calling signature will be Leq_pipe_corr_fun(Leq_pipe).
        """
        if isinstance(L_eq_corr_data, pd.DataFrame):
            L_eq_corr_data = L_eq_corr_data.transpose()
            Leq_corr_fun = interpolate.interp2d(
                x=L_eq_corr_data.columns,
                y=L_eq_corr_data.index,
                z=L_eq_corr_data.to_numpy()
            )
            Leq_min = min(L_eq_corr_data.columns[0], L_eq_corr_data.columns[-1])
            Leq_max = max(L_eq_corr_data.columns[0], L_eq_corr_data.columns[-1])
            CR_min = min(L_eq_corr_data.index[0], L_eq_corr_data.index[-1])
            CR_max = max(L_eq_corr_data.index[0], L_eq_corr_data.index[-1])

            if cf_height != 0.0:

                def f(Leq: float | int, CR: float, h: float | int) -> float:
                    if not (Leq_min <= Leq <= Leq_max):
                        logger.warning(f"Equivalent pipe length {Leq} outside interpolation domain")
                    if not (CR_min <= CR <= CR_max):
                        logger.warning(f"Combination ratio {CR} outside interpolation domain")
                    return Leq_corr_fun(Leq, CR)[0] + cf_height * h

                self.Leq_pipe_corr_fun = f
            else:

                def f(Leq: float | int, CR: float) -> float:
                    if not (Leq_min <= Leq <= Leq_max):
                        logger.warning(f"Equivalent pipe length {Leq} outside interpolation domain")
                    if not (CR_min <= CR <= CR_max):
                        logger.warning(f"Combination ratio {CR} outside interpolation domain")
                    return Leq_corr_fun(Leq, CR)[0]

                self.Leq_pipe_corr_fun = f
        else:
            Leq_corr_fun = interpolate.interp1d(
                x=L_eq_corr_data.index,
                y=L_eq_corr_data.values,
                bounds_error=False,
                fill_value='extrapolate'
            )
            Leq_min = min(L_eq_corr_data.index[0], L_eq_corr_data.index[-1])
            Leq_max = max(L_eq_corr_data.index[0], L_eq_corr_data.index[-1])

            if cf_height != 0.0:

                def f(Leq: float | int, h: float | int) -> float:
                    if not (Leq_min <= Leq <= Leq_max):
                        logger.warning(f"Equivalent pipe length {Leq} outside interpolation domain")
                    return Leq_corr_fun(Leq)[0] + cf_height * h

                self.Leq_pipe_corr_fun = f
            else:

                def f(Leq: float | int) -> float:
                    if not (Leq_min <= Leq <= Leq_max):
                        logger.warning(f"Equivalent pipe length {Leq} outside interpolation domain")
                    return Leq_corr_fun(Leq)

                self.Leq_pipe_corr_fun = f

    def set_HPRTF_fun(self, PLR_min: float = 0.5) -> None:
        """
        Set the function that calculates the heat pump runtime fraction
        when PLR is smaller than PLR_min.

        Parameters
        ----------
        PLR_min: float, default 0.5
            Minimum part-load ratio at which capacity modulation by controlling
            compressor speed is possible. When PLR < PLR_min, the compressor
            will cycle.

        Notes
        -----
        HPRTF_fun is a function object with calling signature HPRTF_fun(PLR),
        wherein PLR is a part-load ratio between 0 and 1, that returns the
        heat pump runtime fraction.
        """
        def HPRTF_fun(PLR: float):
            if PLR < PLR_min:
                CR = PLR / PLR_min  # cycling ratio
                CR_frac = 0.85 + 0.15 * CR  # cycling ratio fraction
                return CR / CR_frac  # heat pump runtime fraction
            else:
                return 1.0  # neutralize HPRTF when PLR >= PLR_min

        self.HPRTF_fun = HPRTF_fun

    def save(self, file_path: str):
        """Save (pickle) the VRF-model to disk."""
        with open(file_path, 'wb') as fh:
            pickle.dump(self, fh)

    @staticmethod
    def load(file_path: str) -> 'VRFModel':
        """Load the VRF-model from disk."""
        with open(file_path, 'rb') as fh:
            vrf_model = pickle.load(fh)
            return vrf_model


class VRFModelCreator:
    """Helper class for the user to create the `VRFModel`."""
    class QTCurve:
        low_limit = -20
        high_limit = 15

        @staticmethod
        def _get_interpolant(points):
            x, y = zip(*points)
            interpolant = interpolate.interp1d(x, y)
            return interpolant

        def __init__(self, points: list[tuple[float, float]]):
            self._interpolant = self._get_interpolant(points)

        def __call__(self, Tao: float) -> float:
            if self.low_limit <= Tao <= self.high_limit:
                return self._interpolant(Tao)
            else:
                raise OverflowError(
                    f'outside air temperature not within limits '
                    f'[{self.low_limit}, {self.high_limit}]'
                )

    class WTCurve(QTCurve):
        pass

    def __init__(
        self,
        Tia_list: list[float | int],
        Tao_limits: tuple[float | int, float | int] = (-20, 15),
        units: Dict[str, str] | None = None
    ):
        """
        Create `VRFModelCreator` instance.

        Parameters
        ----------
        Tia_list:
            List of the indoor air temperatures for which curves of input power
            ratios as function of outdoor air temperature are given.
        Tao_limits:
            The lower and upper limit of the outdoor air temperature axis of
            the curves that give capacity ratio and input power ratio as a
            function of outdoor air temperature.
        units:
            Most input data to the model is dimensionless. Only temperature,
            piping length, and piping height are quantities.
            By default, it is assumed that units of temperature, length, and
            height are degrees Celsius, meter, and meter respectively. If this
            would not be the case, parameter `units` must be given a dictionary
            with keys 'temperature', 'length', and 'height' together with their
            appropriate measuring unit used in the manufacturer's datasheets.
        """
        self.QTCurve.low_limit = self.WTCurve.low_limit = Tao_limits[0]
        self.QTCurve.high_limit = self.WTCurve.high_limit = Tao_limits[1]
        self.QT_curve = None
        self.WT_curves = {float(Tia): None for Tia in Tia_list}
        self.vrf_model = VRFModel(units)
        self.Toa_arr = np.arange(self.QTCurve.low_limit, self.QTCurve.high_limit + 1.0, 1.0)

    def create_QT_curve(self, points: list[tuple[float, float]]):
        self.QT_curve = self.QTCurve(points)

    def create_WT_curve(self, Tia: float | int, points: list[tuple[float, float]]):
        self.WT_curves[float(Tia)] = self.WTCurve(points)

    def create_CAPFT_function(self):
        Q_ratio_list = [self.QT_curve(Toa) for Toa in self.Toa_arr]
        CAP_data = pd.Series(data=Q_ratio_list, index=self.Toa_arr)
        self.vrf_model.set_CAPFT_fun(CAP_data)

    def create_EIRFT_function(self):

        def WT_curve(Tia: float | int, Toa: float | int):
            return self.WT_curves[float(Tia)](Toa)

        EIR_ratio_list = [
            [WT_curve(Tia, Toa) / self.QT_curve(Toa) for Toa in self.Toa_arr]
            for Tia in self.WT_curves.keys()
        ]
        EIR_data = pd.DataFrame(
            data=EIR_ratio_list,
            index=list(self.WT_curves.keys()),
            columns=self.Toa_arr
        )
        self.vrf_model.set_EIRFT_fun(EIR_data)

    def create_defrost_corr_function(self, points: list[tuple[float, float]]):
        Tao_list, defrost_coeff_list = zip(*points)
        defrost_corr_data = pd.Series(index=Tao_list, data=defrost_coeff_list)
        self.vrf_model.set_defrost_corr_fun(defrost_corr_data)

    def create_EIRFPLR_function(self, Qiu_rated: int, points: list[tuple[float, float]]):
        Qiu_list, W_ratio_list = zip(*points)
        PLR_values = [Qiu / Qiu_rated for Qiu in Qiu_list]
        EIRFPLR_data = pd.Series(index=PLR_values, data=W_ratio_list)
        self.vrf_model.set_EIRFPLR_fun(EIRFPLR_data)

    def create_CR_corr_function(self, Qiu_rated: int, points: list[tuple[float, float]]):
        Qiu_list, Q_ratio_list = zip(*points)
        CR_values = [Qiu / Qiu_rated for Qiu in Qiu_list]
        CR_corr_data = pd.Series(index=CR_values, data=Q_ratio_list)
        self.vrf_model.set_CR_corr_fun(CR_corr_data)

    def create_Leq_corr_function(self, points: list[tuple[float, float]]):
        L_eq_list, cf_list = zip(*points)
        Leq_pipe_corr_data = pd.Series(index=L_eq_list, data=cf_list)
        self.vrf_model.set_Leq_pipe_corr_fun(Leq_pipe_corr_data)

    def create_HPRTF_function(self, PLR_min: float = 0.5):
        self.vrf_model.set_HPRTF_fun(PLR_min)

    def save_model(self, fp: str):
        self.vrf_model.save(fp)
