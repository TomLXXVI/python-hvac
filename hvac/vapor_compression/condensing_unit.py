import time
import functools
from scipy import optimize
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.vapor_compression import VariableSpeedCompressor, FixedSpeedCompressor
from hvac.heat_exchanger.recuperator.fintube.continuous_fin import PlainFinTubeCounterFlowAirCondenser
from hvac.heat_exchanger.recuperator.fintube.continuous_fin.air_condenser import CondenserError
from hvac.logging import ModuleLogger
from .machine import Output

Q_ = Quantity

Compressor = VariableSpeedCompressor | FixedSpeedCompressor
Condenser = PlainFinTubeCounterFlowAirCondenser

logger = ModuleLogger.get_logger(__name__)


class CondensingUnit:
    """
    Models a condensing unit that can be used to analyze the steady-state
    performance of the condensing unit when the cooling capacity or heat
    absorption rate at the evaporator is given.
    """

    def __init__(
        self,
        compressor: Compressor,
        condenser: Condenser
    ) -> None:
        """
        Creates an instance of the `CondensingUnit` class.

        Parameters
        ----------
        compressor:
            Compressor model.
        condenser:
            Condenser model.
        """
        self.compressor = compressor
        self.condenser = condenser
        self.Rfg = compressor.refrigerant

        self.T_evp: Quantity | None = None
        self.n_cmp: Quantity | None = None
        self.dT_sh: Quantity | None = None
        self.air_in: HumidAir | None = None
        self.air_m_dot: Quantity | None = None
        self.evp_Q_dot: Quantity | None = None

        self.output: Output | None = None

    @staticmethod
    def time_it(fun):
        @functools.wraps(fun)
        def wrapper(self, *args, **kwargs):
            t_start = time.perf_counter()

            output = fun(self, *args, **kwargs)

            t_finish = time.perf_counter()
            dt_exec = t_finish - t_start

            logger.info(
                f"Execution time: {dt_exec} seconds"
            )

            if output.success:
                eb_err = self._check_energy_balance()
                logger.info(
                    "Error energy balance: "
                    f"absolute error = {eb_err[0].to('kW'):~P.3f}, "
                    f"relative error = {eb_err[1].to('pct'):~P.2f}"
                )
            return output

        return wrapper

    @time_it
    def rate(
        self,
        T_evp: Quantity,
        evp_Q_dot: Quantity,
        n_cmp: Quantity,
        dT_sh: Quantity,
        cnd_air_in: HumidAir,
        cnd_air_m_dot: Quantity,
        x_tol: float | None = 0.001,
        r_tol: float | None = 0.01,
        i_max: int = 20
    ) -> Output:
        """
        Determines for a given (required) cooling capacity or heat absorption
        rate `evp_Q_dot` at the evaporator at a given evaporation temperature
        `T_evp`, what the heat rejection rate at the condenser will be under the
        given operating conditions.

        The rating routine searches for the condensation temperature so that
        the heat rejection capacity of the condenser matches with the required
        heat rejection rate that follows from the given cooling capacity at the
        evaporator and the compressor power.

        Parameters
        ----------
        T_evp:
            Evaporation temperature (saturated suction gas temperature).
        evp_Q_dot:
            Required cooling capacity or heat absorption rate of the evaporator.
        n_cmp:
            Compressor speed at which the condensing unit will be analyzed.
        dT_sh:
            Degree of refrigerant superheating at the compressor inlet, which is
            determined by the setting of the expansion device.
        cnd_air_in:
            State of air entering the condenser.
        cnd_air_m_dot:
            Mass flow rate of air through the condenser.
        x_tol: optional
            Absolute tolerance for root-finding algorithm
            (`scipy.optimize.root_scalar` with 'brentq' method).
        r_tol: optional
            Relative tolerance for root-finding algorithm.
        i_max: optional
            Maximum number of iterations for root-finding algorithm

        Returns
        -------
        Object of class `Output` (see module `machine.py`).
        """
        self.T_evp = T_evp.to('degC')
        self.n_cmp = n_cmp.to('1 / min')
        self.dT_sh = dT_sh.to('K')
        self.air_in = cnd_air_in
        self.air_m_dot = cnd_air_m_dot
        self.evp_Q_dot = evp_Q_dot

        # Set the maximum value that is possible for the condensing temperature:
        T_cnd_max = self.Rfg.critical_point.T.to('degC').m - 1.0

        # Find the condensing temperature for which the heat rejection capacity
        # of the condenser balances with the required heat rejection rate
        # (i.e., the required evaporator cooling capacity, plus the mechanical
        # power added by the compressor to the refrigerant).
        counter = [0]
        try:
            res = optimize.root_scalar(
                self.__fun__,
                args=(counter,),
                method='secant',
                x0=T_cnd_max - 20.0,
                x1=T_cnd_max,
                xtol=x_tol,
                rtol=r_tol,
                maxiter=i_max
            )
        except Exception as err:
            if isinstance(err, CondenserError):
                logger.error("Operating state could not be determined.")
            else:
                logger.error(
                    f'Analysis failed due to "{type(err).__name__}: {err}".'
                )
            self.output = Output(
                evp_air_m_dot=None,
                cnd_air_m_dot=self.air_m_dot,
                evp_air_in=None,
                cnd_air_in=self.air_in,
                n_cmp=self.n_cmp,
                dT_sh=self.dT_sh
            )
            self.output.success = False
            return self.output
        else:
            logger.info(
                f"Analysis finished after {counter[0]} iterations with "
                f"message: {res.flag}"
            )
            evp_Q_dot = self.condenser.Q_dot - self.compressor.W_dot
            P_evp = self.Rfg(T=self.T_evp, x=Q_(1, 'frac')).P
            self.output = Output(
                evp_air_m_dot=None,
                cnd_air_m_dot=self.air_m_dot,
                evp_air_in=None,
                cnd_air_in=self.air_in,
                n_cmp=self.n_cmp,
                dT_sh=self.dT_sh,
                cnd_air_out=self.condenser.air_out,
                evp_Q_dot=evp_Q_dot,
                cnd_Q_dot=self.condenser.Q_dot,
                cmp_W_dot=self.compressor.W_dot,
                rfg_m_dot=self.compressor.m_dot,
                T_evp=self.T_evp,
                P_evp=P_evp,
                T_cnd=self.condenser.T_cnd,
                P_cnd=self.condenser.P_cnd,
                dT_sc=self.condenser.dT_sc,
                suction_gas=self.compressor.suction_gas,
                discharge_gas=self.condenser.rfg_in,
                liquid=self.condenser.rfg_out,
                mixture=None,
                COP=self.condenser.Q_dot / self.compressor.W_dot,
                EER=evp_Q_dot / self.compressor.W_dot,
                evp_eps=None,
                cnd_eps=self.condenser.eps,
                evp_air_dP=None,
                cnd_air_dP=self.condenser.air_dP
            )
            self.output.success = True
            return self.output

    def __fun__(self, T_cnd: float, counter: list[int]) -> float:
        """
        For the given condensation temperature `T_cnd`, determines the deviation
        between the heat rejection capacity of the condenser under the given
        operating conditions and the heat rejection rate that is required for
        turning the liquid leaving the condenser again into suction gas at the
        compressor inlet.
        """
        i = counter[0]
        T_cnd = Q_(T_cnd, 'degC')

        logger.info(
            f"Iteration {i + 1}: "
            f"Try with: T_evp = {self.T_evp:~P.3f}, T_cnd = {T_cnd:~P.3f}"
        )
        self.compressor.T_evp = self.T_evp
        self.compressor.T_cnd = T_cnd
        self.compressor.speed = self.n_cmp
        self.compressor.dT_sh = self.dT_sh

        # Determine the heat rejection capacity of the condenser:
        try:
            self.condenser.solve(
                air_m_dot=self.air_m_dot,
                rfg_m_dot=self.compressor.m_dot,
                air_in=self.air_in,
                rfg_in=self.compressor.discharge_gas
            )
        except CondenserError as err:
            logger.error(
                f"Iteration {i + 1}: "
                f"{type(err).__name__}: {err}"
            )
            raise err
        cnd_Q_dot_cap = self.condenser.Q_dot

        logger.info(
            f"Iteration {i + 1}: "
            f"Heat rejection capacity = {cnd_Q_dot_cap.to('kW'):~P.3f}."
        )

        # Determine the required heat rejection rate for the condenser:
        cnd_Q_dot_req = self.evp_Q_dot + self.compressor.W_dot

        logger.info(
            f"Iteration {i + 1}: "
            f"Required heat rejection rate = {cnd_Q_dot_req.to('kW'):~P.3f}."
        )

        dev = (cnd_Q_dot_cap - cnd_Q_dot_req).to('kW')

        logger.info(
            f"Iteration {i + 1}: "
            f"Heat rejection deviation = {dev:~P.3f}."
        )
        counter[0] += 1
        return dev.magnitude

    def _check_energy_balance(self) -> tuple[Quantity, Quantity]:
        """
        Returns the absolute error and the relative error between the heat
        rejection capacity of the condenser and the required heat rejection rate.
        """
        cnd_Q_dot_cap = self.condenser.Q_dot
        cnd_Q_dot_req = self.evp_Q_dot + self.compressor.W_dot
        abs_err = abs(cnd_Q_dot_cap - cnd_Q_dot_req)
        rel_err = abs_err / abs(cnd_Q_dot_req)
        return abs_err, rel_err.to('pct')
