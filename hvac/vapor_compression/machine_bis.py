import warnings
import logging
import numpy as np
from scipy import optimize
from hvac import Quantity
from hvac.logging import Logger
from hvac.fluids import HumidAir, Fluid, FluidState, CP_HUMID_AIR
from hvac.heat_transfer.heat_exchanger.fin_tube import air_evaporator
from hvac.heat_transfer.heat_exchanger.fin_tube import air_condenser
from hvac.vapor_compression.real_compressor import FixedSpeedCompressor, VariableSpeedCompressor

Q_ = Quantity

Evaporator = air_evaporator.rating.PFT_CO_EVP
Condenser = air_condenser.rating.PFT_CO_CND
Compressor = FixedSpeedCompressor | VariableSpeedCompressor


CondenserWarning = air_condenser.rating.plain_fin_tube_condenser.CondenserWarning


class SingleStageVaporCompressionMachine:
    """Model of a single-stage vapor compression machine."""

    def __init__(
        self,
        evaporator: Evaporator,
        condenser: Condenser,
        compressor: Compressor,
        Refrigerant: Fluid,
        n_cmp_min: Quantity | None = None,
        n_cmp_max: Quantity | None = None,
        logger: Logger | None = None
    ) -> None:
        """Creates a `SingleStageVaporCompressionMachine` instance.

        Parameters
        ----------
        evaporator:
            Model instance of evaporator.
        condenser:
            Model instance of condenser.
        compressor:
            Model instance of compressor.
        Refrigerant:
            Refrigerant used in vapor compression machine.
        n_cmp_min: optional
            Minimum compressor speed in case of variable speed compressor.
        n_cmp_max: optional
            Maximum compressor speed in case of variable speed compressor.
        logger:
            Logger for sending messages about the course of the rating routine.
        """
        self.evaporator = evaporator
        self.condenser = condenser
        self.compressor = compressor
        self.Refrigerant = Refrigerant
        self.n_cmp_min = n_cmp_min
        self.n_cmp_max = n_cmp_max
        self._logger = logger

        self.evp_m_dot_air: Quantity | None = None
        self.cnd_m_dot_air: Quantity | None = None
        self.evp_air_in: HumidAir | None = None
        self.cnd_air_in: HumidAir | None = None
        self.dT_sh: Quantity | None = None
        self.n_cmp: Quantity | None = None
        self.evp_rfg_sat_vap: FluidState | None = None
        self.cnd_rfg_sat_liq: FluidState | None = None

    def _log(self, msg: str, level: int) -> None:
        if self._logger is not None:
            match level:
                case logging.DEBUG:
                    self._logger.debug(msg)
                case logging.INFO:
                    self._logger.info(msg)
                case logging.WARNING:
                    self._logger.warning(msg)
                case logging.ERROR:
                    self._logger.error(msg)

    def _warn_handler(self, message, category, filename, lineno, file=None, line=None):
        self._log(message, logging.WARNING)

    def set_operating_conditions(
        self,
        evp_m_dot_air: Quantity,
        cnd_m_dot_air: Quantity,
        evp_air_in: HumidAir,
        cnd_air_in: HumidAir,
        dT_sh: Quantity,
        n_cmp: Quantity | None = None
    ) -> None:
        """Sets the known operating conditions on the vapor compression machine
        that will determine the performance of this machine.

        Parameters
        ----------
        evp_m_dot_air:
            Mass flow rate of air through evaporator.
        cnd_m_dot_air:
            Mass flow rate of air through condenser.
        evp_air_in:
            State of air entering evaporator.
        cnd_air_in:
            State of air entering condenser.
        dT_sh:
            Degree of superheat of refrigerant leaving evaporator, set on
            expansion device of vapor compression machine.
        n_cmp: optional
            Current compressor speed (only for variable speed compressors).

        Returns
        -------
        None
        """
        self.evp_m_dot_air = evp_m_dot_air
        self.cnd_m_dot_air = cnd_m_dot_air
        self.evp_air_in = evp_air_in
        self.cnd_air_in = cnd_air_in
        self.dT_sh = dT_sh.to('K')
        self.n_cmp = n_cmp

    def rate(
        self,
        T_evp_ini: Quantity,
        T_cnd_ini: Quantity,
        i_max: int = 50,
        abs_tol: Quantity = Q_(0.1, 'kg / hr'),
        rel_tol: Quantity = Q_(1, 'pct')
    ) -> None:
        """Determines the steady-state performance of the single-stage vapor
        compression machine for the operating conditions set with
        `set_operating_conditions`.

        Parameters
        ----------
        T_evp_ini:
            Initial guess of the evaporation temperature.
        T_cnd_ini:
            Initial guess of the condensation temperature.
        i_max:
            Maximum number of iterations.
        abs_tol:
            Absolute tolerance on the difference between the refrigerant mass
            flow rate according to the compressor model and according to the
            evaporator model.
        rel_tol:
            Relative tolerance (ratio of the absolute tolerance on the
            refrigerant mass flow rate according to the compressor model).

        Returns
        -------
        None
        """
        if isinstance(self.compressor, VariableSpeedCompressor):
            self.compressor.speed = self.n_cmp
        self.compressor.dT_sh = self.dT_sh
        T_evp = T_evp_ini.to('degC')
        T_cnd = T_cnd_ini.to('degC')
        for i in range(i_max):
            self._log((
                f"Iteration {i + 1}: "
                f"try with T_evp = {T_evp:~P.3f} and "
                f"T_cnd = {T_cnd:~P.3f}"
                ), logging.INFO
            )
            # Check if condenser temperature is below critical temperature of
            # refrigerant:
            T_crit = self.Refrigerant.critical_point.T
            if T_cnd >= T_crit:
                err_message = (
                    "Rating of VCM failed: "
                    f"Condenser temperature ({T_cnd:~P.3f}) exceeds critical "
                    f"temperature of refrigerant ({T_crit.to('degC'):~P.3f})."
                )
                self._log(err_message, logging.ERROR)
                raise ValueError(err_message)
            # Check if condenser temperature is above temperature of air
            # entering the condenser:
            if T_cnd <= self.cnd_air_in.Tdb:
                err_message = (
                    "Rating of VCM failed: "
                    f"Condenser temperature ({T_cnd:~P.3f}) is lower than "
                    f"temperature of air entering condenser "
                    f"({self.cnd_air_in.Tdb.to('degC'):~P.3f})."
                )
                self._log(err_message, logging.ERROR)
                raise ValueError(err_message)
            # Check if evaporation temperature is lower than temperature of
            # air entering evaporator:
            if T_evp >= self.evp_air_in.Tdb:
                err_message = (
                    "Rating of VCM failed: "
                    f"Evaporator temperature ({T_evp:~P.3f}) is higher than "
                    f"temperature of air entering evaporator "
                    f"({self.evp_air_in.Tdb.to('degC'):~P.3f})."
                )
                self._log(err_message, logging.ERROR)
                raise ValueError(err_message)
            # Set compressor parameters and get mass flow rate of compressor:
            self.compressor.Te = T_evp
            self.compressor.Tc = T_cnd
            cmp_m_dot_rfg = self.compressor.m_dot.to('kg / hr')
            self._log((
                f"Iteration {i + 1}: "
                f"compressor mass flow rate: {cmp_m_dot_rfg:~P.3f}, "
                f"compressor power: {self.compressor.Wc_dot.to('kW'):~P.3f}"
                ), logging.DEBUG
            )
            cnd_rfg_in = self.compressor.discharge_gas
            self._log((
                f"Iteration {i + 1}: "
                "refrigerant at condenser entry: "
                f"{cnd_rfg_in.T.to('degC'):~P.3f}, "
                f"{cnd_rfg_in.P.to('bar'):~P.3f}"
                ), logging.DEBUG
            )
            # Determine performance of condenser with mass flow rate of
            # compressor:
            with warnings.catch_warnings(category=CondenserWarning):
                warnings.showwarning = self._warn_handler
                # Note: we use this to send `CondenserWarnings` to the logger
                r_cnd = self.condenser(
                    m_dot_air=self.cnd_m_dot_air,
                    m_dot_rfg=cmp_m_dot_rfg,
                    air_in=self.cnd_air_in,
                    rfg_in=cnd_rfg_in,
                )
            self._log((
                f"Iteration {i + 1}: "
                "refrigerant at condenser exit: "
                f"{r_cnd.rfg_out.T.to('degC'):~P.3f}, "
                f"{r_cnd.rfg_out.P.to('bar'):~P.3f}, "
                f"{r_cnd.rfg_out.h.to('kJ / kg'):~P.3f}, "
                f"{r_cnd.rfg_out.x.to('frac'):~P.3f}, "
                f"{r_cnd.rfg_out.phase}"
                ), logging.DEBUG
            )
            self._log((
                f"Iteration {i + 1}: "
                "heat rejection rate: "
                f"{r_cnd.Q_dot.to('kW'):~P.3f}"
                ), logging.DEBUG
            )
            # Determine state of refrigerant entering evaporator with mass
            # flow rate of compressor:
            evp_rfg_in = self.Refrigerant(
                h=r_cnd.rfg_out.h,
                P=self.compressor.Pe
            )
            self._log((
                f"Iteration {i + 1}: "
                f"refrigerant at evaporator entry: "
                f"{evp_rfg_in.T.to('degC'):~P.3f}, "
                f"{evp_rfg_in.P.to('bar'):~P.3f}, "
                f"{evp_rfg_in.h.to('kJ / kg'):~P.3f}, "
                f"{evp_rfg_in.x.to('frac'):~P.3f}, "
                f"{evp_rfg_in.phase}"
                ), logging.DEBUG
            )
            # Determine the mass flow rate let through by the expansion device
            # at the given operating conditions of the evaporator in order
            # to maintain the set degree of superheat:
            r_evp = self.evaporator(
                air_in=self.evp_air_in,
                m_dot_air=self.evp_m_dot_air,
                rfg_in=evp_rfg_in,
                dT_rfg_sh=self.dT_sh,
                m_dot_rfg_ini=cmp_m_dot_rfg
            )
            self._log((
                f"Iteration {i + 1}: "
                f"refrigerant at evaporator exit: "
                f"{r_evp.rfg_out.T.to('degC'):~P.3f}, "
                f"{r_evp.rfg_out.P.to('bar'):~P.3f}"
                ), logging.DEBUG
            )
            self._log((
                f"Iteration {i + 1}: "
                f"heat absorption rate: "
                f"{r_evp.Q_dot.to('kW'):~P.3f}"
                ), logging.DEBUG
            )
            evp_m_dot_rfg = r_evp.m_dot_rfg.to('kg / hr')
            self._log((
                f"Iteration {i + 1}: "
                f"evaporator mass flow rate: "
                f"{evp_m_dot_rfg.to('kg / hr'):~P.3f}"
                ), logging.DEBUG
            )
            # Determine new value for evaporation temperature based on the
            # heat rejection rate at the condenser and compressor power.
            Q_evp = r_cnd.Q_dot - self.compressor.Wc_dot
            # `Q_evp` is to be considered as the refrigeration capacity of the
            # condensing unit (= compressor + condenser); from the heat transfer
            # equation of the evaporator (based on a wet surface) the
            # evaporation temperature can be deduced using the heat transfer
            # effectiveness of the evaporator:
            h_air_out_sat = (
                self.evaporator.air_in.h
                - Q_evp / (r_evp.eps * self.evaporator.m_dot_air)
            )
            air_out_sat = HumidAir(h=h_air_out_sat, RH=Q_(100, 'pct'))
            T_evp_new = air_out_sat.Tdb
            # Determine new value for condensation temperature based on the
            # heat absorption rate at the evaporator (determined with the
            # mass flow rate let through by the expansion device) and compressor
            # power:
            Q_cnd = self.evaporator.Q_dot + self.compressor.Wc_dot
            # From the heat transfer equation of the condenser the condensation
            # temperature can be deduced using the heat transfer effectiveness
            # of the condenser:
            T_cnd_new = (
                self.condenser.air_in.Tdb +
                Q_cnd / (r_cnd.eps * CP_HUMID_AIR * self.condenser.m_dot_air)
            )
            # Note: the law of conservation of energy requires that at
            # steady-state the sum of heat absorption rate at the evaporator and
            # compressor power equals the heat rejection rate at the condenser.
            # Determine deviation between mass flow rate let through by
            # expansion device (to maintain set degree of superheat) and mass
            # flow rate displaced by the compressor; the law of continuity
            # requires that both mass flow rates are equal at steady-state
            # operation:
            abs_dev = abs(evp_m_dot_rfg - cmp_m_dot_rfg)
            rel_dev = abs_dev / abs(cmp_m_dot_rfg)
            self._log((
                f"Iteration {i + 1}: deviation "
                f"with {T_evp:~P.3f} and {T_cnd:~P.3f}: "
                f"{abs_dev.to('kg / hr'):~P.3f}, {rel_dev.to('pct'):~P3f}"
                ), logging.INFO
            )
            if (rel_dev <= rel_tol) or (abs_dev <= abs_tol):
                self._log((
                    f"Rating finished after {i + 1} iterations."
                    ), logging.INFO
                )
                self.compressor.Te = T_evp_new
                self.compressor.Tc = T_cnd_new
                break
                # Note: when mass flow rates are (approximately) balanced, the
                # law of conservation of energy should also (approximately) be
                # fulfilled (by making `tol` smaller, the approximations should
                # also become smaller; however, the number of iterations and
                # so calculation time will also increase)
            T_evp = T_evp_new.to('degC')
            T_cnd = T_cnd_new.to('degC')
        else:
            self._log((
                'Steady-state performance could not be determined.'
                ), logging.ERROR
            )
            raise ValueError('Steady-state performance could not be determined.')

    def _eq_root(self, unknowns: np.ndarray, counter: list) -> np.ndarray:
        i = counter[0]
        T_evp = Q_(unknowns[0], 'degC')
        T_cnd = Q_(unknowns[1], 'degC')
        self._log((
            f"Iteration {i + 1}: "
            "try with "
            f"T_evp = {T_evp:~P.3f} and "
            f"T_cnd = {T_cnd:~P.3f}"
            ), logging.INFO
        )
        # Check if condenser temperature is below critical temperature of
        # refrigerant
        T_crit = self.Refrigerant.critical_point.T
        if T_cnd >= T_crit:
            err_message = (
                "Rating of VCM failed: "
                f"Condenser temperature ({T_cnd:~P.3f}) exceeds critical "
                f"temperature of refrigerant ({T_crit.to('degC'):~P.3f})."
            )
            self._log(err_message, logging.ERROR)
            raise ValueError(err_message)
        # Check if condenser temperature is above temperature of air
        # entering the condenser:
        if T_cnd <= self.cnd_air_in.Tdb:
            err_message = (
                "Rating of VCM failed: "
                f"Condenser temperature ({T_cnd:~P.3f}) is lower than "
                f"temperature of air entering condenser "
                f"({self.cnd_air_in.Tdb.to('degC'):~P.3f})."
            )
            self._log(err_message, logging.ERROR)
            raise ValueError(err_message)
        # Check if evaporation temperature is lower than temperature of
        # air entering evaporator:
        if T_evp >= self.evp_air_in.Tdb:
            err_message = (
                "Rating of VCM failed: "
                f"Evaporator temperature ({T_evp:~P.3f}) is higher than "
                f"temperature of air entering evaporator "
                f"({self.evp_air_in.Tdb.to('degC'):~P.3f})."
            )
            self._log(err_message, logging.ERROR)
            raise ValueError(err_message)
        self.compressor.Te = T_evp
        self.compressor.Tc = T_cnd
        cmp_m_dot_rfg = self.compressor.m_dot.to('kg / hr')
        self._log((
            f"Iteration {i + 1}: "
            f"compressor mass flow rate: {cmp_m_dot_rfg:~P.3f}, "
            f"compressor power: {self.compressor.Wc_dot.to('kW'):~P.3f}"
            ), logging.DEBUG
        )
        cnd_rfg_in = self.compressor.discharge_gas
        self._log((
            f"Iteration {i + 1}: "
            "refrigerant at condenser entry: "
            f"{cnd_rfg_in.T.to('degC'):~P.3f}, "
            f"{cnd_rfg_in.P.to('bar'):~P.3f}"
            ), logging.DEBUG
        )
        with warnings.catch_warnings(category=CondenserWarning):
            warnings.showwarning = self._warn_handler
            r_cnd = self.condenser(
                m_dot_air=self.cnd_m_dot_air,
                m_dot_rfg=cmp_m_dot_rfg,
                air_in=self.cnd_air_in,
                rfg_in=cnd_rfg_in,
            )
        self._log((
            f"Iteration {i + 1}: "
            "refrigerant at condenser exit: "
            f"{r_cnd.rfg_out.T.to('degC'):~P.3f}, "
            f"{r_cnd.rfg_out.P.to('bar'):~P.3f}, "
            f"{r_cnd.rfg_out.h.to('kJ / kg'):~P.3f}, "
            f"{r_cnd.rfg_out.x.to('frac'):~P.3f}, "
            f"{r_cnd.rfg_out.phase}"
            ), logging.DEBUG
        )
        self._log((
            f"Iteration {i + 1}: "
            "heat rejection rate: "
            f"{r_cnd.Q_dot.to('kW'):~P.3f}"
            ), logging.DEBUG
        )
        evp_rfg_in = self.Refrigerant(
            h=r_cnd.rfg_out.h,
            P=self.compressor.Pe
        )
        self._log((
            f"Iteration {i + 1}: "
            f"refrigerant at evaporator entry: "
            f"{evp_rfg_in.T.to('degC'):~P.3f}, "
            f"{evp_rfg_in.P.to('bar'):~P.3f}, "
            f"{evp_rfg_in.h.to('kJ / kg'):~P.3f}, "
            f"{evp_rfg_in.x.to('frac'):~P.3f}, "
            f"{evp_rfg_in.phase}"
            ), logging.DEBUG
        )
        r_evp = self.evaporator(
            air_in=self.evp_air_in,
            m_dot_air=self.evp_m_dot_air,
            rfg_in=evp_rfg_in,
            dT_rfg_sh=self.dT_sh,
            m_dot_rfg_ini=cmp_m_dot_rfg
        )
        self._log((
            f"Iteration {i + 1}: "
            f"refrigerant at evaporator exit: "
            f"{r_evp.rfg_out.T.to('degC'):~P.3f}, "
            f"{r_evp.rfg_out.P.to('bar'):~P.3f}"
            ), logging.DEBUG
        )
        self._log((
            f"Iteration {i + 1}: "
            f"heat absorption rate: "
            f"{r_evp.Q_dot.to('kW'):~P.3f}"
            ), logging.DEBUG
        )
        self._log((
            f"Iteration {i + 1}: "
            f"evaporator mass flow rate: "
            f"{r_evp.m_dot_rfg.to('kg / hr'):~P.3f}"
            ), logging.DEBUG
        )
        Q_evp = r_cnd.Q_dot - self.compressor.Wc_dot
        h_air_out_sat = (
                self.evaporator.air_in.h
                - Q_evp / (r_evp.eps * self.evaporator.m_dot_air)
        )
        air_out_sat = HumidAir(h=h_air_out_sat, RH=Q_(100, 'pct'))
        T_evp_new = air_out_sat.Tdb
        Q_cnd = self.evaporator.Q_dot + self.compressor.Wc_dot
        T_cnd_new = (
            self.condenser.air_in.Tdb +
            Q_cnd / (r_cnd.eps * CP_HUMID_AIR * self.condenser.m_dot_air)
        )
        dev_T_evp = T_evp_new.to('K') - T_evp.to('K')
        dev_T_cnd = T_cnd_new.to('K') - T_cnd.to('K')
        self._log((
            f"Iteration {i + 1}: deviation "
            f"with {T_evp:~P.3f} and {T_cnd:~P.3f}: "
            f"for T_evp = {dev_T_evp:~P.3f}, T_cnd = {dev_T_cnd:~P.3f}"
        ), logging.INFO
        )
        counter[0] += 1
        return np.array([dev_T_evp.m, dev_T_cnd.m])

    def rate_root(
        self,
        T_evp_ini: Quantity,
        T_cnd_ini: Quantity,
        i_max: int = 50,
        rel_tol: Quantity = Q_(1.0, 'pct')
    ) -> None:
        """Alternative implementation of `rate` using `scipy.optimize.root`.

        Parameters
        ----------
        T_evp_ini:
            Initial guess of the evaporation temperature.
        T_cnd_ini:
            Initial guess of the condensation temperature.
        i_max:
            Maximum number of iterations.
        rel_tol:
            The algorithm will stop if the relative deviation between 2
            consecutive values of `T_evp` and `T_cnd` is at most `rel_tol`
            (percentage or fraction).
        """
        if isinstance(self.compressor, VariableSpeedCompressor):
            self.compressor.speed = self.n_cmp
        self.compressor.dT_sh = self.dT_sh
        i = [0]
        sol = optimize.root(
            self._eq_root,
            args=(i,),
            x0=np.array([T_evp_ini.to('degC').m, T_cnd_ini.to('degC').m]),
            options=dict(xtol=rel_tol.to('frac').m, maxfev=i_max)
        )
        self._log((
            f"Rating finished after {i[0]} iterations: {sol.message}"
            ), logging.INFO
        )
        self.compressor.Te = Q_(sol.x[0], 'degC')
        self.compressor.Tc = Q_(sol.x[1], 'degC')

    def _eq_minimize(self, unknowns: np.ndarray, counter: list) -> float:
        i = counter[0]
        T_evp = Q_(unknowns[0], 'degC')
        T_cnd = Q_(unknowns[1], 'degC')
        self._log((
            f"Iteration {i + 1}: "
            "try with "
            f"T_evp = {T_evp:~P.3f} and "
            f"T_cnd = {T_cnd:~P.3f}"
            ), logging.INFO
        )
        self.compressor.Te = T_evp
        self.compressor.Tc = T_cnd
        cmp_m_dot_rfg = self.compressor.m_dot.to('kg / hr')
        cnd_rfg_in = self.compressor.discharge_gas
        self._log((
            f"Iteration {i + 1}: "
            f"compressor mass flow rate: {cmp_m_dot_rfg:~P.3f}, "
            f"compressor power: {self.compressor.Wc_dot.to('kW'):~P.3f}"
            ), logging.DEBUG
        )
        self._log((
            f"Iteration {i + 1}: "
            "refrigerant at condenser entry: "
            f"{cnd_rfg_in.T.to('degC'):~P.3f}, "
            f"{cnd_rfg_in.P.to('bar'):~P.3f}"
            ), logging.DEBUG
        )
        with warnings.catch_warnings(category=CondenserWarning):
            warnings.showwarning = self._warn_handler
            r_cnd = self.condenser(
                m_dot_air=self.cnd_m_dot_air,
                m_dot_rfg=cmp_m_dot_rfg,
                air_in=self.cnd_air_in,
                rfg_in=cnd_rfg_in,
            )
        self._log((
            f"Iteration {i + 1}: "
            "refrigerant at condenser exit: "
            f"{r_cnd.rfg_out.T.to('degC'):~P.3f}, "
            f"{r_cnd.rfg_out.P.to('bar'):~P.3f}, "
            f"{r_cnd.rfg_out.h.to('kJ / kg'):~P.3f}, "
            f"{r_cnd.rfg_out.x.to('frac'):~P.3f}, "
            f"{r_cnd.rfg_out.phase}"
            ), logging.DEBUG
        )
        self._log((
            f"Iteration {i + 1}: "
            "heat rejection rate: "
            f"{r_cnd.Q_dot.to('kW'):~P.3f}"
            ), logging.DEBUG
        )
        evp_rfg_in = self.Refrigerant(
            h=r_cnd.rfg_out.h,
            P=self.compressor.Pe
        )
        self._log((
            f"Iteration {i + 1}: "
            f"refrigerant at evaporator entry: "
            f"{evp_rfg_in.T.to('degC'):~P.3f}, "
            f"{evp_rfg_in.P.to('bar'):~P.3f}, "
            f"{evp_rfg_in.h.to('kJ / kg'):~P.3f}, "
            f"{evp_rfg_in.x.to('frac'):~P.3f}, "
            f"{evp_rfg_in.phase}"
            ), logging.DEBUG
        )
        r_evp = self.evaporator(
            air_in=self.evp_air_in,
            m_dot_air=self.evp_m_dot_air,
            rfg_in=evp_rfg_in,
            dT_rfg_sh=self.dT_sh,
            m_dot_rfg_ini=cmp_m_dot_rfg
        )
        self._log((
            f"Iteration {i + 1}: "
            f"refrigerant at evaporator exit: "
            f"{r_evp.rfg_out.T.to('degC'):~P.3f}, "
            f"{r_evp.rfg_out.P.to('bar'):~P.3f}"
        ), logging.DEBUG
        )
        self._log((
            f"Iteration {i + 1}: "
            f"heat absorption rate: "
            f"{r_evp.Q_dot.to('kW'):~P.3f}"
            ), logging.DEBUG
        )
        evp_m_dot_rfg = r_evp.m_dot_rfg.to('kg / hr')
        self._log((
            f"Iteration {i + 1}: "
            f"evaporator mass flow rate: "
            f"{evp_m_dot_rfg:~P.3f}"
            ), logging.DEBUG
        )
        dev = abs(evp_m_dot_rfg - cmp_m_dot_rfg)
        self._log((
            f"Iteration {i + 1}: deviation "
            f"with {T_evp:~P.3f} and {T_cnd:~P.3f}: "
            f"{dev:~P.3f}"
            ), logging.INFO
        )
        counter[0] += 1
        return dev.m

    def rate_min(
        self,
        T_evp_ini: Quantity,
        T_cnd_ini: Quantity,
        i_max: int = 50,
        T_tol: Quantity = Q_(0.1, 'K'),
        m_dot_tol: Quantity = Q_(0.1, 'kg / hr')
    ) -> None:
        """Alternative implementation of `rate` using `scipy.optimize.minimize`.

        The algorithm tries to find the evaporation and condensation temperature
        for which the difference between the mass flow rate of refrigerant
        according to the compressor model and the mass flow rate of refrigerant
        according to the evaporator model is minimal.

        Parameters
        ----------
        T_evp_ini:
            Initial guess of the evaporation temperature.
        T_cnd_ini:
            Initial guess of the condensation temperature.
        i_max:
            Maximum number of iterations.
        T_tol:
            Absolute tolerance on the deviation of evaporation and condensation
            temperature between two consecutive iterations.
        m_dot_tol:
            Absolute tolerance on the refrigerant mass flow rate difference
            between two consecutive iterations.

        Notes
        -----
        Both tolerances must be fulfilled for the method to stop. However, the
        method will always terminate when the maximum number of iterations
        `i_max` has been reached before the tolerance conditions are met.

        Returns
        -------
        None
        """
        if isinstance(self.compressor, VariableSpeedCompressor):
            self.compressor.speed = self.n_cmp
        self.compressor.dT_sh = self.dT_sh
        i = [0]
        # Set physical bounds on the values for evaporation and condensation
        # temperature:
        T_crit = self.Refrigerant.critical_point.T
        bounds = (
            (-np.inf, self.evp_air_in.Tdb.to('degC').m),  # T_evp
            (self.cnd_air_in.Tdb.to('degC').m, T_crit.to('degC').m)  # T_cnd
        )
        sol = optimize.minimize(
            self._eq_minimize,
            args=(i,),
            x0=np.array([T_evp_ini.to('degC').m, T_cnd_ini.to('degC').m]),
            method='Nelder-Mead',
            bounds=bounds,
            options=dict(
                maxiter=i_max,
                maxfev=i_max,
                xatol=T_tol.to('K').m,
                fatol=m_dot_tol.to('kg / hr').m
            )
        )
        self._log((
            f"Rating finished after {i[0]} iterations: {sol.message}"
            ), logging.INFO
        )
        self.compressor.Te = Q_(sol.x[0], 'degC')
        self.compressor.Tc = Q_(sol.x[1], 'degC')

    def balance_by_speed(
        self,
        T_evp: Quantity,
        T_cnd: Quantity,
        i_max: int = 50,
        tol: Quantity = Q_(0.001, 'kg / s')
    ) -> Quantity | None:
        """Finds the compressor speed for which the mass flow rate of refrigerant
        displaced by the compressor balances the mass flow rate of refrigerant
        let through by the expansion device in order to maintain the set
        degree of superheat at the given evaporation and condensing temperature
        and at the operating conditions set with `set_operating_conditions`.

        Parameters
        ----------
        T_evp:
            Evaporation temperature at which the mass flow rates must be
            balanced.
        T_cnd:
            Condensing temperature at which the mass flow rates must be
            balanced.
        i_max:
            Maximum number of iterations.
        tol:
            Acceptable deviation between the mass flow rates.

        Notes
        -----
        This method can only be used with a variable speed compressor. Also,
        the minimum and maximum compressor speed must have been set when
        creating the `SingleStageVaporCompressionMachine` object.

        Returns
        -------
        Compressor speed at balance point.
        """
        c1 = isinstance(self.compressor, VariableSpeedCompressor)
        c2 = (self.n_cmp_min is not None) and (self.n_cmp_max is not None)
        if c1 and c2:
            def _eq(n_cmp: float) -> float:
                self._log((
                    f"Try with: "
                    f"{n_cmp:.3f} rpm"
                    ), logging.INFO
                )
                self.compressor.speed = Q_(n_cmp, '1 / min')
                cmp_m_dot_rfg = self.compressor.m_dot.to('kg / hr')
                self._log((
                    "Mass flow rate compressor: "
                    f"{cmp_m_dot_rfg.to('kg / hr'):~P.3f}"
                    ), logging.INFO
                )
                r_cnd = self.condenser(
                    m_dot_air=self.cnd_m_dot_air,
                    m_dot_rfg=cmp_m_dot_rfg,
                    air_in=self.cnd_air_in,
                    rfg_in=self.compressor.discharge_gas
                )
                evp_rfg_in = self.Refrigerant(
                    h=r_cnd.rfg_out.h,
                    P=self.compressor.Pe
                )
                self._log((
                    "Temperature refrigerant entering evaporator: "
                    f"{evp_rfg_in.T.to('degC'):~P.3f}"
                    ), logging.INFO
                )
                self._log((
                    "Vapor quality refrigerant entering evaporator: "
                    f"{evp_rfg_in.x.to('frac'):~P.5f}"
                    ), logging.INFO
                )
                r_evp = self.evaporator(
                    air_in=self.evp_air_in,
                    m_dot_air=self.evp_m_dot_air,
                    rfg_in=evp_rfg_in,
                    dT_rfg_sh=self.dT_sh,
                    m_dot_rfg_ini=cmp_m_dot_rfg
                )
                evp_m_dot_rfg = r_evp.m_dot_rfg.to('kg / hr')
                self._log((
                        "Mass flow rate evaporator: "
                        f"{evp_m_dot_rfg.to('kg / hr'):~P.3f}"
                    ), logging.INFO
                )
                dev = cmp_m_dot_rfg - evp_m_dot_rfg
                return dev.m

            # initialize compressor:
            self.compressor.Te = T_evp
            self.compressor.Tc = T_cnd
            self.compressor.dT_sh = self.dT_sh
            # find balanced speed between min and max compressor speed:
            mag_n_cmp_min = self.n_cmp_min.to('1 / min').m
            mag_n_cmp_max = self.n_cmp_max.to('1 / min').m
            sol = optimize.root_scalar(
                _eq,
                bracket=[mag_n_cmp_min, mag_n_cmp_max],
                xtol=tol.to('kg / hr').m,
                maxiter=i_max
            )
            n_cmp = Q_(sol.root, '1 / min')
            return n_cmp
        return None

    @property
    def Qc_dot(self) -> Quantity:
        """Get cooling capacity of vapor compression system under the conditions
        set at instantiation of the `VaporCompressionSystem`-object.
        """
        return self.evaporator.Q_dot

    @property
    def Wc_dot(self) -> Quantity:
        """Get compressor input power under the conditions set at instantiation
        of the `VaporCompressionSystem`-object.
        """
        return self.compressor.Wc_dot

    @property
    def m_dot(self) -> Quantity:
        """Get mass flow rate of refrigerant according to compressor under the
        conditions set at instantiation of the `VaporCompressionSystem`-object.
        """
        return self.compressor.m_dot

    @property
    def Qh_dot(self) -> Quantity:
        """Get condenser heat rejection under the conditions
        set at instantiation of the `VaporCompressionSystem`-object.
        """
        return self.condenser.Q_dot

    @property
    def COP(self) -> Quantity:
        """Get COP of vapor compression system under the conditions
        set at instantiation of the `VaporCompressionSystem`-object.
        """
        return self.compressor.COP

    @property
    def Te(self) -> Quantity:
        """Get evaporator temperature of vapor compression system under the
        conditions set at instantiation of the `VaporCompressionSystem`-object.
        """
        return self.evaporator.boiling_part.rfg_sat_vap_out.T

    @property
    def Tc(self) -> Quantity:
        """Get condenser temperature of vapor compression system under the
        conditions set at instantiation of the `VaporCompressionSystem`-object.
        """
        return self.condenser.condensing_part.rfg_sat_liq_out.T

    @property
    def suction_gas(self) -> FluidState:
        """Get state of suction gas under the conditions set at instantiation
        of the `VaporCompressionSystem`-object.
        """
        return self.evaporator.rfg_out

    @property
    def discharge_gas(self) -> FluidState:
        """Get state of discharge gas under the conditions set at instantiation
        of the `VaporCompressionSystem`-object.
        """
        return self.condenser.rfg_in

    @property
    def liquid(self) -> FluidState:
        """Get state of liquid under the conditions set at instantiation
        of the `VaporCompressionSystem`-object.
        """
        return self.condenser.rfg_out

    @property
    def mixture(self) -> FluidState:
        """Get state of mixture under the conditions set at instantiation
        of the `VaporCompressionSystem`-object.
        """
        return self.evaporator.rfg_in

    @property
    def isentropic_efficiency(self) -> Quantity:
        """Get the isentropic efficiency of the compressor."""
        return self.compressor.eta_is

    @property
    def super_heat(self) -> Quantity:
        """Get the degree of superheat of the suction gas."""
        return self.evaporator.dT_rfg_sh

    @property
    def sub_cooling(self) -> Quantity:
        """Get the degree of subcooling of the liquid at the inlet of
        the expansion device.
        """
        return self.condenser.dT_sc

    @property
    def Pe(self) -> Quantity:
        return self.compressor.Pe

    @property
    def Pc(self) -> Quantity:
        return self.compressor.Pc

    def check_mass_balance(self) -> tuple[Quantity, Quantity]:
        """Returns the absolute error and relative error between the mass flow
        rate of refrigerant according to the compressor model (computed by the
        correlation in the compressor model with the resulting evaporation and
        condensation temperature) and the mass flow rate of refrigerant
        according to the evaporator model (computed by the evaporator model in
        order to maintain the degree of superheat set on the expansion device).
        The relative error is the ratio of the absolute error to the mass flow
        rate of refrigerant according to the compressor model. It expresses
        the deviation of the refrigerant mass flow rate according to the
        evaporator model with respect to the refrigerant mass flow rate according
        to the compressor model as a fraction (percentage).
        """
        abs_err = abs(self.evaporator.m_dot_rfg - self.compressor.m_dot)
        rel_err = abs_err / abs(self.compressor.m_dot)
        return abs_err, rel_err.to('pct')

    def check_energy_balance(self) -> tuple[Quantity, Quantity]:
        """Returns the absolute error and relative error between the
        refrigeration capacity (heat absorption rate) according to the
        condensing unit model (computed by taking the difference between the
        heat rejection rate of the condenser and the compressor power) and
        the refrigeration capacity according to the evaporator model.
        The relative error is the ratio of the absolute error to the
        refrigeration capacity according to the compressor model. It expresses
        the deviation of the refrigerant capacity according to the evaporator
        model with respect to the refrigerant capacity according to the
        compressor model as a fraction (percentage).
        """
        Q_evp_cmp_model = self.condenser.Q_dot - self.compressor.Wc_dot
        Q_evp_evp_model = self.evaporator.Q_dot
        abs_err = abs(Q_evp_evp_model - Q_evp_cmp_model)
        rel_err = abs_err / abs(Q_evp_cmp_model)
        return abs_err, rel_err.to('pct')
