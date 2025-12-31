from pathlib import Path
import pandas as pd
import numpy as np
from scipy.optimize import fsolve
from .. import Quantity
from ..fluids import Fluid, FluidState, CoolPropMixtureError


Q_ = Quantity


class FixedSpeedCorrelation:

    def __init__(self, param: str, file: Path):
        """
        Create `Correlation` object.

        Parameters
        ----------
        param: str
            The quantity to which the correlation applies.
        file: Path
            The file path to the csv-file with polynomial coefficients.
        """
        self.param = param
        self.df = pd.read_csv(file)
        self.df.set_index('Unnamed: 0', drop=True, inplace=True)
        self.C: list[float] = self._get_coefficients()

    def _get_coefficients(self) -> list[float]:
        """Returns the polynomial coefficients for the given quantity indicated
        by `param` at instantiation of the `Correlation`-object.
        """
        return self.df.loc[self.param].values.tolist()

    def __call__(self, *args, **kwargs) -> float | None:
        """Returns the value of quantity "X" at the given evaporation temperature
        `T_evp` and the given condensing temperature `T_cnd`.
        """
        T_evp = args[0] or kwargs.get('T_evp', 0.0)
        T_cnd = args[1] or kwargs.get('T_cnd', 0.0)
        if len(self.C) == 10:
            return self.correlationC10(T_evp, T_cnd)
        return None

    def correlationC10(self, T_evp: float, T_cnd: float) -> float:
        X = self.C[0]
        X += self.C[1] * T_evp
        X += self.C[2] * T_cnd
        X += self.C[3] * (T_evp ** 2)
        X += self.C[4] * T_evp * T_cnd
        X += self.C[5] * (T_cnd ** 2)
        X += self.C[6] * (T_evp ** 3)
        X += self.C[7] * T_cnd * (T_evp ** 2)
        X += self.C[8] * T_evp * (T_cnd ** 2)
        X += self.C[9] * (T_cnd ** 3)
        return X


class VariableSpeedCorrelation(FixedSpeedCorrelation):

    def __call__(self, *args, **kwargs) -> float | None:
        """Returns the value of the quantity "X" at the given evaporation
        temperature `T_evp`, given condensing temperature `T_cnd`, and given
        compressor n `n`."""
        T_evp = args[0] or kwargs.get('T_evp', 0.0)
        T_cnd = args[1] or kwargs.get('T_cnd', 0.0)
        n = args[2] or kwargs.get('n', 0.0)
        if len(self.C) == 20:
            return self.correlationC20(T_evp, T_cnd, n)
        if len(self.C) == 30:
            return self.correlationC30(T_evp, T_cnd, n)
        return None

    def correlationC20(self, T_evp: float, T_cnd: float, n: float) -> float:
        X = self.C[0]
        X += self.C[1] * T_evp
        X += self.C[2] * T_cnd
        X += self.C[3] * n
        X += self.C[4] * T_evp * T_cnd
        X += self.C[5] * T_evp * n
        X += self.C[6] * T_cnd * n
        X += self.C[7] * T_evp ** 2
        X += self.C[8] * T_cnd ** 2
        X += self.C[9] * n ** 2
        X += self.C[10] * T_evp * T_cnd * n
        X += self.C[11] * T_evp ** 2 * T_cnd
        X += self.C[12] * T_evp ** 2 * n
        X += self.C[13] * T_evp ** 3
        X += self.C[14] * T_evp * T_cnd ** 2
        X += self.C[15] * T_cnd ** 2 * n
        X += self.C[16] * T_cnd ** 3
        X += self.C[17] * T_evp * n ** 2
        X += self.C[18] * T_cnd * n ** 2
        X += self.C[19] * n ** 3
        return X

    def correlationC30(self, T_evp: float, T_cnd: float, n: float) -> float:
        X = self.C[0]
        X += self.C[1] * T_evp
        X += self.C[2] * T_cnd
        X += self.C[3] * T_evp ** 2
        X += self.C[4] * T_cnd ** 2
        X += self.C[5] * T_evp * T_cnd * n ** 2
        X += self.C[6] * T_evp ** 2 * T_cnd * n ** 2
        X += self.C[7] * T_evp * T_cnd ** 2 * n ** 2
        X += self.C[8] * T_evp * T_cnd * n
        X += self.C[9] * T_evp ** 2 * T_cnd * n
        X += self.C[10] * T_evp * T_cnd ** 2 * n
        X += self.C[11] * T_evp * T_cnd
        X += self.C[12] * T_evp ** 2 * T_cnd
        X += self.C[13] * T_evp * T_cnd ** 2
        X += self.C[14] * T_evp ** 3
        X += self.C[15] * T_cnd ** 3
        X += self.C[16] * n
        X += self.C[17] * T_evp * n
        X += self.C[18] * T_cnd * n
        X += self.C[19] * T_evp ** 2 * n
        X += self.C[20] * T_cnd ** 2 * n
        X += self.C[21] * T_evp ** 3 * n
        X += self.C[22] * T_cnd ** 3 * n
        X += self.C[23] * n ** 2
        X += self.C[24] * T_evp * n ** 2
        X += self.C[25] * T_cnd * n ** 2
        X += self.C[26] * T_evp ** 2 * n ** 2
        X += self.C[27] * T_cnd ** 2 * n ** 2
        X += self.C[28] * T_evp ** 3 * n ** 2
        X += self.C[29] * T_cnd ** 3 * n ** 2
        return X


class FixedSpeedCompressor:
    params = ['Q_dot_evp', 'W_dot', 'm_dot', 'T_dis']
    units = {
        'Q_dot': 'kW',
        'W_dot': 'kW',
        'm_dot': 'kg / hr',
        'T_dis': 'degC',
        'n': 'rpm'
    }

    def __init__(
        self,
        coeff_file: Path | str,
        refrigerant: Fluid,
        dT_sh: Quantity = Q_(8, 'K'),
        dT_sc: Quantity = Q_(2, 'K'),
        units: dict[str, str] | None = None
    ) -> None:
        """Creates a `FixedSpeedCompressor` model using the polynomial
        coefficients in `coeff_file` for the correlations that return cooling
        capacity (`Q_dot_evp`), compressor input power (`W_dot`), refrigerant
        mass flow rate (`m_dot`), and discharge temperature (`T_dis`) as a
        function of evaporation temperature and condensing temperature.

        The polynomial coefficients of a compressor are to be retrieved through
        the selection software of compressor manufacturers. The values of the
        polynomial coefficients depend also on the amount of superheat and the
        amount of subcooling that were set in the selection program. These
        amounts need to be assigned to parameters `dT_sh` and `dT_sc`.

        Parameters
        ----------
        coeff_file:
            Path to the csv-file with the polynomial coefficients of the
            compressor.
        dT_sh:
            Amount of superheat for which the polynomial coefficients are valid.
        dT_sc:
            Amount of subcooling for which the polynomial coefficients are valid.
        refrigerant:
            Refrigerant for which the polynomial coefficients are valid.
        units: optional
            Dictionary with the units used in the coefficient file. The default
            units are: {'Q_dot': 'kW', 'W_dot': 'kW', 'm_dot': 'kg / hr', 
            'T_dis': 'degC'} where 'Q_dot' stands for cooling capacity, 'W_dot'
            for compressor power, 'm_dot' for refrigerant mass flow rate, and
            'T_dis' for discharge temperature.
        """
        self.dT_sh = dT_sh
        self.dT_sc = dT_sc
        self.refrigerant = refrigerant
        self._T_evp: float = float('nan')
        self._T_cnd: float = float('nan')
        self._set_correlations(coeff_file)
        if units is not None:
            self.units.update(units)

    def _set_correlations(self, coeff_file: Path):
        # Creates a dictionary that holds the `Correlation`-object for each of 
        # the possible quantities (`Q_dot_evp`, `W_dot`, `m_dot`, and `T_dis`).
        self.correlations = {}
        for param in self.params:
            try:
                self.correlations[param] = FixedSpeedCorrelation(param, coeff_file)
            except KeyError:
                pass

    @property
    def T_evp(self) -> Quantity:
        """Evaporation temperature."""
        return Q_(self._T_evp, 'degC')

    @T_evp.setter
    def T_evp(self, value: Quantity) -> None:
        self._T_evp = value.to('degC').m

    @property
    def T_cnd(self) -> Quantity:
        """Condensing temperature."""
        return Q_(self._T_cnd, 'degC')

    @T_cnd.setter
    def T_cnd(self, value: Quantity) -> None:
        self._T_cnd = value.to('degC').m

    @property
    def P_evp(self) -> Quantity:
        """Evaporation pressure."""
        return self.refrigerant(T=self.T_evp, x=Q_(1.0, 'frac')).P

    @property
    def P_cnd(self) -> Quantity:
        """Condensing pressure."""
        return self.refrigerant(T=self.T_cnd, x=Q_(0.0, 'frac')).P

    @property
    def Q_dot_evp(self) -> Quantity:
        """Cooling capacity at the set evaporation temperature and condensing
        temperature."""
        Q_dot_evp = self.correlations['Q_dot_evp'](self._T_evp, self._T_cnd)
        return Q_(Q_dot_evp, self.units['Q_dot'])

    @property
    def W_dot(self) -> Quantity:
        """Compressor input power at the set evaporation temperature and
        condensing temperature."""
        Wc_dot = self.correlations['W_dot'](self._T_evp, self._T_cnd)
        return Q_(Wc_dot, self.units['W_dot'])

    @property
    def m_dot(self) -> Quantity:
        """Mass flow rate at the set evaporation temperature and condensing
        temperature."""
        m_dot = self.correlations['m_dot'](self._T_evp, self._T_cnd)
        return Q_(m_dot, self.units['m_dot'])

    @property
    def T_dis(self) -> Quantity:
        """Discharge temperature at condenser entrance."""
        T_dis = self.correlations['T_dis'](self._T_evp, self._T_cnd)
        return Q_(T_dis, self.units['T_dis'])

    @property
    def q_rfg(self) -> Quantity:
        """Refrigeration effect."""
        q = self.Q_dot_evp / self.m_dot
        return q

    @property
    def suction_gas(self) -> FluidState:
        """Suction gas at evaporator outlet."""
        if self.dT_sh.m == 0:
            suction_gas = self.refrigerant(T=self.T_evp, x=Q_(100, 'pct'))
        else:
            P_evp = self.refrigerant(T=self.T_evp, x=Q_(100, 'pct')).P
            T_suc = self.T_evp.to('K') + self.dT_sh.to('K')
            suction_gas = self.refrigerant(T=T_suc, P=P_evp)
        return suction_gas

    @property
    def mixture(self) -> FluidState:
        """Liquid/vapor mixture at evaporator inlet."""
        P_evp = self.suction_gas.P
        h_mix = self.liquid.h
        try:
            mixture = self.refrigerant(P=P_evp, h=h_mix)
        except CoolPropMixtureError:  
            # refrigerant is a blend
            mixture = self.refrigerant(P=P_evp, h=h_mix, x=Q_(0, 'pct'))
        return mixture

    @property
    def liquid(self) -> FluidState:
        """Liquid at condenser outlet."""
        P_cnd = self.refrigerant(T=self.T_cnd, x=Q_(0, 'pct')).P
        T_liq = self.T_cnd.to('K') - self.dT_sc.to('K')
        if self.dT_sc == 0:
            liquid = self.refrigerant(P=P_cnd, x=Q_(0, 'pct'))
        else:
            liquid = self.refrigerant(P=P_cnd, T=T_liq)
        return liquid

    @property
    def discharge_gas(self) -> FluidState:
        """Discharge gas at condenser inlet."""
        P_cnd = self.refrigerant(T=self.T_cnd, x=Q_(100, 'pct')).P
        wc = self.W_dot / self.m_dot
        h_dis = self.suction_gas.h + wc
        try:
            discharge_gas = self.refrigerant(P=P_cnd, h=h_dis)
        except CoolPropMixtureError:  
            # refrigerant is a blend
            discharge_gas = self.refrigerant(P=P_cnd, h=h_dis, T=self.T_cnd)
        return discharge_gas

    @property
    def Wis_dot(self) -> Quantity:
        """Isentropic compressor power at the set evaporation temperature and
        condensing temperature.
        """
        P_cnd = self.refrigerant(T=self.T_cnd, x=Q_(0, 'pct')).P
        h_suc = self.suction_gas.h
        s_suc = self.suction_gas.s
        try:
            h_is_dis = self.refrigerant(P=P_cnd, s=s_suc).h
        except CoolPropMixtureError:  
            # refrigerant is mixture
            h_is_dis = self.refrigerant(P=P_cnd, s=s_suc, T=self.T_cnd).h
        Wis_dot = self.m_dot * (h_is_dis - h_suc)
        return Wis_dot

    @property
    def eta_is(self) -> Quantity:
        """Isentropic efficiency at the set evaporation temperature and
        condensing temperature.
        """
        return self.Wis_dot / self.W_dot

    @property
    def Q_dot_cnd(self) -> Quantity:
        """Condenser heat rejection rate at the set evaporation temperature and
        condensing temperature.
        """
        return self.Q_dot_evp + self.W_dot

    @property
    def COP(self) -> Quantity:
        """COP at the set evaporation temperature and condensing temperature."""
        return self.Q_dot_evp / self.W_dot

    def get_refrigerant_cycle(
        self, 
        units: dict[str, str] | None = None
    ) -> pd.DataFrame:
        """Returns a Pandas DataFrame with state properties at the evaporator 
        inlet, the evaporator outlet, the condenser inlet, and the condenser 
        outlet, assuming a standard vapor compression cycle without pressure 
        losses.
        
        Parameters
        ----------
        units:
            Dictionary with the units for displaying the quantities in the
            dataframe. The default units are: {'T': 'degC', 'P': 'bar', 
            'rho': 'kg / m**3', 'h': 'kJ / kg', 's': 'kJ / kg / K'}.
        """
        units_ = {
            'T': 'degC',
            'P': 'bar',
            'rho': 'kg / m**3',
            'h': 'kJ / kg',
            's': 'kJ / kg / K'
        }
        if units is not None:
            units_.update(units)
        units = units_
        columns = [
            f"T [{units['T']}]",
            f"P [{units['P']}]",
            f"rho [{units['rho']}]",
            f"h [{units['h']}]",
            f"s [{units['s']}]"
        ]
        index = [
            "evaporator inlet",
            "evaporator outlet",
            "condenser inlet",
            "condenser outlet"
        ]
        data = []
        for x in (self.mixture, self.suction_gas, self.discharge_gas, self.liquid):
            data.append([
                x.T.to(units['T']).m,
                x.P.to(units['P']).m,
                x.rho.to(units['rho']).m,
                x.h.to(units['h']).m,
                x.s.to(units['s']).m
            ])
        data = np.array(data)
        df = pd.DataFrame(data=data, columns=columns, index=index)
        return df


class VariableSpeedCompressor(FixedSpeedCompressor):

    def __init__(
        self,
        coeff_file: Path | str,
        refrigerant: Fluid,
        dT_sh: Quantity = Q_(0, 'K'),
        dT_sc: Quantity = Q_(0, 'K'),
        units: dict[str, str] | None = None
    ) -> None:
        """Creates a `VariableSpeedCompressor` model using the polynomial
        coefficients in `coeff_file` for the correlations that return cooling
        capacity (`Q_dot_evp`), compressor input power (`W_dot`), refrigerant
        mass flow rate (`m_dot`), and discharge temperature (`T_dis`) as a
        function of evaporation temperature, condensing temperature and
        compressor speed.

        The polynomial coefficients of a compressor are to be retrieved through
        the selection software of compressor manufacturers. The values of the
        polynomial coefficients depend also on the amount of superheat and the
        amount of subcooling that were set in the selection program. These
        amounts need to be assigned to parameters `dT_sh` and `dT_sc`.

        Parameters
        ----------
        coeff_file:
            File path to csv-file with the polynomial coefficients of the
            compressor.
        dT_sh:
            Amount of superheat for which the polynomial coefficients are valid.
        dT_sc:
            Amount of subcooling for which the polynomial coefficients are valid.
        refrigerant:
            Refrigerant for which the polynomial coefficients are valid.
        units: optional
            Dictionary with the units used in the coefficient file. The default
            units are: {'Q_dot': 'kW', 'W_dot': 'kW', 'm_dot': 'kg / hr', 
            'T_dis': 'degC', 'n': 'rpm'} where 'Q_dot' stands for cooling 
            capacity, 'W_dot' for compressor power, 'm_dot' for refrigerant mass
            flow rate, 'T_dis' for discharge temperature, and 'n' for compressor
            speed.
        """
        super().__init__(coeff_file, refrigerant, dT_sh, dT_sc, units)
        self._n = float('nan')

    def _set_correlations(self, coeff_file: Path):
        # creates a dictionary that holds the `Correlation`-objects for each of 
        # the possible quantities (`Q_dot_evp`, `W_dot`, `m_dot` and `T_dis`).
        self.correlations = {}
        for param in self.params:
            try:
                self.correlations[param] = VariableSpeedCorrelation(param, coeff_file)
            except KeyError:
                pass

    @property
    def speed(self) -> Quantity:
        """Compressor speed."""
        return Q_(self._n, self.units['n'])

    @speed.setter
    def speed(self, value: Quantity) -> None:
        self._n = value.to(self.units['n']).m

    @property
    def Q_dot_evp(self) -> Quantity:
        """Cooling capacity at the set evaporation temperature, condensing
        temperature, and compressor speed.
        """
        Q_dot_evp = self.correlations['Q_dot_evp'](self._T_evp, self._T_cnd, self._n)
        return Q_(Q_dot_evp, self.units['Q_dot'])

    @property
    def W_dot(self) -> Quantity:
        """Compressor input power at the set evaporation temperature, condensing
        temperature, and compressor speed.
        """
        W_dot = self.correlations['W_dot'](self._T_evp, self._T_cnd, self._n)
        return Q_(W_dot, self.units['W_dot'])

    @property
    def m_dot(self) -> Quantity:
        """Refrigerant mass flow rate at the set evaporation temperature,
        condensing temperature, and compressor speed.
        """
        try:
            m_dot = self.correlations['m_dot'](self._T_evp, self._T_cnd, self._n)
            return Q_(m_dot, self.units['m_dot'])
        except KeyError:
            q_re = self.suction_gas.h - self.liquid.h
            m_dot = self.Q_dot_evp / q_re
            return m_dot

    @property
    def T_dis(self) -> Quantity:
        """Discharge temperature at the set evaporation temperature, condensing
        temperature, and compressor speed.
        """
        T_dis = self.correlations['T_dis'](self._T_evp, self._T_cnd, self._n)
        return Q_(T_dis, 'degC')

    def compressor_speed(
        self,
        evp_Q_dot: Quantity | None = None,
        m_dot: Quantity | None = None
    ) -> Quantity | None:
        """Returns the compressor speed at which the cooling capacity or the
        refrigerant mass flow rate equals the value of `Q_dot_evp` or `m_dot`
        at the set evaporation temperature and condensing temperature.
        """
        if evp_Q_dot is not None:
            n = self.__solve1__(evp_Q_dot)
            return n
        if m_dot is not None:
            n = self.__solve2__(m_dot)
            return n
        return None

    def __solve1__(self, evp_Q_dot: Quantity) -> Quantity:
        _evp_Q_dot = evp_Q_dot.to(self.units['Q_dot']).m

        def eq(unknowns: np.ndarray) -> np.ndarray:
            n_ = unknowns[0]
            lhs = _evp_Q_dot
            rhs = self.correlations['Q_dot_evp'](self._T_evp, self._T_cnd, n_)
            out = lhs - rhs
            return np.array([out])

        n_ini = 0.0
        roots = fsolve(eq, np.array([n_ini]))
        n = Q_(roots[0], self.units['n'])
        return n

    def __solve2__(self, m_dot: Quantity) -> Quantity:
        _m_dot = m_dot.to(self.units['m_dot']).m

        def eq(unknowns: np.ndarray) -> np.ndarray:
            n_ = unknowns[0]
            lhs = _m_dot
            rhs = self.correlations['m_dot'](self._T_evp, self._T_cnd, n_)
            out = lhs - rhs
            return np.array([out])

        n_ini = 0.0
        roots = fsolve(eq, np.array([n_ini]))
        n = Q_(roots[0], self.units['n'])
        return n
