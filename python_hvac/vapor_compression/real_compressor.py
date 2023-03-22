from typing import List, Optional, Dict
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.optimize import fsolve
from .. import Quantity
from ..fluids import Fluid, FluidState


Q_ = Quantity


class FixedSpeedCorrelation:

    def __init__(self, param: str, file: Path):
        """
        Create `Correlation` object.

        Parameters
        ----------
        param: str ['Qc_dot', 'Wc_dot', 'I', 'm_dot']
            The quantity to which the correlation applies.
        file: Path
            The file path to the csv-file with polynomial coefficients.
        """
        self.param = param
        self.df = pd.read_csv(file)
        self.df.set_index('Unnamed: 0', drop=True, inplace=True)
        self.C: List[float] = self._get_coefficients()

    def _get_coefficients(self) -> List[float]:
        """Returns the polynomial coefficients for the given quantity set by
        `param` at instantiation of the `Correlation`-object."""
        return self.df.loc[self.param].values

    def __call__(self, *args, **kwargs) -> float:
        """Get the value of the quantity X at the given evaporator
        temperature `Te` and given condenser temperature `Tc`."""
        Te = args[0] or kwargs.get('Te', 0.0)
        Tc = args[1] or kwargs.get('Tc', 0.0)
        if len(self.C) == 10:
            return self.correlationC10(Te, Tc)

    def correlationC10(self, Te: float, Tc: float) -> float:
        X = self.C[0]
        X += self.C[1] * Te
        X += self.C[2] * Tc
        X += self.C[3] * (Te ** 2)
        X += self.C[4] * Te * Tc
        X += self.C[5] * (Tc ** 2)
        X += self.C[6] * (Te ** 3)
        X += self.C[7] * Tc * (Te ** 2)
        X += self.C[8] * Te * (Tc ** 2)
        X += self.C[9] * (Tc ** 3)
        return X


class VariableSpeedCorrelation(FixedSpeedCorrelation):

    def __call__(self, *args, **kwargs) -> float:
        """Get the value of the quantity X at the given evaporator
        temperature `Te`, given condenser temperature `Tc`, and given
        compressor speed `speed`."""
        Te = args[0] or kwargs.get('Te', 0.0)
        Tc = args[1] or kwargs.get('Tc', 0.0)
        speed = args[2] or kwargs.get('speed', 0.0)
        if len(self.C) == 20:
            return self.correlationC20(Te, Tc, speed)
        if len(self.C) == 30:
            return self.correlationC30(Te, Tc, speed)

    def correlationC20(self, Te: float, Tc: float, speed: float) -> float:
        X = self.C[0]
        X += self.C[1] * Te
        X += self.C[2] * Tc
        X += self.C[3] * speed
        X += self.C[4] * Te * Tc
        X += self.C[5] * Te * speed
        X += self.C[6] * Tc * speed
        X += self.C[7] * Te ** 2
        X += self.C[8] * Tc ** 2
        X += self.C[9] * speed ** 2
        X += self.C[10] * Te * Tc * speed
        X += self.C[11] * Te ** 2 * Tc
        X += self.C[12] * Te ** 2 * speed
        X += self.C[13] * Te ** 3
        X += self.C[14] * Te * Tc ** 2
        X += self.C[15] * Tc ** 2 * speed
        X += self.C[16] * Tc ** 3
        X += self.C[17] * Te * speed ** 2
        X += self.C[18] * Tc * speed ** 2
        X += self.C[19] * speed ** 3
        return X

    def correlationC30(self, Te: float, Tc: float, speed: float) -> float:
        X = self.C[0]
        X += self.C[1] * Te
        X += self.C[2] * Tc
        X += self.C[3] * Te ** 2
        X += self.C[4] * Tc ** 2
        X += self.C[5] * Te * Tc * speed ** 2
        X += self.C[6] * Te ** 2 * Tc * speed ** 2
        X += self.C[7] * Te * Tc ** 2 * speed ** 2
        X += self.C[8] * Te * Tc * speed
        X += self.C[9] * Te ** 2 * Tc * speed
        X += self.C[10] * Te * Tc ** 2 * speed
        X += self.C[11] * Te * Tc
        X += self.C[12] * Te ** 2 * Tc
        X += self.C[13] * Te * Tc ** 2
        X += self.C[14] * Te ** 3
        X += self.C[15] * Tc ** 3
        X += self.C[16] * speed
        X += self.C[17] * Te * speed
        X += self.C[18] * Tc * speed
        X += self.C[19] * Te ** 2 * speed
        X += self.C[20] * Tc ** 2 * speed
        X += self.C[21] * Te ** 3 * speed
        X += self.C[22] * Tc ** 3 * speed
        X += self.C[23] * speed ** 2
        X += self.C[24] * Te * speed ** 2
        X += self.C[25] * Tc * speed ** 2
        X += self.C[26] * Te ** 2 * speed ** 2
        X += self.C[27] * Tc ** 2 * speed ** 2
        X += self.C[28] * Te ** 3 * speed ** 2
        X += self.C[29] * Tc ** 3 * speed ** 2
        return X


class FixedSpeedCompressor:
    params = ['Qc_dot', 'Wc_dot', 'm_dot', 'T_dis']
    units = {
        'Qc_dot': 'kW',
        'Wc_dot': 'kW',
        'm_dot': 'g / s',
        'speed': '1 / min'
    }

    def __init__(
        self,
        coeff_file: Path,
        dT_sh: Quantity,
        dT_sc: Quantity,
        refrigerant_type: Fluid,
        units: Dict[str, str] | None = None
    ) -> None:
        """
        Create `FixedSpeedCompressor`-model from polynomial coefficients that define
        the correlations for cooling capacity (`Qc_dot`), compressor input power
        (`Wc_dot`) and refrigerant mass flow rate (`m_dot`) as a function of
        evaporator temperature and condenser temperature.

        Parameters
        ----------
        coeff_file: Path
            File path to csv-file with polynomial coefficients of compressor.
        dT_sh: Quantity
            Amount of superheat for which the polynomial coefficients are valid.
        dT_sc: Quantity
            Amount of subcooling for which the polynomial coefficients are valid.
        refrigerant_type: Fluid
            Refrigerant for which the polynomial coefficients are valid.
        units: Dict[str, str], optional
            The units used in the coefficients file.
        """
        self.dT_sh = dT_sh
        self.dT_sc = dT_sc
        self.refrigerant_type = refrigerant_type
        self._Te: float = float('nan')
        self._Tc: float = float('nan')
        self._set_correlations(coeff_file)
        if units is not None:
            self.units.update(units)

    def _set_correlations(self, coeff_file: Path):
        # create a dictionary to hold the `Correlation`-object for each of the
        # possible quantities (`Qc_dot`, `Wc_dot`, `m_dot`, `T_dis`).
        self._correlations = {}
        for param in self.params:
            try:
                self._correlations[param] = FixedSpeedCorrelation(param, coeff_file)
            except KeyError:
                pass

    @property
    def Te(self) -> Quantity:
        return Q_(self._Te, 'degC')

    @Te.setter
    def Te(self, Te: Quantity) -> None:
        """Set evaporator temperature."""
        self._Te = Te.to('degC').m

    @property
    def Tc(self) -> Quantity:
        return Q_(self._Tc, 'degC')

    @Tc.setter
    def Tc(self, Tc: Quantity) -> None:
        """Set condenser temperature."""
        self._Tc = Tc.to('degC').m

    @property
    def Qc_dot(self) -> Quantity:
        """Get cooling capacity at the set evaporator temperature and condenser
        temperature."""
        Qc_dot = self._correlations['Qc_dot'](self._Te, self._Tc)
        return Q_(Qc_dot, self.units['Qc_dot'])

    @property
    def Wc_dot(self) -> Quantity:
        """Get compressor input power at the set evaporator temperature and
        condenser temperature."""
        Wc_dot = self._correlations['Wc_dot'](self._Te, self._Tc)
        return Q_(Wc_dot, self.units['Wc_dot'])

    @property
    def m_dot(self) -> Quantity:
        """Get mass flow rate at the set evaporator temperature and condenser
        temperature."""
        m_dot = self._correlations['m_dot'](self._Te, self._Tc)
        return Q_(m_dot, self.units['m_dot'])

    @property
    def T_dis(self) -> Quantity:
        """Get discharge temperature at condenser entrance."""
        T_dis = self._correlations['T_dis'](self._Te, self._Tc)
        return Q_(T_dis, 'degC')

    @property
    def refrigeration_effect(self) -> Quantity:
        q = self.Qc_dot / self.m_dot
        return q

    @property
    def suction_gas(self) -> FluidState:
        """Get state of suction gas at evaporator exit."""
        if self.dT_sh == 0:
            suction_gas = self.refrigerant_type(T=self.Te, x=Q_(100, 'pct'))
        else:
            P_eva = self.refrigerant_type(T=self.Te, x=Q_(100, 'pct')).P
            T_suc = self.Te.to('K') + self.dT_sh.to('K')
            suction_gas = self.refrigerant_type(T=T_suc, P=P_eva)
        return suction_gas

    @property
    def mixture(self) -> FluidState:
        """Get state of mixture at evaporator entrance."""
        P_eva = self.suction_gas.P
        h_mix = self.liquid.h
        try:
            mixture = self.refrigerant_type(P=P_eva, h=h_mix)
        except IndexError:  # refrigerant is mixture
            mixture = self.refrigerant_type(P=P_eva, h=h_mix, x=Q_(0, 'pct'))
        return mixture

    @property
    def liquid(self) -> FluidState:
        """Get state of liquid at condenser exit."""
        P_con = self.refrigerant_type(T=self.Tc, x=Q_(0, 'pct')).P
        T_liq = self.Tc.to('K') - self.dT_sc.to('K')
        if self.dT_sc == 0:
            liquid = self.refrigerant_type(P=P_con, x=Q_(0, 'pct'))
        else:
            liquid = self.refrigerant_type(P=P_con, T=T_liq)
        return liquid

    @property
    def discharge_gas(self) -> FluidState:
        """Get state of discharge gas at condenser entrance."""
        P_con = self.refrigerant_type(T=self.Tc, x=Q_(0, 'pct')).P
        wc = self.Wc_dot / self.m_dot
        h_dis = self.suction_gas.h + wc
        try:
            discharge_gas = self.refrigerant_type(P=P_con, h=h_dis)
        except IndexError:  # refrigerant is mixture
            discharge_gas = self.refrigerant_type(P=P_con, h=h_dis, T=self.Tc)
        return discharge_gas

    @property
    def Wis_dot(self) -> Quantity:
        """Get isentropic compressor power at the set evaporator temperature and
        condenser temperature."""
        P_con = self.refrigerant_type(T=self.Tc, x=Q_(0, 'pct')).P
        h_suc = self.suction_gas.h
        s_suc = self.suction_gas.s
        try:
            h_is_dis = self.refrigerant_type(P=P_con, s=s_suc).h
        except IndexError:  # refrigerant is mixture
            h_is_dis = self.refrigerant_type(P=P_con, s=s_suc, T=self.Tc).h
        Wis_dot = self.m_dot * (h_is_dis - h_suc)
        return Wis_dot

    @property
    def eta_is(self) -> Quantity:
        """Get isentropic efficiency at the set evaporator temperature and
        condenser temperature."""
        return self.Wis_dot / self.Wc_dot

    @property
    def Qh_dot(self) -> Quantity:
        """Get condenser heat rejection at the set evaporator temperature and
        condenser temperature."""
        return self.Qc_dot + self.Wc_dot

    @property
    def COP(self) -> Quantity:
        """Get COP at the set evaporator temperature and condenser temperature."""
        return self.Qc_dot / self.Wc_dot

    def get_refrigerant_cycle(self, units: Optional[Dict[str, str]] = None) -> pd.DataFrame:
        """Returns a Pandas DataFrame with state properties at evaporator inlet, evaporator outlet
        condenser inlet and condenser outlet, assuming a standard vapor compression cycle without
        pressure losses."""
        if units is None:
            units = {
                'T': 'degC',
                'P': 'bar',
                'rho': 'kg / m**3',
                'h': 'kJ / kg',
                's': 'kJ / kg / K'
            }
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
        coeff_file: Path,
        dT_sh: Quantity,
        dT_sc: Quantity,
        refrigerant_type: Fluid,
        units: Dict[str, str] | None = None
    ) -> None:
        """
        Create `VariableSpeedCompressor`-model from polynomial coefficients that define
        the correlations for cooling capacity (`Qc_dot`), compressor input power
        (`Wc_dot`) and refrigerant mass flow rate (`m_dot`) as a function of
        evaporator temperature, condenser temperature and compressor speed.

        Parameters
        ----------
        coeff_file: Path
            File path to csv-file with polynomial coefficients of compressor.
        dT_sh: Quantity
            Amount of superheat for which the polynomial coefficients are valid.
        dT_sc: Quantity
            Amount of subcooling for which the polynomial coefficients are valid.
        refrigerant_type: Fluid
            Refrigerant for which the polynomial coefficients are valid.
        """
        super().__init__(coeff_file, dT_sh, dT_sc, refrigerant_type, units)
        self._speed = float('nan')

    def _set_correlations(self, coeff_file: Path):
        # create a dictionary to hold the `Correlation`-object for each of the
        # possible quantities (`Qc_dot`, `Wc_dot` and `m_dot`).
        self._correlations = {}
        for param in self.params:
            try:
                self._correlations[param] = VariableSpeedCorrelation(param, coeff_file)
            except KeyError:
                pass

    @property
    def speed(self) -> Quantity:
        return Q_(self._speed, self.units['speed'])

    @speed.setter
    def speed(self, s: Quantity) -> None:
        """Set compressor speed."""
        self._speed = s.to(self.units['speed']).m

    @property
    def Qc_dot(self) -> Quantity:
        """Get cooling capacity at the set evaporator temperature, condenser
        temperature, and compressor speed."""
        Qc_dot = self._correlations['Qc_dot'](self._Te, self._Tc, self._speed)
        return Q_(Qc_dot, self.units['Qc_dot'])

    @property
    def Wc_dot(self) -> Quantity:
        """Get compressor input power at the set evaporator temperature, condenser
        temperature, and compressor speed."""
        Wc_dot = self._correlations['Wc_dot'](self._Te, self._Tc, self._speed)
        return Q_(Wc_dot, self.units['Wc_dot'])

    @property
    def m_dot(self) -> Quantity:
        """Get mass flow rate at the set evaporator temperature, condenser
        temperature, and compressor speed."""
        try:
            m_dot = self._correlations['m_dot'](self._Te, self._Tc, self._speed)
            return Q_(m_dot, self.units['m_dot'])
        except KeyError:
            q_re = self.suction_gas.h - self.liquid.h
            m_dot = self.Qc_dot / q_re
            return m_dot

    @property
    def T_dis(self) -> Quantity:
        """Get discharge temperature at the set evaporator temperature, condenser
        temperature, and compressor speed."""
        T_dis = self._correlations['T_dis'](self._Te, self._Tc, self._speed)
        return Q_(T_dis, 'degC')

    def get_compressor_speed(self, Qc_dot: Quantity) -> Quantity:
        """Returns the compressor speed at which the cooling capacity is
        equal to `Qc_dot` at the set evaporator temperature and condenser
        temperature."""
        _Qc_dot = Qc_dot.to(self.units['Qc_dot']).m

        def eq(unknowns: np.ndarray) -> np.ndarray:
            speed = unknowns[0]
            lhs = _Qc_dot
            rhs = self._correlations['Qc_dot'](self._Te, self._Tc, speed)
            out = lhs - rhs
            return np.array([out])

        speed_ini = 0.0
        roots = fsolve(eq, np.array([speed_ini]))
        speed = Q_(roots[0], self.units['speed'])
        return speed
