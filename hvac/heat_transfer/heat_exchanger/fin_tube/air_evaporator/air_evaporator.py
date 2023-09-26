"""
AIR EVAPORATOR MODEL
Plain fin-tube counter-flow heat exchanger.

Known inputs:
- State of air entering the evaporator.
- Mass flow rate of air entering the evaporator.
- State of refrigerant entering the evaporator.
- Required degree of refrigerant superheating at the evaporator outlet.

Purpose:
Determine the mass flow rate of refrigerant through the evaporator so that the
required degree of refrigerant superheating is reached. Once this mass flow
rate has been determined, other operational characteristics are also determined.
"""
import numpy as np
from scipy import optimize
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState, CP_HUMID_AIR
from hvac.heat_transfer.heat_exchanger.fin_tube.core import PlainFinTubeHeatExchangerCore
from hvac.heat_transfer.heat_exchanger import eps_ntu as dry
from hvac.heat_transfer.heat_exchanger import eps_ntu_wet as wet
from hvac.logging import ModuleLogger

Q_ = Quantity
logger = ModuleLogger.get_logger(__name__)
logger.setLevel(ModuleLogger.ERROR)


class EvaporatorError(Exception):
    pass


class SuperheatingError(EvaporatorError):
    pass


class BoilingError(EvaporatorError):
    pass


class SuperheatingRegion:
    
    def __init__(
        self,
        W_fro: Quantity,
        H_fro: Quantity,
        S_trv: Quantity,
        S_lon: Quantity,
        D_int: Quantity,
        D_ext: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity
    ) -> None:
        """
        Creates the superheating region of the plain fin-tube evaporator.

        Parameters
        ----------
        W_fro:
            Width of the frontal area of the evaporator.
        H_fro:
            Height of the frontal area.
        S_trv:
            Transversal pitch, i.e., the spacing between tubes in a single row.
        S_lon:
            Longitudinal pitch, i.e., the spacing between tubes of adjacent rows.
        D_int:
            Inside diameter of the tubes.
        D_ext:
            Outside diameter of the tubes.
        t_fin:
            Thickness of the plain fins.
        N_fin:
            Number of fins per unit length of tube (i.e., the inverse of fin
            density).
        k_fin:
            Thermal conductivity of the fin material.
        """
        self.core = PlainFinTubeHeatExchangerCore(
            L1=W_fro,
            L3=H_fro,
            S_t=S_trv,
            S_l=S_lon,
            D_i=D_int,
            D_o=D_ext,
            t_f=t_fin,
            N_f=N_fin,
            k_fin=k_fin
        )
        # Known parameters:
        self.air_in: HumidAir | None = None
        self.air_m_dot: Quantity | None = None
        self.rfg_in: FluidState | None = None
        self.rfg_out: FluidState | None = None

        # Parameters to be solved for:
        self.rfg_m_dot: Quantity | None = None
        self.air_out: HumidAir | None = None
        self.L_flow: Quantity | None = None

    def __fun__(
        self,
        L_flow: float,
        rfg_m_dot: Quantity,
        Q_dot: Quantity,
        air_mean: HumidAir,
        rfg_mean: FluidState,
        counter: list[int]
    ) -> float:
        """
        Calculates the deviation between the heat transfer rate through the
        heat exchanger and the heat absorbed by the refrigerant stream, which
        is also the heat rejected by the air stream.
        """
        # Set the parameters on the heat exchanger core needed to determine
        # the overall heat transfer conductance of the superheated region:
        self.core.L2 = Q_(L_flow, 'mm')
        self.core.ext.m_dot = self.air_m_dot
        self.core.int.m_dot = rfg_m_dot
        self.core.ext.fluid_mean = air_mean
        self.core.int.fluid_mean = rfg_mean

        # Determine the heat transfer rate through the heat exchanger core:
        cnt_flow_hex = dry.CounterFlowHeatExchanger(
            C_cold=rfg_mean.cp * rfg_m_dot,
            C_hot=air_mean.cp * self.air_m_dot,
            T_cold_in=self.rfg_in.T,
            T_hot_in=self.air_in.Tdb,
            UA=self.core.UA
        )
        dev = (cnt_flow_hex.Q - Q_dot).to('kW')

        i = counter[0]
        logger.debug(
            f"Superheating region/Iteration {i + 1}: "
            f"Refrigerant mass flow rate = {rfg_m_dot.to('kg / hr'):~P.3f}. "
            f"Try flow length {self.core.L2.to('mm'):~P.3f}. "
            f"Heat transfer deviation = {dev:~P.3f}."
        )
        counter[0] += 1
        return dev.m

    def solve(
        self,
        rfg_m_dot: Quantity,
        L_flow_max: Quantity
    ) -> tuple[Quantity, HumidAir]:
        """
        Finds the flow length of the superheating region for a given mass flow
        rate of refrigerant so that the refrigerant leaving the evaporator has
        the degree of superheating set on the expansion device.

        Parameters
        ----------
        rfg_m_dot:
            Mass flow rate of refrigerant.
        L_flow_max:
            The maximum flow length possible of the superheating region.

        Returns
        -------
        The flow length of the evaporator's superheating region and the state
        of air leaving the superheating region.
        """
        # Heat absorbed by the refrigerant stream:
        Q_dot = rfg_m_dot * (self.rfg_out.h - self.rfg_in.h)

        # As the heat absorbed by the refrigerant stream is also the heat
        # rejected by the air stream, the state of air leaving the superheating
        # region can be determined:
        T_air_out = self.air_in.Tdb - Q_dot / (CP_HUMID_AIR * self.air_m_dot)
        # Assume only sensible air cooling in the superheating region:
        self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)

        # The mean states of air and refrigerant in the superheating region:
        air_mean, rfg_mean = self._get_fluid_mean_states(self.air_out, rfg_m_dot)

        # Find the flow length of the superheating region for which the transfer
        # of heat through the heat exchanger balances the heat absorbed by the
        # refrigerant stream and the heat rejected by the air stream:
        counter = [0]
        try:
            sol = optimize.root_scalar(
                self.__fun__,
                args=(rfg_m_dot, Q_dot, air_mean, rfg_mean, counter),
                bracket=(Q_(1.0, 'mm').m, L_flow_max.to('mm').m),
                xtol=0.01,  # mm
                rtol=0.001,  # 0.1 %
                maxiter=20
            )
        except ValueError:
            raise SuperheatingError(
                "The required degree of refrigerant superheating "
                "cannot be reached under the current operating "
                "conditions of the evaporator."
            ) from None
        else:
            self.L_flow = Q_(sol.root, 'mm')
            return self.L_flow, self.air_out

    def _get_fluid_mean_states(
        self,
        air_out: HumidAir,
        rfg_m_dot: FluidState
    ) -> tuple[HumidAir, FluidState]:
        """
        Determines the mean state of air and refrigerant in the superheating
        region of the evaporator.

        We need this to calculate the overall heat transfer conductance of the
        superheating region in method `__fun__`.
        """
        T_air_avg = (self.air_in.Tdb + air_out.Tdb) / 2
        W_air_avg = (self.air_in.W + air_out.W) / 2
        air_avg = HumidAir(Tdb=T_air_avg, W=W_air_avg)

        Rfg = self.rfg_in.fluid
        P_evp = self.rfg_in.P
        T_rfg_avg = (self.rfg_in.T + self.rfg_out.T) / 2
        rfg_avg = Rfg(P=P_evp, T=T_rfg_avg)

        C_air = self.air_m_dot * air_avg.cp
        C_rfg = rfg_m_dot * rfg_avg.cp
        C_max = max(C_air, C_rfg)
        C_min = min(C_air, C_rfg)
        C_rat = C_min / C_max

        if C_rat >= 0.5:
            return air_avg, rfg_avg
        else:
            lmtd = self._get_lmtd(air_out)
            if C_max == C_rfg:
                T_air_avg = rfg_avg.T + lmtd
                air_avg = HumidAir(Tdb=T_air_avg, W=W_air_avg)
                return air_avg, rfg_avg
            else:
                T_rfg_avg = air_avg.Tdb - lmtd
                rfg_avg = Rfg(T=T_rfg_avg, P=P_evp)
                return air_avg, rfg_avg
    
    def _get_lmtd(self, air_out: HumidAir) -> Quantity:
        """
        Calculates the LMTD in the superheating region.

        This is just a helper method inside method `_get_fluid_mean_states`.
        """
        dT_in = self.air_in.Tdb - self.rfg_out.T
        dT_out = air_out.Tdb - self.rfg_in.T
        dT_max = max(dT_in, dT_out)
        dT_min = min(dT_in, dT_out)
        if dT_min.m <= 0.0:
            logger.debug(
                f"Calculated LMTD. `dT_min` was {dT_min.to('K'):~P.3g}, but cannot "
                f"be zero or negative. It has been changed to a positive value near "
                f"zero."
            )
            dT_min = Q_(1.e-12, 'K')
        lmtd = (dT_max - dT_min) / np.log(dT_max / dT_min)
        return lmtd

    @property
    def Q_dot(self) -> Quantity:
        """
        Returns the heat transfer rate in the superheating region of the
        evaporator.
        """
        Q_dot = self.rfg_m_dot * (self.rfg_out.h - self.rfg_in.h)
        return Q_dot

    @property
    def Q_dot_max(self) -> Quantity:
        """
        Returns the theoretically maximum heat transfer rate between the
        air and refrigerant stream in the superheating region.
        """
        air_mean, rfg_mean = self._get_fluid_mean_states(self.air_out, self.rfg_m_dot)
        C_air = self.air_m_dot * air_mean.cp
        C_rfg = self.rfg_m_dot * rfg_mean.cp
        C_min = min(C_air, C_rfg)
        Q_dot_max = C_min * (self.air_in.Tdb - self.rfg_in.T)
        return Q_dot_max

    @property
    def dP_air(self) -> Quantity:
        """
        Returns the air-side pressure drop along the superheating region of the
        evaporator.
        """
        dP_air = self.core.ext.get_pressure_drop(self.air_in, self.air_out)
        return dP_air


class BoilingRegion:

    def __init__(
        self,
        W_fro: Quantity,
        H_fro: Quantity,
        S_trv: Quantity,
        S_lon: Quantity,
        D_int: Quantity,
        D_ext: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity
    ) -> None:
        """
        Creates the boiling region of the plain fin-tube evaporator.

        Parameters
        ----------
        W_fro:
            Width of the frontal area of the evaporator.
        H_fro:
            Height of the frontal area.
        S_trv:
            Transversal pitch, i.e., the spacing between tubes in a single row.
        S_lon:
            Longitudinal pitch, i.e., the spacing between tubes of adjacent rows.
        D_int:
            Inside diameter of the tubes.
        D_ext:
            Outside diameter of the tubes.
        t_fin:
            Thickness of the plain fins.
        N_fin:
            Number of fins per unit length of tube (i.e., the inverse of fin
            density).
        k_fin:
            Thermal conductivity of the fin material.
        """
        self.core = PlainFinTubeHeatExchangerCore(
            L1=W_fro,
            L3=H_fro,
            S_t=S_trv,
            S_l=S_lon,
            D_i=D_int,
            D_o=D_ext,
            t_f=t_fin,
            N_f=N_fin,
            k_fin=k_fin,
            boiling=True
        )
        # Known parameters:
        self.air_m_dot: Quantity | None = None
        self.rfg_in: FluidState | None = None
        self.rfg_out: FluidState | None = None

        # Unknown parameters to be solved for:
        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.rfg_m_dot: Quantity | None = None
        self.L_flow: Quantity | None = None

    def solve(
        self,
        air_in: HumidAir,
        rfg_m_dot: Quantity,
        L_flow: Quantity,
        tol: Quantity = Q_(0.1, 'kg / hr'),
        i_max: int = 10
    ) -> Quantity:
        """
        Solves for the refrigerant mass flow rate that is needed to boil the
        entering liquid/vapor mixture completely into saturated vapor, given
        the available flow length `L_flow` of the boiling region.

        Parameters
        ----------
        air_in:
            State of air entering the boiling region (this is also the state
            of air leaving the superheating region).
        rfg_m_dot:
            Initial guess for the refrigerant mass flow rate through the
            evaporator.
        L_flow:
            The available flow length of the boiling region (after the flow
            length of the superheating region has been determined).
        tol:
            Allowable difference between the last and previous calculated
            refrigerant mass flow rate for which iteration is terminated.
        i_max:
            The maximum number of iterations to find an acceptable solution
            for the refrigerant mass flow rate.

        Returns
        -------
        The calculated mass flow rate of refrigerant.
        """
        # Get the heat flow rate absorbed by the refrigerant stream:
        Q_dot = rfg_m_dot * (self.rfg_out.h - self.rfg_in.h)

        # Initially, assume saturated air leaving the evaporator:
        self.air_in = air_in
        h_air_out = air_in.h - Q_dot / self.air_m_dot
        air_out = HumidAir(h=h_air_out, RH=Q_(100, 'pct'))

        # Determine by iteration the refrigerant mass flow rate that is needed
        # to balance the heat transfer rate through the heat exchanger with the
        # heat the refrigerant stream must absorb to become a saturated
        # vapor in the boiling region of the evaporator:
        rfg_mean = self._get_rfg_mean()
        for i in range(i_max):
            air_mean = self._get_air_mean(air_in, air_out)

            # Set the parameters on the heat exchanger core needed to determine
            # the air-side and refrigerant-side heat transfer coefficients, and
            # the air-side finned surface efficiency in the superheated region:
            self.core.L2 = L_flow
            self.core.ext.m_dot = self.air_m_dot
            self.core.int.fluid_mean = rfg_mean
            self.core.int.m_dot = rfg_m_dot
            self.core.ext.fluid_mean = air_mean
            self.core.int.Q = Q_dot

            # Determine the heat transfer rate through heat exchanger core:
            cnt_flow_hex = wet.CounterFlowHeatExchanger(
                m_dot_r=rfg_m_dot,
                m_dot_a=self.air_m_dot,
                T_r_in=self.rfg_in.T,
                T_r_out=self.rfg_out.T,
                P_r=self.rfg_in.P,
                refrigerant=self.rfg_in.fluid,
                air_in=air_in,
                h_ext=self.core.ext.h,
                h_int=self.core.int.h,
                eta_surf_wet=self.core.ext.eta,
                A_ext_to_A_int=self.core.ext.geo.A / self.core.int.geo.A,
                A_ext=self.core.ext.geo.A
            )
            # Note: It is assumed that in the boiling region the air-side heat
            # transfer surface is fully wet.

            # Determine a new value for the refrigerant mass flow rate:
            rfg_m_dot_new = cnt_flow_hex.Q / (self.rfg_out.h - self.rfg_in.h)

            # Check the deviation between the current and previous calculated
            # mass flow rate:
            dev = (rfg_m_dot_new - rfg_m_dot).to('kg / hr')

            logger.debug(
                f"Boiling region/Iteration {i + 1}: "
                f"Flow length = {L_flow.to('mm'):~P.3f}. "
                f"Mass flow rate deviation {dev:~P.3f}."
            )

            if abs(dev) < tol.to('kg / hr'):
                self.rfg_m_dot = rfg_m_dot
                self.air_out = air_out
                self.L_flow = L_flow
                return rfg_m_dot

            # No final solution: Assign the new values and repeat the loop
            # calculations:
            rfg_m_dot = rfg_m_dot_new
            Q_dot = cnt_flow_hex.Q
            air_out = cnt_flow_hex.air_out
        else:
            raise BoilingError(
                f"No acceptable solution found after {i_max} iterations."
            ) from None

    def _get_air_mean(self, air_in: HumidAir, air_out: HumidAir) -> HumidAir:
        """
        Calculates the mean state of air in the boiling region.
        We need this to calculate the heat transfer coefficients in method
        `solve`.
        """
        # Calculate LMED to determine average air enthalpy:
        sat_air_in = HumidAir(Tdb=self.rfg_out.T, RH=Q_(100, 'pct'))
        sat_air_out = HumidAir(Tdb=self.rfg_in.T, RH=Q_(100, 'pct'))
        dh_in = air_in.h - sat_air_in.h
        dh_out = air_out.h - sat_air_out.h
        dh_max = max(dh_in, dh_out)
        dh_min = min(dh_in, dh_out)
        if dh_min.m <= 0.0:
            logger.debug(
                f"Calculated LMED. `dh_min` was {dh_min.to('kJ / kg'):~P.3g}, "
                f"but cannot be zero or negative. It has been changed to a "
                f"positive value near zero."
            )
            dh_min = Q_(1.e-12, 'kJ / kg')
        lmed = (dh_max - dh_min) / np.log(dh_max / dh_min)

        h_sat_air_avg = (sat_air_in.h + sat_air_out.h) / 2
        h_air_avg = h_sat_air_avg + lmed

        # Calculate LMTD to determine average air temperature:
        dT_in = air_in.Tdb - self.rfg_out.T
        dT_out = air_out.Tdb - self.rfg_in.T
        dT_max = max(dT_in, dT_out)
        dT_min = min(dT_in, dT_out)
        if dT_min.m <= 0.0:
            logger.debug(
                f"Calculated LMTD. `dT_min` was {dT_min.to('K'):~P.3g}, but "
                f"cannot be zero or negative. It has been changed to a positive "
                f"value near zero."
            )
            dT_min = Q_(1.e-12, 'K')
        lmtd = (dT_max - dT_min) / np.log(dT_max / dT_min)

        T_rfg_avg = (self.rfg_in.T + self.rfg_out.T) / 2
        T_air_avg = T_rfg_avg + lmtd

        # Mean temperature and mean enthalpy determine the mean air state:
        air_mean = HumidAir(Tdb=T_air_avg, h=h_air_avg)
        return air_mean

    def _get_rfg_mean(self) -> FluidState:
        """
        Calculates the mean state of refrigerant in the boiling region.
        We need this to calculate the heat transfer coefficients in method
        `solve`.
        """
        Rfg = self.rfg_in.fluid

        # Determine average enthalpy of refrigerant:
        h_rfg_avg = (self.rfg_in.h + self.rfg_out.h) / 2

        # Evaporation pressure is assumed constant in evaporator:
        P_evp = self.rfg_in.P
        rfg_mean = Rfg(P=P_evp, h=h_rfg_avg)
        return rfg_mean

    @property
    def Q_dot(self) -> Quantity:
        """
        Returns the heat transfer rate in the boiling region of the evaporator.
        """
        Q_dot = self.rfg_m_dot * (self.rfg_out.h - self.rfg_in.h)
        return Q_dot

    @property
    def Q_dot_max(self) -> Quantity:
        """
        Returns the theoretically maximum heat transfer rate between the
        air and refrigerant stream in the boiling region.
        """
        h_air_in = self.air_in.h
        h_sat_air_out = HumidAir(Tdb=self.rfg_in.T, RH=Q_(100, 'pct')).h
        Q_max = self.air_m_dot * (h_air_in - h_sat_air_out)
        return Q_max

    @property
    def dP_air(self) -> Quantity:
        """
        Returns the air-side pressure drop along the boiling region of the
        evaporator.
        """
        dP_air = self.core.ext.get_pressure_drop(self.air_in, self.air_out)
        return dP_air


class PlainFinTubeCounterFlowAirEvaporator:
    """
    Model class for a plain fin-tube counter-flow air evaporator (DX coil).

    Attributes
    ----------
    air_in: HumidAir:
        State of air entering the evaporator (known).
    air_m_dot: Quantity
        Mass flow rate of air through the evaporator (known).
    rfg_in: FluidState
        State of refrigerant entering the evaporator (known).
    dT_sh: Quantity:
        Required degree of refrigerant superheating at the evaporator outlet
        (known).
    rfg_sat_vap: FluidState
        State of refrigerant as a saturated vapor.
    rfg_out: FluidState
        State of refrigerant leaving the evaporator.
    P_evp: Quantity
        Evaporation pressure.
    T_evp: Quantity
        Evaporation temperature.
    rfg_m_dot: Quantity
        Mass flow rate of refrigerant (calculated after calling method `solve`).
    air_out: HumidAir
        State of air leaving the evaporator (calculated after calling method
        `solve`).
    Q_dot: Quantity
        Refrigeration capacity of the evaporator under the current operating
        conditions (calculated after calling method `solve`).
    Q_dot_max: Quantity
        Theoretical maximum refrigeration capacity of the evaporator under the
        current operating conditions (calculated after calling method `solve`).
    eps: Quantity
        Heat transfer effectiveness of the evaporator under the current
        operating conditions (calculated after calling method `solve`).
    air_dP: Quantity
        Air-side pressure drop under the current operating conditions
        (calculated after calling method `solve`).
    """
    
    def __init__(
        self,
        W_fro: Quantity,
        H_fro: Quantity,
        N_rows: int,
        S_trv: Quantity,
        S_lon: Quantity,
        D_int: Quantity,
        D_ext: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)')
    ) -> None:
        """
        Creates the plain fin-tube counter-flow air evaporator.

        Parameters
        ----------
        W_fro:
            Width of the frontal area of the evaporator.
        H_fro:
            Height of the frontal area.
        N_rows:
            The number of rows in the evaporator.
        S_trv:
            Transversal pitch, i.e., the spacing between tubes in a single row.
        S_lon:
            Longitudinal pitch, i.e., the spacing between tubes of adjacent rows.
        D_int:
            Inside diameter of the tubes.
        D_ext:
            Outside diameter of the tubes.
        t_fin:
            Thickness of the plain fins.
        N_fin:
            Number of fins per unit length of tube (i.e., the inverse of fin
            density).
        k_fin:
            Thermal conductivity of the fin material. The default value applies
            to aluminum.
        """
        self.superheating_region = SuperheatingRegion(
            W_fro, H_fro, S_trv, S_lon, D_int,
            D_ext, t_fin, N_fin, k_fin
        )
        self.boiling_region = BoilingRegion(
            W_fro, H_fro, S_trv, S_lon, D_int,
            D_ext, t_fin, N_fin, k_fin
        )
        self.L_flow = N_rows * S_lon

        # Known parameters:
        self.air_in: HumidAir | None = None
        self.air_m_dot: Quantity | None = None
        self.rfg_in: FluidState | None = None
        self.dT_sh: Quantity | None = None
        self.rfg_sat_vap: FluidState | None = None
        self.rfg_out: FluidState | None = None
        self.P_evp: Quantity | None = None
        self.T_evp: Quantity | None = None

        # Unknown parameters to be solved for:
        self.rfg_m_dot: Quantity | None = None
        self.air_out: HumidAir | None = None

        # Parameters that can be determined when the unknown parameters are
        # solved:
        self.Q_dot: Quantity | None = None
        self.Q_dot_max: Quantity | None = None
        self.eps: Quantity | None = None
        self.air_dP: Quantity | None = None

    def _guess_rfg_m_dot_max(self) -> Quantity:
        """
        Returns a guess for the maximum refrigerant mass flow rate, assuming
        that the air is cooled to the temperature of the refrigerant entering the
        evaporator and is fully saturated.
        """
        sat_air_out = HumidAir(Tdb=self.rfg_in.T, RH=Q_(100, 'pct'))
        Q_dot_max = self.air_m_dot * (self.air_in.h - sat_air_out.h)
        rfg_m_dot_max = Q_dot_max / (self.rfg_out.h - self.rfg_in.h)
        return rfg_m_dot_max

    def __fun__(self, rfg_m_dot: float, counter: list[int]) -> float:
        """
        Returns the deviation between the mass flow rate in the boiling region
        to completely boil the refrigerant to saturated vapor and the mass
        flow rate in the superheating region to superheat the refrigerant to
        the degree set on the expansion device.
        """
        rfg_m_dot = Q_(rfg_m_dot, 'kg / hr')
        i = counter[0]

        # Superheating region: determine the superheating flow length for
        # the given refrigerant mass flow rate.
        L_flow_superheat, air_out_shr = self.superheating_region.solve(
            rfg_m_dot=rfg_m_dot,
            L_flow_max=self.L_flow
        )

        logger.debug(
            f"Superheating region/Iteration {i + 1}: "
            f"Refrigerant mass flow rate = {rfg_m_dot.to('kg / hr'):~P.3f}. "
            f"Superheating flow length = {L_flow_superheat.to('mm'):~P.3f}. "
            f"Determine mass flow rate in boiling region..."
        )

        # Boiling region: get the mass flow rate of refrigerant needed
        # to boil the refrigerant completely within the available boiling flow
        # length:
        L_flow_boil = self.L_flow - L_flow_superheat
        rfg_m_dot_new = self.boiling_region.solve(
            air_in=air_out_shr,
            rfg_m_dot=rfg_m_dot,
            L_flow=L_flow_boil,
        )

        # Determine the deviation between the mass flow rate in the boiling
        # region and the mass flow rate in the superheating region. Ultimately,
        # the deviation should become zero.
        dev = (rfg_m_dot_new - rfg_m_dot).to('kg / hr')

        logger.debug(
            f"Boiling region/Iteration {i + 1}: "
            f"Boiling flow length = {L_flow_boil.to('mm'):~P.3f}. "
            "Deviation between refrigerant mass flow rate in boiling and "
            f"superheating region = {dev:~P.3f}."
        )

        counter[0] += 1
        return dev.m

    def solve(
        self,
        air_in: HumidAir,
        air_m_dot: Quantity,
        rfg_in: FluidState,
        dT_sh: Quantity,
        rfg_m_dot_ini: Quantity | None = None
    ) -> Quantity:
        """
        Solves for the refrigerant mass flow rate that is needed to transform
        the liquid/vapor mixture entering the evaporator into superheated vapor
        at the evaporator's outlet under the given operating conditions.

        Parameters
        ----------
        air_in:
            State of the air entering the evaporator.
        air_m_dot:
            Mass flow rate of air through the evaporator.
        rfg_in:
            State of the refrigerant entering the evaporator.
        dT_sh:
            The required degree of refrigerant superheating set on the expansion
            device.
        rfg_m_dot_ini:
            Optional, initial guess for the refrigerant mass flow rate through
            the evaporator.

        Returns
        -------
        The calculated mass flow rate of refrigerant.

        Notes
        -----
        The evaporator has two regions. First, the refrigerant is boiled into
        saturated vapor. Next, the refrigerant is superheated to the degree
        being set on the expansion device.
        Under steady-state operation, the mass flow rate of refrigerant through
        the evaporator is such that the degree of refrigerant superheating is
        maintained to the setting on the expansion device.
        The flow length needed to superheat the refrigerant depends on this mass
        flow rate.
        We make an initial guess for the refrigerant mass flow rate. With this,
        we determine the superheating flow length so that the heat transfer
        rate from air to refrigerant through the superheating part of the
        heat exchanger balances with the heat rate needed to superheat the
        refrigerant vapor.
        Since the total flow length of the evaporator is fixed, we also retrieve
        the flow length that remains available for refrigerant boiling. Now, we
        calculate the mass flow rate of refrigerant required to boil the
        refrigerant from the known inlet state into fully saturated vapor.
        Using a root-finding algorithm, these steps are repeated until the
        deviation between the mass flow rate in the boiling region and in the
        superheating region has become small enough.
        """
        # Assign or calculate what is known or can be known directly:
        self.air_in = air_in
        self.superheating_region.air_in = air_in
        self.air_m_dot = air_m_dot
        self.superheating_region.air_m_dot = air_m_dot
        self.boiling_region.air_m_dot = air_m_dot
        self.rfg_in = rfg_in
        self.boiling_region.rfg_in = rfg_in
        self.dT_sh = dT_sh.to('K')
        Rfg = self.rfg_in.fluid
        self.P_evp = self.rfg_in.P
        self.rfg_sat_vap = Rfg(P=self.P_evp, x=Q_(1, 'frac'))
        self.boiling_region.rfg_out = self.rfg_sat_vap
        self.superheating_region.rfg_in = self.rfg_sat_vap
        self.T_evp = self.rfg_sat_vap.T
        T_rfg_out = self.T_evp + self.dT_sh
        self.rfg_out = Rfg(P=self.P_evp, T=T_rfg_out)
        self.superheating_region.rfg_out = self.rfg_out

        # Solve for the unknown mass flow rate of refrigerant through the
        # evaporator so that the required degree of refrigerant superheating
        # is attained under the current operating conditions:
        logger.debug(
            "Solving for the operating state of the evaporator..."
        )
        counter = [0]
        if rfg_m_dot_ini is not None:
            rfg_m_dot_max = max(
                self._guess_rfg_m_dot_max().to('kg / hr').m,
                rfg_m_dot_ini.to('kg / hr').m
            )
        else:
            rfg_m_dot_max = self._guess_rfg_m_dot_max().to('kg / hr').m
        try:
            sol = optimize.root_scalar(
                self.__fun__,
                args=(counter,),
                method='secant',
                x0=0.5 * rfg_m_dot_max,
                x1=rfg_m_dot_max,
                xtol=0.01,   # kg/hr
                rtol=0.001,  # 0.1 %
                maxiter=20
            )
        except ValueError:
            message = (
                "The refrigerant cannot be superheated to the required degree "
                "in the evaporator under the current operating conditions."
            )
            logger.error(message)
            raise EvaporatorError(message) from None
        except (BoilingError, SuperheatingError) as err:
            logger.error(f"{type(err).__name__}: {err}")
            raise err
        else:
            logger.debug(
                f"Calculation finished: {sol.flag}"
            )
            self.rfg_m_dot = Q_(sol.root, 'kg / hr')
            self.boiling_region.rfg_m_dot = self.rfg_m_dot
            self.superheating_region.rfg_m_dot = self.rfg_m_dot
            self.Q_dot = self.boiling_region.Q_dot + self.superheating_region.Q_dot
            self.air_out = self.boiling_region.air_out
            self.Q_dot_max = (
                self.boiling_region.Q_dot_max
                + self.superheating_region.Q_dot_max
            )
            self.eps = self.Q_dot / self.Q_dot_max
            self.air_dP = (
                self.boiling_region.dP_air
                + self.superheating_region.dP_air
            )
            return self.rfg_m_dot
