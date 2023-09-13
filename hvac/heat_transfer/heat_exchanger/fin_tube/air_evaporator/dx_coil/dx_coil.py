import numpy as np
from scipy import optimize
from hvac import Quantity
from hvac.fluids import Fluid, HumidAir, FluidState
from hvac.heat_transfer.heat_exchanger.fin_tube.core import PlainFinTubeHeatExchangerCore
import hvac.heat_transfer.heat_exchanger.eps_ntu as eps_ntu_dry
import hvac.heat_transfer.heat_exchanger.eps_ntu_wet as eps_ntu_wet
from hvac.logging import ModuleLogger
from hvac.air_conditioning import AirConditioningProcess


Q_ = Quantity
logger = ModuleLogger.get_logger(__name__)
logger.setLevel(ModuleLogger.ERROR)


class SuperheatedRegion:

    def __init__(
        self,
        Rfg: Fluid,
        W_fro: Quantity,
        H_fro: Quantity,
        S_trv: Quantity,
        S_lon: Quantity,
        D_in: Quantity,
        D_out: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
    ) -> None:
        """Creates the superheating region of the DX-coil."""
        self.Rfg = Rfg
        self.core = PlainFinTubeHeatExchangerCore(
            W_fro, H_fro,
            S_trv, S_lon,
            D_in, D_out,
            t_fin, N_fin, k_fin
        )

        self.air_in: HumidAir | None = None
        self.m_dot_air: Quantity | None = None
        self.T_evp: Quantity | None = None
        self.dT_sh: Quantity | None = None

        self.P_evp: Quantity | None = None
        self.rfg_sat_vap: FluidState | None = None
        self.rfg_out: FluidState | None = None

        self._m_dot_rfg: Quantity | None = None
        self.Q_dot: Quantity | None = None
        self.air_out: HumidAir | None = None

    def set_fixed_operating_conditions(
        self,
        air_in: HumidAir,
        m_dot_air: Quantity,
        T_evp: Quantity,
        dT_sh: Quantity
    ) -> None:
        """Sets the known operating conditions on the superheating region of the
        DX-coil.

        Parameters
        ----------
        air_in:
            State of the air entering the DX-coil.
        m_dot_air:
            The mass flow rate of air through the DX-coil.
        T_evp:
            The evaporation temperature of the refrigerant in the DX-coil.
        dT_sh:
            The degree of refrigerant superheating, which is a setting of the
            expansion device.
        """
        self.air_in = air_in
        self.m_dot_air = m_dot_air
        self.T_evp = T_evp.to('K')
        self.dT_sh = dT_sh.to('K')

        self.rfg_sat_vap = self.Rfg(T=self.T_evp, x=Q_(1, 'frac'))
        self.P_evp = self.rfg_sat_vap.P
        self.rfg_out = self.Rfg(P=self.P_evp, T=self.T_evp + self.dT_sh)

        self._m_dot_rfg = None
        self.Q_dot = None
        self.air_out = None

    @property
    def m_dot_rfg(self) -> Quantity:
        """Returns the mass flow rate of refrigerant through the DX-coil."""
        return self._m_dot_rfg

    @m_dot_rfg.setter
    def m_dot_rfg(self, v: Quantity) -> None:
        """Sets the mass flow rate of refrigerant through the DX-coil.

        Once the mass flow rate of refrigerant is set, the heat absorption rate
        of refrigerant in the superheated region of the DX-coil can be determined.
        It is assumed that the heat transfer surface on the air side is dry.
        So, the humidity ratio of air leaving the superheated region is assumed
        to be equal to the humidity ratio of the entering air.
        """
        self._m_dot_rfg = v
        self.Q_dot = self._m_dot_rfg * (self.rfg_out.h - self.rfg_sat_vap.h)
        h_air_out = self.air_in.h - self.Q_dot / self.m_dot_air
        W_air_out = self.air_in.W
        self.air_out = HumidAir(h=h_air_out, W=W_air_out)

    def _get_mean_fluid_states(self) -> tuple[HumidAir, FluidState]:
        """Calculates the mean state of air and the mean state of refrigerant
        along the superheating region of the DX-coil. Only to be used inside
        method `L_flow`.
        """
        T_air_avg = (self.air_in.Tdb.to('K') + self.air_out.Tdb.to('K')) / 2
        W_air_avg = self.air_in.W
        air_avg = HumidAir(Tdb=T_air_avg, W=W_air_avg)
        cp_air = air_avg.cp
        C_air = self.m_dot_air * cp_air

        T_rfg_avg = (self.rfg_sat_vap.T.to('K') + self.rfg_out.T.to('K')) / 2
        P_rfg_avg = self.P_evp
        rfg_avg = self.Rfg(T=T_rfg_avg, P=P_rfg_avg)
        cp_rfg = rfg_avg.cp
        C_rfg = self.m_dot_rfg * cp_rfg

        C_max = max(C_air, C_rfg)
        C_min = min(C_air, C_rfg)
        C_r = C_min / C_max

        if C_r >= 0.5:
            return air_avg, rfg_avg
        else:
            DT = (
                self.air_in.Tdb.to('K') - self.rfg_out.T.to('K'),
                self.air_out.Tdb.to('K') - self.rfg_sat_vap.T.to('K')
            )
            DT_max = max(DT)
            DT_min = min(DT)
            if DT_min <= 0:
                logger.debug('DT_min <= 0 -> set to 1e-12 K')
                DT_min = Q_(1e-12, 'K')
            LMTD = (DT_max - DT_min) / np.log(DT_max / DT_min)
            if C_max == C_rfg:
                rfg_mean = rfg_avg
                air_mean = HumidAir(
                    Tdb=rfg_mean.T + LMTD,
                    W=self.air_in.W
                )
            else:
                air_mean = air_avg
                rfg_mean = self.Rfg(
                    T=air_mean.Tdb - LMTD,
                    P=self.rfg_sat_vap.P
                )
            return air_mean, rfg_mean

    def __dev_Q_dot__(
        self,
        L_flow: float,
        air_mean: HumidAir,
        rfg_mean: HumidAir
    ) -> float:
        """Calculates the deviation between the heat transfer rate through the
        heat exchanger for a given flow length `L_flow` and the heat absorption
        rate of the refrigerant. Only to be used inside method `L_flow`.
        """
        self.core.L2 = Q_(L_flow, 'm')
        self.core.m_dot_ext = self.m_dot_air
        self.core.m_dot_int = self.m_dot_rfg
        self.core.ext.fluid_mean = air_mean
        self.core.int.fluid_mean = rfg_mean

        cnt_hex = eps_ntu_dry.CounterFlowHeatExchanger(
            C_cold=rfg_mean.cp * self.m_dot_rfg,
            C_hot=air_mean.cp * self.m_dot_air,
            T_cold_in=self.rfg_sat_vap.T,
            T_hot_in=self.air_in.Tdb,
            UA=self.core.UA
        )

        dev_Q_dot = (cnt_hex.Q - self.Q_dot).to('W').m
        return dev_Q_dot

    def L_flow(self, L_flow_ini: Quantity) -> Quantity:
        """Finds the flow length of the superheated region for which the heat
        transfer rate through the heat exchanger equals the heat absorption
        rate of the refrigerant and the heat extraction rate of the air.
        """
        air_mean, rfg_mean = self._get_mean_fluid_states()
        L_flow_ini = L_flow_ini.to('m').magnitude
        try:
            L_flow = optimize.root_scalar(
                self.__dev_Q_dot__,
                args=(air_mean, rfg_mean),
                bracket=(1.e-6, L_flow_ini)
            ).root
        except ValueError:
            raise ValueError(
                "Unable to reach the required degree of refrigerant "
                "superheating."
            )
        return Q_(L_flow, 'm')


class BoilingRegion:

    def __init__(
        self,
        Rfg: Fluid,
        W_fro: Quantity,
        H_fro: Quantity,
        S_trv: Quantity,
        S_lon: Quantity,
        D_in: Quantity,
        D_out: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)')
    ) -> None:
        """Creates the boiling region of the evaporator."""
        self.Rfg = Rfg
        self.core = PlainFinTubeHeatExchangerCore(
            W_fro, H_fro,
            S_trv, S_lon,
            D_in, D_out,
            t_fin, N_fin, k_fin,
            boiling=True
        )

        self.m_dot_air: Quantity | None = None
        self.T_air_out: Quantity | None = None
        self.T_evp: Quantity | None = None

        self.rfg_sat_vap: FluidState | None = None
        self.P_evp: Quantity | None = None

        self.air_in: HumidAir | None = None
        self.L_flow: Quantity | None = None
        self.W_air_out: Quantity | None = None
        self.m_dot_rfg: Quantity | None = None

    def set_fixed_operating_conditions(
        self,
        m_dot_air: Quantity,
        T_air_out: Quantity,
        T_evp: Quantity,
    ) -> None:
        """Sets the known operating conditions on the boiling part of the
        DX-coil.

        Parameters
        ----------
        m_dot_air:
            The mass flow rate of air through the DX-coil.
        T_air_out:
            The temperature of the air leaving the DX-coil.
        T_evp:
            The evaporation temperature of the refrigerant in the DX-coil.
        """
        self.m_dot_air = m_dot_air
        self.T_air_out = T_air_out
        self.T_evp = T_evp.to('K')

        self.rfg_sat_vap = self.Rfg(
            T=self.T_evp,
            x=Q_(1, 'frac')
        )
        self.P_evp = self.rfg_sat_vap.P

        self.air_in = None
        self.L_flow = None
        self.W_air_out = None
        self.m_dot_rfg = None

    @property
    def Q_dot(self) -> Quantity:
        """Returns the heat rate extracted from the air along the boiling region
        of the evaporator.
        """
        air_out = HumidAir(
            Tdb=self.T_air_out,
            W=self.W_air_out
        )
        Q_dot = self.m_dot_air * (self.air_in.h - air_out.h)
        return Q_dot

    @property
    def rfg_in(self) -> FluidState:
        """Returns the state of refrigerant entering the DX-coil."""
        h_rfg_in = self.rfg_sat_vap.h - self.Q_dot / self.m_dot_rfg
        rfg_in = self.Rfg(P=self.P_evp, h=h_rfg_in)
        return rfg_in

    @property
    def __air_mean__(self) -> HumidAir:
        """Returns the mean state of air along the boiling region of the
        DX-coil. Only to be used inside method `calculate_heat_transfer`.
        """
        air_out = HumidAir(
            Tdb=self.T_air_out,
            W=self.W_air_out
        )
        air_sat = HumidAir(
            Tdb=self.T_evp,
            RH=Q_(100, 'pct')
        )
        dh_in = self.air_in.h - air_sat.h
        dh_out = air_out.h - air_sat.h
        if dh_out.magnitude <= 0.0:
            logger.debug('dh_out <= 0 -> set to 1e-12 K')
            dh_out = Q_(1.e-12, 'J / kg')
        dh_max = max(dh_in, dh_out)
        dh_min = min(dh_in, dh_out)
        lmed = (dh_max - dh_min) / np.log(dh_max / dh_min)
        h_air_avg = air_sat.h + lmed

        dT_in = self.air_in.Tdb - self.T_evp
        dT_out = air_out.Tdb - self.T_evp
        if dT_out.magnitude <= 0.0:
            logger.debug('dT_out <= 0 -> set to 1e-12 K')
            dT_out = Q_(1.e-12, 'K')
        dT_max = max(dT_in, dT_out)
        dT_min = min(dT_in, dT_out)
        lmtd = (dT_max - dT_min) / np.log(dT_max / dT_min)
        T_air_avg = self.T_evp + lmtd

        air_avg = HumidAir(Tdb=T_air_avg, h=h_air_avg)
        return air_avg

    @property
    def __rfg_mean__(self) -> FluidState:
        """Returns the mean state of refrigerant along the boiling region of the
        DX-coil. Only to be used inside method `calculate_heat_transfer`.
        """
        h_rfg_avg = (self.rfg_in.h + self.rfg_sat_vap.h) / 2
        rfg_mean = self.Rfg(
            P=self.P_evp,
            h=h_rfg_avg
        )
        return rfg_mean

    def calculate_heat_transfer(self) -> Quantity:
        """Calculates the heat transfer rate through the heat exchanger in the
        boiling region of the DX-coil.
        """
        self.core.L2 = self.L_flow
        self.core.m_dot_int = self.m_dot_rfg
        self.core.m_dot_ext = self.m_dot_air
        self.core.int.fluid_mean = self.__rfg_mean__
        self.core.ext.fluid_mean = self.__air_mean__
        self.core.int.Q = self.Q_dot

        cnt_hex = eps_ntu_wet.CounterFlowHeatExchanger(
            m_dot_r=self.m_dot_rfg,
            m_dot_a=self.m_dot_air,
            T_r_in=self.rfg_in.T,
            T_r_out=self.rfg_sat_vap.T,
            P_r=self.P_evp,
            refrigerant=self.Rfg,
            air_in=self.air_in,
            h_ext=self.core.ext.h,
            h_int=self.core.int.h,
            eta_surf_wet=self.core.ext.eta,
            A_ext_to_A_int=(
                self.core.ext.geo.alpha
                / self.core.int.geo.alpha
            ).to('m / m').m,
            A_ext=self.core.ext.geo.A
        )
        return cnt_hex.Q


class DXAirCoolingCoil:

    def __init__(
        self,
        Rfg: Fluid,
        W_fro: Quantity,
        H_fro: Quantity,
        N_rows: int,
        S_trv: Quantity,
        S_lon: Quantity,
        D_in: Quantity,
        D_out: Quantity,
        t_fin: Quantity,
        N_fin: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)')
    ) -> None:
        """Creates the model of the DX air-cooling coil.

        Parameters
        ----------
        Rfg:
            The type of refrigerant that flows in the DX coil.
        W_fro:
            The width of the frontal area of the DX-coil.
        H_fro:
            The height of the frontal area.
        N_rows:
            The number of rows of the DX-coil.
        S_trv:
            The transversal pitch, i.e., the vertical distance between adjacent
            tubes in the same row.
        S_lon:
            The longitudinal pitch, i.e., the horizontal distance between tubes
            of adjacent rows.
        D_in:
            Inner diameter of the tubes.
        D_out:
            Outer diameter of the tubes.
        t_fin:
            The thickness of fins.
        N_fin:
            The number of fins per unit length of tube.
        k_fin:
            The thermal conductivity of the fin material. The default value
            applies to aluminum.
        """
        self.Rfg = Rfg
        self.N_rows = N_rows
        self.L_flow = N_rows * S_lon
        self.geometry = {
            'W_fro': W_fro,
            'H_fro': H_fro,
            'S_trv': S_trv,
            'S_lon': S_lon,
            'D_in': D_in,
            'D_out': D_out,
            't_fin': t_fin,
            'N_fin': N_fin,
            'k_fin': k_fin
        }
        self.superheating_region = SuperheatedRegion(Rfg, **self.geometry)
        self.boiling_region = BoilingRegion(Rfg, **self.geometry)

        self.m_dot_air: Quantity | None = None
        self.air_in: HumidAir | None = None
        self.T_air_out: Quantity | None = None
        self.T_evp: Quantity | None = None
        self.dT_sh: Quantity | None = None

        self.P_evp: Quantity | None = None
        self.rfg_out: FluidState | None = None

        self.rfg_in: FluidState | None = None
        self.m_dot_rfg: Quantity | None = None
        self.air_out: HumidAir | None = None
        self.Q_dot: Quantity | None = None

    def set_operating_conditions(
        self,
        m_dot_air: Quantity,
        air_in: HumidAir,
        T_air_out: Quantity,
        T_evp: Quantity,
        dT_sh: Quantity,
    ) -> None:
        """Sets the known operating parameters on the DX-coil.

        Parameters
        ----------
        m_dot_air:
            The mass flow rate of air through the DX-coil.
        air_in:
            The state of the air entering the DX-coil.
        T_air_out:
            The wanted temperature of the air leaving the DX-coil.
        T_evp:
            The evaporation temperature of the refrigerant.
        dT_sh:
            The degree of refrigerant superheating, which is a setting of the
            expansion device.
        """
        self.m_dot_air = m_dot_air
        self.air_in = air_in
        self.T_air_out = T_air_out
        self.T_evp = T_evp.to('K')
        self.P_evp = self.Rfg(T=self.T_evp, x=Q_(1, 'frac')).P
        self.dT_sh = dT_sh.to('K')
        self.rfg_out = self.Rfg(P=self.P_evp, T=self.T_evp + self.dT_sh)

        self.rfg_in = None
        self.m_dot_rfg = None
        self.air_out = None
        self.Q_dot = None

        self.superheating_region.set_fixed_operating_conditions(
            air_in=self.air_in,
            m_dot_air=self.m_dot_air,
            T_evp=self.T_evp,
            dT_sh=self.dT_sh
        )

        self.boiling_region.set_fixed_operating_conditions(
            m_dot_air=self.m_dot_air,
            T_air_out=self.T_air_out,
            T_evp=self.T_evp,
        )

    def solve(
        self,
        x_rfg_in_ini: Quantity = Q_(15, 'pct'),
        tol_Q_dot: Quantity = Q_(1, 'W'),
        i_max: int = 20
    ) -> HumidAir:
        """Solves for the mass flow rate of refrigerant that needs to flow
        through the DX-coil at the given evaporation temperature to cool the
        air to the wanted temperature at the cooling coil outlet.

        Parameters
        ----------
        x_rfg_in_ini:
            The solving method uses iteration and needs an initial guess of the
            vapor quality of the refrigerant entering the DX-coil.
        tol_Q_dot:
            Tolerance on the heat transfer rate. The heat transfer rate through
            the heat exchanger must be equal to the rate that heat is extracted
            from the air and the rate that heat is absorbed by the refrigerant.
            The iterations stop when the deviation between the heat transfer rate
            through the boiling region of the heat exchanger and the rate that
            heat is extracted from the air stream is less than the tolerance.
        i_max:
            The maximum number of iterations. If the deviation is still greater
            than the tolerance after this number of iterations, a ValueError will
            be raised to indicate that no solution within the set tolerance has
            been found.
        """
        # We guess an initial value for the humidity ratio of the air leaving
        # the cooling coil:
        air_adp = HumidAir(Tdb=self.T_evp, RH=Q_(100, 'pct'))
        cooling_coil = AirConditioningProcess(
            air_in=self.air_in,
            m_da=self.m_dot_air,
            T_ao=self.T_air_out,
            ADP=air_adp
        )
        air_out_ini = cooling_coil.air_out
        W_air_out_ini = air_out_ini.W

        # With this, we determine an initial value for the rate of heat extracted
        # from the air stream:
        Q_dot_ini = self.m_dot_air * (self.air_in.h - air_out_ini.h)

        # We also guess an initial value for the vapor quality of the refrigerant
        # entering the cooling coil:
        rfg_in_ini = self.Rfg(P=self.P_evp, x=x_rfg_in_ini)

        # With this, we determine an initial value for the mass flow rate of
        # refrigerant:
        m_dot_rfg_ini = Q_dot_ini / (self.rfg_out.h - rfg_in_ini.h)

        W_air_out = W_air_out_ini
        m_dot_rfg = m_dot_rfg_ini
        for i in range(i_max):

            logger.debug(
                f"Iteration {i}. Try with: "
                f"`W_air_out` = {W_air_out.to('g / kg'):~P.3f}, "
                f"`m_dot_rfg` = {m_dot_rfg.to('kg / hr'):~P.3f}"
            )

            # We determine the flow length of the superheating region with the
            # current value of the refrigerant mass flow rate:
            self.superheating_region.m_dot_rfg = m_dot_rfg
            L_flow_super = self.superheating_region.L_flow(self.L_flow)

            # We set the current parameters of the boiling region:
            self.boiling_region.m_dot_rfg = m_dot_rfg
            self.boiling_region.L_flow = self.L_flow - L_flow_super
            self.boiling_region.air_in = self.superheating_region.air_out
            self.boiling_region.W_air_out = W_air_out

            # We determine the heat transfer rate through the boiling region
            # of the DX-coil:
            Q_dot_boil = self.boiling_region.calculate_heat_transfer()

            # With this, we can determine a new value for the mass flow rate
            # of refrigerant:
            m_dot_rfg_new = (
                Q_dot_boil
                / (self.boiling_region.rfg_sat_vap.h - self.boiling_region.rfg_in.h)
            )

            # Also, we can determine a new value for the humidity ratio of the
            # air leaving the cooling coil:
            h_air_out_new = (
                self.boiling_region.air_in.h
                - Q_dot_boil / self.m_dot_air
            )
            air_out_new = HumidAir(Tdb=self.T_air_out, h=h_air_out_new)
            W_air_out_new = air_out_new.W

            # The heat rate extracted from the air in the boiling part of the
            # DX-coil must be equal to the heat transfer rate through the
            # heat exchanger core to the refrigerant.
            dev = Q_dot_boil.to('kW') - self.boiling_region.Q_dot.to('kW')

            logger.debug(
                f"Deviation on `Q_dot`: {dev:~P.3f}"
            )

            # If the deviation is within the tolerance, stop the iterations and
            # gather the results:
            if abs(dev) < tol_Q_dot.to('kW'):
                self.m_dot_rfg = m_dot_rfg_new
                self.air_out = air_out_new
                self.Q_dot = self.superheating_region.Q_dot + self.boiling_region.Q_dot
                h_rfg_in = self.boiling_region.rfg_sat_vap.h - Q_dot_boil / m_dot_rfg_new
                self.rfg_in = self.Rfg(P=self.P_evp, h=h_rfg_in)
                return air_out_new

            # Else, repeat the calculations with the new values for the mass
            # flow rate of refrigerant and the leaving air humidity ratio.
            m_dot_rfg = m_dot_rfg_new
            W_air_out = W_air_out_new
        else:
            raise ValueError(
                f"No convergence for `W_air_out` after {i_max} iterations."
            )
