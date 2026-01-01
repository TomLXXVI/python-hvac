import math
from abc import ABC, abstractmethod

import numpy as np
from scipy import optimize

from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import FluidState, HumidAir, CP_HUMID_AIR
from hvac.fluids.exceptions import CoolPropMixtureError
import hvac.heat_transfer.forced_convection.internal_flow as single_phase_flow
import hvac.heat_transfer.boiling.flow_boiling as boiling_flow
from ...general import eps_ntu as dry
from ...general import eps_ntu_wet as wet
from .geometry import ContinuousFinStaggeredTubeBank
from .air_to_water import ExternalSurface

Q_ = Quantity

logger = ModuleLogger.get_logger(__name__)
logger.setLevel(ModuleLogger.ERROR)


class EvaporatorError(Exception):
    pass


class SuperheatingError(EvaporatorError):
    pass


class BoilingError(EvaporatorError):
    pass


class HeatExchangerCore:
    """
    Represents the heat exchanger core of a region of the air evaporator
    (boiling region or superheating region). Its task is to calculate the heat
    transfer characteristics under specified operating conditions on the
    refrigerant side and air side of the heat exchanger core. It contains three
    main attributes.

    Attributes
    ----------
    geometry: ContinuousFinStaggeredTubeBank
        Groups all the geometric characteristics of the heat exchanger core.
    internal: BoilingPhaseInternalSurface | SinglePhaseInternalSurface
        Represents the internal, refrigerant-side heat exchanging surface of the
        heat exchanger core.
    external: ExternalSurface
        Represents the external, air-side heat exchanging surface of the heat
        exchanger core.
    """
    def __init__(
        self,
        width: Quantity,
        height: Quantity,
        num_rows: int | None,
        pitch_trv: Quantity,
        pitch_lon: Quantity,
        d_o: Quantity,
        d_i: Quantity,
        t_fin: Quantity,
        fin_density: Quantity,
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
        d_r: Quantity | None = None,
        boiling: bool = False,
        num_circuits: int | None = None
    ) -> None:
        """
        Creates a `HeatExchangerCore` object.

        Parameters
        ----------
        width: Quantity
            Frontal area width.
        height: Quantity
            Frontal area height.
        num_rows: int | None
            Number of rows.
        pitch_trv: Quantity
            Vertical spacing between tubes in one row (transversal pitch).
        pitch_lon: Quantity
            Horizontal spacing between the tubes in two adjacent rows
            (longitudinal pitch).
        d_o: Quantity
            Tube external diameter.
        d_i: Quantity
            Tube internal diameter.
        t_fin: Quantity
            Fin thickness.
        fin_density: Quantity
            Number of fins per unit length of tube (= inverse of fin spacing).
        k_fin: Quantity, default Q_(237, 'W / (m * K)')
            Thermal conductivity of fin material. Default value applies to
            aluminum.
        d_r: Quantity | None = None
            TODO: add description.
        boiling: bool, default False
            If True, indicates refrigerant boiling heat transfer.
        num_circuits: int | None
            Number of refrigerant circuits.
        """
        # Create the geometry of the heat exchanger core.
        self.geometry = ContinuousFinStaggeredTubeBank(
            width, height, num_rows, pitch_trv,
            pitch_lon, d_o, d_i, t_fin, fin_density,
            k_fin, d_r
        )
        self.num_circuits = num_circuits

        # Create the appropriate internal, refrigerant-side heat exchanging
        # surface.
        if boiling:
            self.internal = BoilingPhaseInternalSurface(self)
        else:
            self.internal = SinglePhaseInternalSurface(self)

        # Create the external, air-side heat exchanging surface.
        self.external = ExternalSurface(self)

        self.T_wall: Quantity | None = None
        self.h_int: Quantity | None = None
        self.R_int: Quantity | None = None
        self.h_ext: Quantity | None = None
        self.R_ext: Quantity | None = None
        self.R_tot: Quantity | None = None
        self.UA: Quantity | None = None
        self.eta_surf: float | None = None
        self.dP_ext: Quantity | None = None

    def hex_properties(
        self,
        air_mean: HumidAir,
        air_in: HumidAir,
        air_out: HumidAir,
        rfg_mean: FluidState,
        air_m_dot: Quantity,
        rfg_m_dot: Quantity,
        Q_dot: Quantity
    ) -> dict[str, Quantity | float]:
        """
        Returns the heat exchanger properties under the given operating
        conditions.

        Parameters
        ----------
        air_mean: HumidAir
            Mean state of the air flow along the external heat transfer
            surface.
        air_in: HumidAir
            State of air at the entrance of the external heat transfer surface.
        air_out: HumidAir
            State of air at the exit of the external heat transfer surface.
        rfg_mean: FluidState
            Mean state of the refrigerant flow along the internal heat transfer
            surface.
        air_m_dot: Quantity
            Mass flow rate of air along the external heat transfer surface.
        rfg_m_dot: Quantity
            Mass flow rate of refrigerant along the internal heat transfer
            surface.
        Q_dot: Quantity
            Heat transfer between refrigerant and air.

        Returns
        -------
        dict[str, Quantity | float]
            h_int: Quantity
                Convection heat transfer coefficient on the internal refrigerant
                side.
            h_ext: Quantity
                Convection heat transfer coefficient on the external air side.
            R_int: Quantity
                Thermal resistance on the internal refrigerant side.
            R_ext: Quantity
                Thermal resistance on the external air side.
            R_tot: Quantity
                Total thermal resistance between internal, refrigerant side and
                external, air side.
            UA: Quantity
                Total thermal conductance between internal, refrigerant side and
                external, air side.
            eta_surf:
                Finned surface efficiency of the external, air-side heat transfer
                surface.
            dP_ext:
                Pressure drop of air across the external heat transfer surface.
        """
        self.internal.m_dot = rfg_m_dot
        self.internal.rfg = rfg_mean
        self.external.m_dot = air_m_dot
        self.external.air = air_mean
        self.T_wall = self.__find_wall_temperature(Q_dot)
        self.h_int = self.internal.heat_transfer_coeff(Q_dot)
        self.R_int = self.internal.thermal_resistance(self.h_int)
        self.h_ext = self.external.heat_transfer_coeff(self.T_wall)
        self.R_ext = self.external.thermal_resistance(self.h_ext)
        self.R_tot = self.R_int + self.R_ext
        # Note: The thermal conduction resistance of the heat exchanger body
        # is ignored, and fouling is not being taken into account.
        self.UA = 1 / self.R_tot
        self.dP_ext = self.external.pressure_drop(air_in.rho, air_out.rho, self.T_wall)
        self.eta_surf = self.external.eta_surf(self.h_ext)
        return {
            'h_int': self.h_int,
            'h_ext': self.h_ext,
            'R_int': self.R_int,
            'R_ext': self.R_ext,
            'R_tot': self.R_tot,
            'UA': self.UA,
            'eta_surf': self.eta_surf,
            'dP_ext': self.dP_ext
        }

    def __find_wall_temperature(self, Q_dot: Quantity) -> Quantity:

        def __eq__(T_wall_: float) -> float:
            T_wall_ = Q_(T_wall_, 'K')
            h_int = self.internal.heat_transfer_coeff(Q_dot)
            R_int = self.internal.thermal_resistance(h_int)
            h_ext = self.external.heat_transfer_coeff(T_wall_)
            R_ext = self.external.thermal_resistance(h_ext)
            T_int_ = self.internal.rfg.T.to('K')
            T_ext_ = self.external.air.Tdb.to('K')
            n = T_int_ / R_int + T_ext_ / R_ext
            d = 1 / R_int + 1 / R_ext
            T_wall_new = n / d
            dev = (T_wall_new - T_wall_).to('K').m
            return dev

        T_int = self.internal.rfg.T.to('K').m    # cold fluid
        T_ext = self.external.air.Tdb.to('K').m  # hot fluid
        if T_ext > T_int:
            try:
                sol = optimize.root_scalar(__eq__, bracket=(T_int, T_ext))
            except ValueError as err:
                raise BoilingError(
                    f"Wall temperature cannot be determined: {err}"
                )
            else:
                T_wall = Q_(sol.root, 'K')
        else:
            T_wall = Q_(T_int, 'K')
        return T_wall


class InternalSurface(ABC):
    """
    Represents the internal, refrigerant side of the heat exchanger core.
    """
    def __init__(self, parent: HeatExchangerCore):
        self.parent = parent
        self.rfg: FluidState | None = None
        self.m_dot: Quantity | None = None
        self.tube = None

    def _get_m_dot_tube(self) -> Quantity:
        # It is assumed that every tube in the first row is connected to
        # the supply header.
        if self.parent.num_circuits is None:
            n = self.parent.geometry.num_tubes_1st_row
        else:
            n = self.parent.num_circuits
        m_dot_tube = self.m_dot / n
        return m_dot_tube

    def thermal_resistance(self, h: Quantity) -> Quantity:
        A_tot = self.parent.geometry.internal.A_tot
        R = 1 / (h * A_tot)
        return R

    @abstractmethod
    def heat_transfer_coeff(self, *args, **kwargs) -> Quantity:
        ...


class SinglePhaseInternalSurface(InternalSurface):
    """
    Represents the internal, refrigerant side of the heat exchanger core in the
    superheating region of the evaporator.
    """
    def heat_transfer_coeff(self, *args, **kwargs) -> Quantity:
        self.tube = single_phase_flow.CircularTube(
            Di=self.parent.geometry.d_i,
            L=self.parent.geometry.width,
            fluid=self.rfg
        )
        self.tube.m_dot = self._get_m_dot_tube()
        h = self.tube.avg_heat_transfer_coefficient()
        if isinstance(h, tuple):
            h = h[0]  # if laminar flow, assume constant heat flux
        return h.to('W / (m ** 2 * K)')


class BoilingPhaseInternalSurface(InternalSurface):
    """
    Represents the internal, refrigerant side of the heat exchanger core in the
    boiling region of the evaporator.
    """
    def heat_transfer_coeff(self, Q_dot: Quantity) -> Quantity:
        d_i = self.parent.geometry.d_i
        A_tube = math.pi * d_i ** 2 / 4
        A_tot = self.parent.geometry.internal.A_tot
        self.tube = boiling_flow.Tube(
            Dh=self.parent.geometry.internal.d_h,
            A=A_tube,
            fluid=self.rfg.fluid,
            tube_orient='horizontal'
        )
        self.tube.m_dot = self._get_m_dot_tube()
        h = self.tube.heat_trf_coeff(
            x=self.rfg.x,
            T_sat=self.rfg.T,
            q_s=Q_dot / A_tot
        )
        return h.to('W / (m ** 2 * K)')


class SuperheatingRegion:

    def __init__(
        self,
        geometry: ContinuousFinStaggeredTubeBank,
        num_circuits: int | None = None
    ) -> None:
        """
        Creates the superheating region of the plain fin-tube evaporator.

        Parameters
        ----------
        geometry: ContinuousFinStaggeredTubeBank
            Geometric properties of the heat exchanger.
        num_circuits: int | None
            Number of refrigerant circuits.
        """
        self.core = HeatExchangerCore(
            width=geometry.width,
            height=geometry.height,
            num_rows=None,
            pitch_trv=geometry.pitch_trv,
            pitch_lon=geometry.pitch_lon,
            d_i=geometry.d_i,
            d_o=geometry.d_o,
            t_fin=geometry.t_fin,
            fin_density=geometry.fin_density,
            k_fin=geometry.k_fin,
            num_circuits=num_circuits
        )
        # Known parameters:
        self.air_in: HumidAir | None = None
        self.air_m_dot: Quantity | None = None
        self.rfg_in: FluidState | None = None
        self.rfg_out: FluidState | None = None

        # Parameters to be solved for:
        self.rfg_m_dot: Quantity | None = None
        self.air_out: HumidAir | None = None
        self.Q_dot: Quantity | None = None
        self.L_flow: Quantity | None = None

    def __fun__(
        self,
        L_flow: float,
        counter: list[int]
    ) -> float:
        """
        Calculates the deviation between the heat transfer rate through the
        heat exchanger and the heat absorbed by the refrigerant stream, which
        is also the heat rejected by the air stream.
        """
        # Set the current flow length of the superheating region.
        self.core.geometry.length = Q_(L_flow, 'mm')
        # The mean states of air and refrigerant in the superheating region.
        air_mean, rfg_mean = self._get_fluid_mean_states()
        # Determine the overall heat transfer conductance of the superheated 
        # region.
        hex_props = self.core.hex_properties(
            air_mean=air_mean,
            air_in=self.air_in,
            air_out=self.air_out,
            rfg_mean=rfg_mean,
            air_m_dot=self.air_m_dot,
            rfg_m_dot=self.rfg_m_dot,
            Q_dot=self.Q_dot
        )
        # Determine the heat transfer rate through the heat exchanger core 
        # in the superheated region.
        self.cnt_flow_hex = dry.CounterFlowHeatExchanger(
            C_cold=rfg_mean.cp * self.rfg_m_dot,
            C_hot=air_mean.cp * self.air_m_dot,
            T_cold_in=self.rfg_in.T,
            T_hot_in=self.air_in.Tdb,
            UA=hex_props['UA']
        )
        # Check heat balance. Determine deviation between heat transfer rate 
        # through the heat exchanger core and the heat transfer rate needed to 
        # superheat the refrigerant vapor.
        dev = (self.cnt_flow_hex.Q - self.Q_dot).to('W')
        i = counter[0]
        logger.debug(
            f"{i}. Superheating flow length "
            f"{self.core.geometry.length.to('mm'):~P.3f} -> "
            f"heat transfer rate = {self.cnt_flow_hex.Q.to('W'):~P.3f} "
            f"(deviation = {dev.to('W'):~P.3f})"
        )
        counter[0] += 1
        return dev.m

    def solve(
        self,
        rfg_m_dot: Quantity,
        L_flow_max: Quantity
    ) -> tuple[Quantity, HumidAir]:
        """
        Solves for the required flow length of the superheating region when mass
        flow rate of refrigerant is given so that refrigerant leaving the 
        evaporator has obtained the degree of superheating set on the expansion
        device.

        Parameters
        ----------
        rfg_m_dot:
            Mass flow rate of refrigerant.
        L_flow_max:
            The maximum flow length possible of the superheating region.

        Returns
        -------
        tuple[Quantity, HumidAir]
            The flow length of the evaporator's superheating region and the state
            of air leaving the superheating region.
        """
        self.rfg_m_dot = rfg_m_dot
        # Heat that must be absorbed by the refrigerant stream.
        self.Q_dot = self.rfg_m_dot * (self.rfg_out.h - self.rfg_in.h)
        logger.debug(
            f"Heat transfer required in superheating region = "
            f"{self.Q_dot.to('W'):~P.3f}"
        )
        # As the heat absorbed by the refrigerant is also the heat rejected by 
        # the air, the state of air leaving the superheating region can be 
        # determined. Only sensible air cooling in the superheating region is 
        # assumed.
        T_air_out = self.air_in.Tdb - self.Q_dot / (CP_HUMID_AIR * self.air_m_dot)
        self.air_out = HumidAir(Tdb=T_air_out, W=self.air_in.W)
        # Find the flow length of the superheating region so that the heat 
        # transfer rate from air to refrigerant balances the heat rate that
        # must be absorbed by the refrigerant.
        counter = [0]
        try:
            L_flow_min = Q_(1e-3, 'mm').m
            L_flow_max = L_flow_max.to('mm').m
            sol = optimize.root_scalar(
                self.__fun__,
                args=(counter,),
                bracket=(L_flow_min, L_flow_max),
                xtol=0.01,    # mm
                rtol=0.0001,
                maxiter=20
            )
        except ValueError:
            raise SuperheatingError
        else:
            self.L_flow = Q_(sol.root, 'mm')
            self.__fun__(self.L_flow, counter)
            return self.L_flow, self.air_out

    def _get_fluid_mean_states(self) -> tuple[HumidAir, FluidState]:
        """
        Determines the mean state of air and refrigerant in the superheating
        region of the evaporator.

        We need this to calculate the overall heat transfer conductance of the
        superheating region in method `__fun__`.
        """
        T_air_avg = (self.air_in.Tdb + self.air_out.Tdb) / 2
        W_air_avg = (self.air_in.W + self.air_out.W) / 2
        air_avg = HumidAir(Tdb=T_air_avg, W=W_air_avg)

        Rfg = self.rfg_in.fluid
        P_evp = self.rfg_in.P
        T_rfg_avg = (self.rfg_in.T + self.rfg_out.T) / 2
        rfg_avg = Rfg(P=P_evp, T=T_rfg_avg)

        C_air = self.air_m_dot * air_avg.cp
        C_rfg = self.rfg_m_dot * rfg_avg.cp
        C_max = max(C_air, C_rfg)
        C_min = min(C_air, C_rfg)
        C_rat = C_min / C_max

        if C_rat >= 0.5:
            return air_avg, rfg_avg
        else:
            lmtd = self._get_lmtd(self.air_out)
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
            logger.warning(
                f"Calculated LMTD. `dT_min` was {dT_min.to('K'):~P.3g}, but cannot "
                f"be zero or negative. It has been changed to a positive value near "
                f"zero."
            )
            dT_min = Q_(1.e-12, 'K')
        lmtd = (dT_max - dT_min) / math.log(dT_max / dT_min)
        return lmtd

    @property
    def Q_dot_max(self) -> Quantity:
        """
        Returns the theoretically maximum heat transfer rate between the
        air and refrigerant stream in the superheating region.
        """
        air_mean, rfg_mean = self._get_fluid_mean_states()
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
        dP_air = self.core.external.pressure_drop(
            self.air_in.rho, self.air_out.rho,
            self.core.T_wall
        )
        return dP_air


class BoilingRegion:

    def __init__(
        self,
        geometry: ContinuousFinStaggeredTubeBank,
        num_circuits: int | None = None
    ) -> None:
        """
        Creates the boiling region of the plain fin-tube evaporator.

        Parameters
        ----------
        geometry: ContinuousFinStaggeredTubeBank
            Geometric properties of the heat exchanger.
        num_circuits: int | None
            Number of refrigerant circuits.
        """
        self.core = HeatExchangerCore(
            width=geometry.width,
            height=geometry.height,
            num_rows=None,
            pitch_trv=geometry.pitch_trv,
            pitch_lon=geometry.pitch_lon,
            d_i=geometry.d_i,
            d_o=geometry.d_o,
            t_fin=geometry.t_fin,
            fin_density=geometry.fin_density,
            k_fin=geometry.k_fin,
            boiling=True,
            num_circuits=num_circuits
        )
        # Known parameters:
        self.air_m_dot: Quantity | None = None
        self.rfg_in: FluidState | None = None
        self.rfg_out: FluidState | None = None

        # Unknown parameters to be solved for:
        self.air_in: HumidAir | None = None
        self.air_out: HumidAir | None = None
        self.rfg_m_dot: Quantity | None = None
        self.Q_dot: Quantity | None = None
        self.L_flow: Quantity | None = None
    
    def __fun__(self, L_flow: float, counter: list[int]) -> float:
        # Set the current flow length of the boiling region.
        self.core.geometry.length = Q_(L_flow, 'mm')
        # Mean state of refrigerant in the boiling region.
        rfg_mean = self._get_rfg_mean()
        # Initially assume that leaving air is saturated.
        h_air_out = self.air_in.h - self.Q_dot / self.air_m_dot
        self.air_out = HumidAir(h=h_air_out, RH=Q_(100, 'pct'))
        i = 0
        i_max = 10
        Q_dot = None
        while i < i_max:
            # Mean state of air in the boiling region.
            air_mean = self._get_air_mean()
            # Determine overall heat transfer conductance of the boiling region.
            hex_props = self.core.hex_properties(
                air_mean=air_mean,
                air_in=self.air_in,
                air_out=self.air_out,
                rfg_mean=rfg_mean,
                air_m_dot=self.air_m_dot,
                rfg_m_dot=self.rfg_m_dot,
                Q_dot=self.Q_dot
            )
            # Determine heat transfer rate through heat exchanger core.
            # Note: It is assumed that in the boiling region the air-side heat
            # transfer surface is fully wet.
            A_ext = self.core.geometry.external.A_tot.to('m ** 2')
            A_int = self.core.geometry.internal.A_tot.to('m ** 2')
            self.cnt_flow_hex = wet.CounterFlowHeatExchanger(
                m_dot_r=self.rfg_m_dot,
                m_dot_a=self.air_m_dot,
                T_r_in=self.rfg_in.T,
                T_r_out=self.rfg_out.T,
                T_r_out_ini=None,
                P_r=self.rfg_in.P,
                refrigerant=self.rfg_in.fluid,
                air_in=self.air_in,
                h_ext=hex_props['h_ext'],
                h_int=hex_props['h_int'],
                eta_surf_wet=hex_props['eta_surf'],
                A_ext_to_A_int=A_ext.m / A_int.m,
                A_ext=A_ext
            )
            Q_dot = self.cnt_flow_hex.Q
            air_out = self.cnt_flow_hex.air_out
            dev_air_out = air_out.W.to('g / kg').m - self.air_out.W.to('g / kg').m
            if abs(dev_air_out) < 0.1:
                break
            self.air_out = self.cnt_flow_hex.air_out
            i += 1
        # Check heat balance. Determine deviation between heat transfer rate 
        # through the heat exchanger core and heat transfer rate needed to 
        # boil the refrigerant vapor.
        dev = (Q_dot - self.Q_dot).to('W')
        i = counter[0]

        logger.debug(
            f"{i}. Boiling flow length "
            f"{self.core.geometry.length.to('mm'):~P.3f} -> "
            f"heat transfer rate = {Q_dot.to('W'):~P.3f} "
            f"(deviation = {dev.to('W'):~P.3f})"
        )

        counter[0] += 1
        return dev.m
        
    def solve(
        self,
        air_in: HumidAir,
        rfg_m_dot: Quantity,
        L_flow_max: Quantity,
    ) -> Quantity:
        """
        Solves for the required flow length of the boiling region so that the 
        heat transfer rate through the heat exchanger core balances with the 
        heat rate that must be absorbed by the refrigerant to turn into 
        saturated vapor at the exit of the boiling region.

        Parameters
        ----------
        air_in:
            State of air entering the boiling region (this is also the state
            of air leaving the superheating region).
        rfg_m_dot:
            The refrigerant mass flow rate through the evaporator.
        L_flow_max:
            The maximum flow length possible of the boiling region.

        Returns
        -------
        Quantity
            The required flow length of the boiling region.
        """
        # Set current state of air entering the boiling region and current
        # mass flow rate of refrigerant.
        self.air_in = air_in
        self.rfg_m_dot = rfg_m_dot
        # Determine heat rate that must be absorbed by the refrigerant stream.
        self.Q_dot = self.rfg_m_dot * (self.rfg_out.h - self.rfg_in.h)
        logger.debug(
            f"Heat transfer required in boiling region = "
            f"{self.Q_dot.to('W'):~P.3f}"
        )
        # Find the flow length of the boiling region so that the heat transfer
        # rate through the heat exchanger core balances with the heat rate that
        # must be absorbed by the refrigerant.
        counter = [0]
        L_flow_min = Q_(1.0, 'mm').m
        L_flow_max = L_flow_max.to('mm').m
        try:
            sol = optimize.root_scalar(
                self.__fun__,
                args=(counter,),
                bracket=(L_flow_min, L_flow_max),
                xtol=0.01,   # mm
                rtol=0.001,  # 0.1 %
                maxiter=20
            )
        except ValueError:
            raise BoilingError
        else:
            self.L_flow = Q_(sol.root, 'mm')
            self.__fun__(self.L_flow, counter)
            return self.L_flow

    def _get_air_mean(self) -> HumidAir:
        """
        Calculates the mean state of air in the boiling region.
        We need this to calculate the heat transfer coefficients in method
        `solve`.
        """
        # Calculate LMED to determine average air enthalpy:
        sat_air_in = HumidAir(Tdb=self.rfg_out.T, RH=Q_(100, 'pct'))
        sat_air_out = HumidAir(Tdb=self.rfg_in.T, RH=Q_(100, 'pct'))
        dh_in = self.air_in.h - sat_air_in.h
        dh_out = self.air_out.h - sat_air_out.h
        dh_max = max(dh_in, dh_out)
        dh_min = min(dh_in, dh_out)
        if dh_min.m <= 0.0:
            logger.warning(
                f"Calculated LMED. `dh_min` was {dh_min.to('kJ / kg'):~P.3g}, "
                f"but cannot be zero or negative. It has been changed to a "
                f"positive value near zero."
            )
            dh_min = Q_(1.e-12, 'kJ / kg')
        lmed = (dh_max - dh_min) / math.log(dh_max / dh_min)
        h_sat_air_avg = (sat_air_in.h + sat_air_out.h) / 2
        h_air_avg = h_sat_air_avg + lmed
        # Calculate LMTD to determine average air temperature:
        dT_in = self.air_in.Tdb - self.rfg_out.T
        dT_out = self.air_out.Tdb - self.rfg_in.T
        dT_max = max(dT_in, dT_out)
        dT_min = min(dT_in, dT_out)
        if dT_min.m <= 0.0:
            logger.warning(
                f"Calculated LMTD. `dT_min` was {dT_min.to('K'):~P.3g}, but "
                f"cannot be zero or negative. It has been changed to a positive "
                f"value near zero."
            )
            dT_min = Q_(1.e-12, 'K')
        lmtd = (dT_max - dT_min) / math.log(dT_max / dT_min)

        T_rfg_avg = (self.rfg_in.T + self.rfg_out.T) / 2
        T_air_avg = T_rfg_avg + lmtd

        # Mean temperature and mean enthalpy determine the mean air state:
        air_mean = HumidAir(Tdb=T_air_avg, h=h_air_avg)
        return air_mean

    def _get_rfg_mean(self) -> FluidState:
        """
        Calculates the mean state of refrigerant in the boiling region.
        """
        Rfg = self.rfg_in.fluid
        # Determine average enthalpy of refrigerant.
        h_rfg_avg = (self.rfg_in.h + self.rfg_out.h) / 2
        # Evaporation pressure is assumed constant in evaporator.
        P_evp = self.rfg_in.P
        try:
            rfg_mean = Rfg(P=P_evp, h=h_rfg_avg)
        except CoolPropMixtureError:
            rfg_mean = Rfg(P=P_evp, h=h_rfg_avg, x=Q_(0.0, 'frac'))
        return rfg_mean

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
        dP_air = self.core.external.pressure_drop(
            self.air_in.rho, self.air_out.rho,
            self.core.T_wall
        )
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
        k_fin: Quantity = Q_(237, 'W / (m * K)'),
        num_circuits: int | None = None
    ) -> None:
        """
        Creates the plain fin-tube counter-flow air evaporator.

        Parameters
        ----------
        W_fro: Quantity
            Width of the frontal area of the evaporator.
        H_fro: Quantity
            Height of the frontal area.
        N_rows: Quantity
            The number of rows in the evaporator.
        S_trv: Quantity
            Transversal pitch, i.e., the spacing between tubes in a single row.
        S_lon: Quantity
            Longitudinal pitch, i.e., the spacing between tubes of adjacent rows.
        D_int: Quantity
            Inside diameter of the tubes.
        D_ext: Quantity
            Outside diameter of the tubes.
        t_fin: Quantity
            Thickness of the plain fins.
        N_fin: Quantity
            Number of fins per unit length of tube (i.e., the inverse of fin
            spacing).
        k_fin: Quantity
            Thermal conductivity of the fin material. The default value applies
            to aluminum.
        num_circuits: int | None
            Number of refrigerant circuits. If None, the number of circuits is
            set equal to the number of tubes in the first row.
        """
        # Create the heat exchanger geometry of the evaporator.
        self.geometry = ContinuousFinStaggeredTubeBank(
            W_fro, H_fro, N_rows, S_trv, S_lon,
            D_ext, D_int, t_fin, N_fin, k_fin
        )
        # Create the heat exchanger of the superheating region.
        self.superheating_region = SuperheatingRegion(self.geometry, num_circuits)
        # Create the heat exchanger of the boiling region.
        self.boiling_region = BoilingRegion(self.geometry, num_circuits)
        # Determine total flow length of the evaporator.
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
        that air is cooled to the temperature of the refrigerant entering the
        evaporator and is fully saturated.
        """
        sat_air_out = HumidAir(Tdb=self.rfg_in.T, RH=Q_(100, 'pct'))
        Q_dot_max = self.air_m_dot * (self.air_in.h - sat_air_out.h)
        rfg_m_dot_max = Q_dot_max / (self.rfg_out.h - self.rfg_in.h)
        return rfg_m_dot_max
    
    def _calc_flow_length(self, rfg_m_dot: Quantity) -> Quantity:
        """
        Returns the required flow length of the evaporator needed to completely
        boil and to superheat the refrigerant to the required degree for the
        given mass flow rate of refrigerant.
        """
        # Superheating region: determine the required flow length to superheat
        # the refrigerant from saturated to superheated vapor having the
        # required degree of superheating.
        L_flow_superheat, air_out = None, None
        for k in range(4):
            L_flow_max = (5 ** k) * self.L_flow
            try:
                L_flow_superheat, air_out = self.superheating_region.solve(
                    rfg_m_dot=rfg_m_dot,
                    L_flow_max=L_flow_max
                )
            except SuperheatingError:
                logger.debug('Try again...')
                continue
            else:
                break
        if L_flow_superheat is None:
            raise SuperheatingError(
                "Required flow length of superheating region "
                "could not be determined."
            ) from None
        logger.debug(
            f"Required superheating flow length = "
            f"{L_flow_superheat.to('mm'):~P.3f}"
        )
        # Boiling region: determine the required flow length to completely boil
        # the refrigerant from the given inlet mixture state to saturated vapor.
        L_flow_boil = None
        for k in range(4):
            L_flow_max = (5 ** k) * self.L_flow
            try:
                L_flow_boil = self.boiling_region.solve(
                    air_in=air_out,
                    rfg_m_dot=rfg_m_dot,
                    L_flow_max=L_flow_max
                )
            except BoilingError:
                logger.debug('Try again...')
                continue
            else:
                break
        if L_flow_boil is None:
            raise BoilingError(
                "Required flow length of boiling region "
                "could not be determined"
            ) from None
        logger.debug(
            f"Required boiling flow length = "
            f"{L_flow_boil.to('mm'):~P.3f}"
        )
        L_flow = L_flow_superheat + L_flow_boil
        return L_flow
    
    def __fun__(self, rfg_m_dot: float, counter: list[int]) -> float:
        """
        Returns the deviation between the required flow length and the actual
        flow length of the evaporator.

        The required flow length is the flow length the evaporator should have
        to completely boil, and to superheat the refrigerant to the required
        degree of superheat at the evaporator exit.
        """
        rfg_m_dot = Q_(rfg_m_dot, 'kg / hr')
        i = counter[0]
        logger.debug(
            f"Iteration {i + 1}"
        )
        logger.debug(
            f"Try with refrigerant mass flow rate = {rfg_m_dot:~P.3f}"
        )
        # Given the mass flow rate of refrigerant, calculate the flow length
        # the evaporator should have to completely boil and superheat the
        # refrigerant to the required degree of superheat at the evaporator
        # exit.
        L_flow = self._calc_flow_length(rfg_m_dot)
        # Determine the deviation between the required flow length and the
        # actual flow length of the evaporator.
        dev = (L_flow - self.L_flow).to('mm')
        logger.debug(
            f"Total flow length = {L_flow.to('mm'):~P.3f} "
            f"(deviation = {dev:~P.3f})"
        )
        counter[0] += 1
        return dev.m

    def solve(
        self,
        air_in: HumidAir,
        air_m_dot: Quantity,
        rfg_in: FluidState,
        dT_sh: Quantity,
        rfg_m_dot: Quantity | None = None
    ) -> Quantity | tuple[Quantity, int]:
        """
        This function can have two different applications:
        
        1.  Analysis:
        If the mass flow rate of refrigerant is not specified (`rfg_m_dot` is
        `None`), solves for the refrigerant mass flow rate that is required
        under the given operating conditions to turn the liquid/vapor mixture
        entering the evaporator into superheated vapor with the required degree
        of superheating leaving the evaporator.
        
        2.  Design:
        If the mass flow rate of refrigerant is specified, solves for the
        required flow length of the evaporator so that the refrigerant leaves 
        the evaporator at the required degree of superheat.

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
        rfg_m_dot: optional
            The mass flow rate of refrigerant.

        Returns
        -------
        In case of analysis:
            The required mass flow rate of refrigerant.
        In case of design:
            The required flow length of the evaporator and the corresponding 
            number of rows.
        
        Raises
        ------
        SuperheatingError:
            If the required flow length for superheating the refrigerant to the
            required degree could not be determined.
        BoilingError:
            If the required flow length for boiling the refrigerant could not be
            determined.
        EvaporatorError:
            If the required mass flow rate of refrigerant for superheating the 
            refrigerant to the required degree could not be determined.  
        
        Notes
        -----
        The evaporator has two regions. First, the refrigerant is boiled into
        saturated vapor. Next, the refrigerant is superheated to the degree
        set on the expansion device.

        Under steady-state operation, the mass flow rate of refrigerant through
        the evaporator must be such that the degree of refrigerant superheating
        is maintained to the setting on the expansion device.

        The required flow length needed to superheat the saturated refrigerant
        vapor to the required degree depends on the refrigerant mass flow rate.
        The required flow length needed to boil the refrigerant from the
        specified inlet state into saturated vapor also depends on this mass
        flow rate.

        The required flow length of the superheating region and of the boiling
        region for a given mass flow rate of refrigerant is found when the heat
        transfer rate through the heat exchanger core in the superheating region
        and in the boiling region balances with the heat rate absorbed by the
        refrigerant in these regions.
        The entering and leaving states of the refrigerant in the superheating
        and in the boiling region are all known, and for a given mass flow rate
        of refrigerant it can then be determined what heat rate the refrigerant
        must absorb to reach these states.

        To solve the first application problem, the mass flow rate of
        refrigerant needs to be found for which the sum of the required
        superheating flow length and required boiling flow length is equal to
        the given, actual flow length of the evaporator.
        
        To determine the required flow length of the superheating region and of
        the boiling region in order to balance the heat transfer rate and the
        heat absorption rate of the refrigerant, Scipy's root-finding algorithm
        `root_scalar()` is used.
        This algorithm needs a minimum and a maximum limit between which it 
        searches for a flow length so that the deviation between the required
        heat rate that the refrigerant must absorb in the superheating/boiling
        region and the heat transfer rate through the heat exchanger core 
        becomes (near) zero. The function uses the number of rows specified by
        the user to determine a maximum limit for the flow length. When a
        `SuperheatingError` or `BoilingError` is raised, it could indicate that
        the given number of rows should be increased.
        
        To determine the refrigerant mass flow rate for which the required flow
        length of the evaporator is equal to its actual flow length, Scipy's
        root-finding algorithm `root_scalar` is also used. The maximum limit for
        the refrigerant mass flow rate is determined as 90 % of the
        theoretically maximum heat rate that can be transferred from the air to
        the refrigerant in the evaporator, and the minimum limit as 10 % of the 
        maximum limit.
        """
        # Assign/calculate what is known or can be calculated directly:
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
        # Mass flow rate of refrigerant is unknown, actual flow length of the
        # evaporator is given.
        # Determine the mass flow rate of refrigerant so that the required flow
        # length to boil and superheat the refrigerant is equal to the actual
        # flow length of the evaporator.
        if rfg_m_dot is None:
            rfg_m_dot_max = 0.9 * self._guess_rfg_m_dot_max().to('kg / hr').m
            rfg_m_dot_min = 0.1 * rfg_m_dot_max
            counter = [0]
            try:
                sol = optimize.root_scalar(
                    self.__fun__,
                    args=(counter,),
                    bracket=(rfg_m_dot_min, rfg_m_dot_max),
                    xtol=0.01,
                    rtol=0.001,
                    maxiter=20
                )
            except ValueError:
                message = (
                    f"The required refrigerant mass flow rate to superheat "
                    f"refrigerant with {dT_sh.to('K'):~P.1f} could not be "
                    f"determined."
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
                self.air_out = self.boiling_region.air_out
                self.Q_dot = (
                    self.boiling_region.Q_dot 
                    + self.superheating_region.Q_dot
                )
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
        # Mass flow rate of refrigerant is given, and the evaporator flow length
        # must be determined so that refrigerant leaves the evaporator at the
        # required degree of superheating.
        elif isinstance(rfg_m_dot, Quantity):
            L_flow = self._calc_flow_length(rfg_m_dot)
            N_rows = int(np.ceil(
                L_flow.to('mm').m / self.geometry.pitch_lon.to('mm').m
            ))
            self._set_flow_length(N_rows)

            self.rfg_m_dot = rfg_m_dot
            self.boiling_region.rfg_m_dot = self.rfg_m_dot
            self.superheating_region.rfg_m_dot = self.rfg_m_dot
            self.air_out = self.boiling_region.air_out
            self.Q_dot = (
                self.boiling_region.Q_dot
                + self.superheating_region.Q_dot
            )
            self.Q_dot_max = (
                self.boiling_region.Q_dot_max
                + self.superheating_region.Q_dot_max
            )
            self.eps = self.Q_dot / self.Q_dot_max
            self.air_dP = (
                self.boiling_region.dP_air
                + self.superheating_region.dP_air
            )

            return L_flow, N_rows
        return None

    def _set_flow_length(self, N_rows: int) -> None:
        self.L_flow = N_rows * self.geometry.pitch_lon
        self.geometry.length = self.L_flow

    @property
    def N_rows(self) -> int:
        return self.geometry.num_rows
