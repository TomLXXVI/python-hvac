from typing import Optional, Dict, Union
from abc import ABC, abstractmethod
from enum import Enum
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState, CP_HUMID_AIR, Water, Ice, STANDARD_PRESSURE
from .geometry import CoilGeometry
from .components import AirStream, CoolantStream, Surface, HeatTransferParameters, HeatTransfer
from .fins import FinnedSurface
from .heat_transfer.internal_flow import BoilingHeatTransfer, SinglePhaseHeatTransfer
from .heat_transfer.external_air_flow import WetSurfaceHeatTransfer, DrySurfaceHeatTransfer
from .friction_factor.external_air_flow import WetSurfaceFriction, DrySurfaceFriction


Q_ = Quantity


class AbstractControlVolume(ABC):

    def __init__(
        self,
        coil_geometry: CoilGeometry,
        air_stream: AirStream,
        coolant_stream: CoolantStream,
        delta_Ao: Optional[Quantity]
    ) -> None:
        self.coil_geometry = coil_geometry
        self.air_stream = air_stream
        self.coolant_stream = coolant_stream
        self.delta_Ao = delta_Ao
        self.surface = Surface()
        self.htp = HeatTransferParameters()
        self.ht = HeatTransfer()

    def get_surface_temperature(self) -> Quantity:
        Tsi = self._calculate_surface_temperature(self.air_stream.inlet, self.coolant_stream.outlet)
        return Tsi

    @abstractmethod
    def size(self) -> None:
        self.surface.Ti = self._calculate_surface_temperature(
            self.air_stream.inlet,
            self.coolant_stream.outlet
        )
        self.surface.To = self._calculate_surface_temperature(
            self.air_stream.outlet,
            self.coolant_stream.inlet
        )
        self._calculate_heat_transfer_params()
        self.ht.delta_Q = self.coolant_stream.m * (self.coolant_stream.outlet.h - self.coolant_stream.inlet.h)
        # extended in derived classes

    @abstractmethod
    def analyze(self) -> None:
        self.surface.Ti = self._calculate_surface_temperature(
            self.air_stream.inlet,
            self.coolant_stream.outlet
        )
        self._calculate_heat_transfer_params()
        # extended in derived classes

    @abstractmethod
    def _calculate_surface_temperature(self, air: HumidAir, coolant: FluidState) -> Quantity:
        pass

    @abstractmethod
    def _calculate_he(self, air: HumidAir) -> Quantity:
        pass

    def _calculate_hi(self, Ts: Quantity, Tc: Quantity) -> Quantity:
        if 0 < self.coolant_stream.outlet.x < 1:
            return BoilingHeatTransfer(
                Coolant=self.coolant_stream.coolant_type, mc=self.coolant_stream.m,
                Tc=Tc, Ts=Ts,  # self.coolant_stream.outlet.T
                Di=self.coil_geometry.Di, x=self.coolant_stream.outlet.x
            ).h
        else:
            return SinglePhaseHeatTransfer(
                self.coolant_stream.outlet, self.coil_geometry,
                self.coolant_stream.m_tube, self.air_stream.m,
                self.air_stream.vfa, self.air_stream.inlet,
                self.coolant_stream.e
            ).h

    @abstractmethod
    def _calculate_xi(self, Ts: Quantity, Tc: Quantity) -> Quantity:
        pass

    @abstractmethod
    def _calculate_eta_surf(self, he: Quantity, xi: Optional[Quantity] = None) -> Quantity:
        pass

    @abstractmethod
    def _calculate_Uo(self, he: Quantity, eta_surf: Quantity, hi: Quantity, xi: Optional[Quantity] = None) -> Quantity:
        pass

    def _calculate_R(self, he: Quantity, eta_surf: Quantity, hi: Quantity) -> Quantity:
        return he * eta_surf * self.coil_geometry.Ao_to_Ai / hi

    def _calculate_heat_transfer_params(self):
        # heat transfer parameters are based on air inlet section
        self.htp.he = self._calculate_he(self.air_stream.inlet)
        self.htp.xi = self._calculate_xi(self.surface.Ti, self.coolant_stream.outlet.T)
        self.htp.eta_surf = self._calculate_eta_surf(self.htp.he, self.htp.xi)
        self.htp.hi = self._calculate_hi(self.surface.Ti, self.coolant_stream.outlet.T)
        self.htp.Uo = self._calculate_Uo(self.htp.he, self.htp.eta_surf, self.htp.hi, self.htp.xi)
        self.htp.R = self._calculate_R(self.htp.he, self.htp.eta_surf, self.htp.hi)

    @abstractmethod
    def _calculate_delta_Ao(self) -> None:
        pass

    @abstractmethod
    def _solve_analysis_equations(self) -> Dict[str, Quantity]:
        pass

    @abstractmethod
    def get_fanning_friction_factor_air(self) -> float:
        pass


class WetControlVolume(AbstractControlVolume):

    class SystemOfEquations:

        def __init__(self, **kwargs):
            self.ma = kwargs.get('ma').to('kg / s').m
            self.mc = kwargs.get('mc').to('kg / s').m

            self.cpa = CP_HUMID_AIR.to('J / (kg * K)').m
            self.Tai = kwargs.get('Tai').to('K').m
            self.Wai = kwargs.get('Wai').to('kg / kg').m

            self.Tsi = kwargs.get('Tsi').to('K').m
            self.Wsi = kwargs.get('Wsi').to('kg / kg').m

            self.hco = kwargs.get('hco').to('J / kg').m
            self.Pc = kwargs.get('Pc')

            self.alpha_e = kwargs.get('alpha_e').to('W / (m ** 2 * K)').m
            self.eta_surf = kwargs.get('eta_surf').to('frac').m
            self.R = kwargs.get('R').to('K * kg / J').m

            self.delta_Ao = kwargs.get('delta_Ao').to('m ** 2').m

            self.coolant_type = kwargs.get('coolant_type')

        @staticmethod
        def _get_hfg(Ts: float):
            # get latent heat of vaporization/condensation at surface temperature Ts
            try:
                hf = Water(T=Q_(Ts, 'K'), x=Q_(0, 'frac')).h.to('J /kg').m
                hg = Water(T=Q_(Ts, 'K'), x=Q_(1, 'frac')).h.to('J /kg').m
            except (ValueError, NotImplementedError):
                hf = Ice(T=Q_(273.16, 'K'), P=STANDARD_PRESSURE).h.to('J /kg').m
                hg = Water(T=Q_(273.16, 'K'), P=STANDARD_PRESSURE).h.to('J /kg').m
            return hg - hf

        # air-side sensible heat balance
        def eq1(self, delta_Qs: float, Tao: float) -> float:
            out = delta_Qs - self.ma * self.cpa * (self.Tai - Tao)
            return out

        # air-side latent heat balance
        def eq2(self, delta_Ql: float, Wao: float, Ts_avg: float) -> float:
            hfg = self._get_hfg(Ts_avg)
            out = delta_Ql - self.ma * (self.Wai - Wao) * hfg
            return out

        # air-side sensible heat transfer to surface
        def eq3(self, delta_Qs: float, Ta_avg: float, Ts_avg: float) -> float:
            out = delta_Qs - self.alpha_e * self.eta_surf * self.delta_Ao * (Ta_avg - Ts_avg)
            return out

        # air-side latent heat transfer to surface
        def eq4(self, delta_Ql: float, Wa_avg: float, Ts_avg: float) -> float:
            Ws_avg = HumidAir(Tdb=Q_(Ts_avg, 'K'), RH=Q_(100, 'pct')).W.to('kg / kg').m
            hfg = self._get_hfg(Ts_avg)
            out = delta_Ql - (self.alpha_e / self.cpa) * self.eta_surf * self.delta_Ao * (Wa_avg - Ws_avg) * hfg
            return out

        # surface temperature at air outlet
        def eq5(self, Tao: float, Wao: float, Tso: float, hci: float) -> float:
            try:
                hao = HumidAir(Tdb=Q_(Tao, 'K'), W=Q_(Wao, 'kg / kg')).h.to('J / kg').m
            except ValueError:
                hao = HumidAir(Tdb=Q_(Tao, 'K'), RH=Q_(100, 'pct')).h.to('J / kg').m
            # Note: here we assume that the coolant pressure remains constant
            Tci = self.coolant_type(P=self.Pc, h=Q_(hci, 'J / kg')).T.to('K').m
            hso = HumidAir(Tdb=Q_(Tso, 'K'), RH=Q_(100, 'pct')).h.to('J / kg').m
            out = Tso - (Tci + self.R * (hao - hso))
            return out

        # heat balance: total heat = sensible heat + latent heat
        @staticmethod
        def eq6(delta_Q: float, delta_Qs: float, delta_Ql: float) -> float:
            out = delta_Q - (delta_Qs + delta_Ql)
            return out

        # coolant-side heat balance
        def eq7(self, delta_Q, hci: float) -> float:
            out = delta_Q - self.mc * (self.hco - hci)
            return out

        def eq8(self, Ta_avg: float, Tao: float) -> float:
            out = Ta_avg - 0.5 * (self.Tai + Tao)
            return out

        def eq9(self, Wa_avg: float, Wao: float) -> float:
            out = Wa_avg - 0.5 * (self.Wai + Wao)
            return out

        def eq10(self, Ts_avg: float, Tso: float) -> float:
            out = Ts_avg - 0.5 * (self.Tsi + Tso)
            return out

        def eq11(self, hc_avg: float, hci: float) -> float:
            out = hc_avg - 0.5 * (self.hco + hci)
            return out

        def equations(self, unknowns: np.ndarray) -> np.ndarray:
            delta_Qs = unknowns[0]
            delta_Ql = unknowns[1]
            delta_Q = unknowns[2]
            Tao = unknowns[3]
            Wao = unknowns[4]
            Tso = unknowns[5]
            hci = unknowns[6]
            Ta_avg = unknowns[7]
            Wa_avg = unknowns[8]
            Ts_avg = unknowns[9]
            hc_avg = unknowns[10]

            out1 = self.eq1(delta_Qs, Tao)
            out2 = self.eq2(delta_Ql, Wao, Ts_avg)
            out3 = self.eq3(delta_Qs, Ta_avg, Ts_avg)
            out4 = self.eq4(delta_Ql, Wa_avg, Ts_avg)
            out5 = self.eq5(Tao, Wao, Tso, hci)
            out6 = self.eq6(delta_Q, delta_Qs, delta_Ql)
            out7 = self.eq7(delta_Q, hci)
            out8 = self.eq8(Ta_avg, Tao)
            out9 = self.eq9(Wa_avg, Wao)
            out10 = self.eq10(Ts_avg, Tso)
            out11 = self.eq11(hc_avg, hci)

            return np.array(
                [out1, out2, out3, out4, out5, out6,
                 out7, out8, out9, out10, out11]
            )

        def guess(self) -> np.ndarray:
            Tao = self.Tai - 1.0
            Tso = self.Tsi - 0.01
            Ta_avg = 0.5 * (self.Tai + Tao)
            Ts_avg = 0.5 * (self.Tsi + Tso)
            delta_Qs = (
                    self.alpha_e *
                    self.eta_surf *
                    self.delta_Ao *
                    (Ta_avg - Ts_avg)
            )
            Wao = self.Wai - 0.0005
            Wso = self.Wsi - 0.0005
            Wa_avg = 0.5 * (self.Wai + Wao)
            Ws_avg = 0.5 * (self.Wsi + Wso)
            hfg = self._get_hfg(Ts_avg)
            delta_Ql = (
                    (self.alpha_e / self.cpa) * self.eta_surf *
                    self.delta_Ao * (Wa_avg - Ws_avg) * hfg
            )
            delta_Q = delta_Qs + delta_Ql
            hci = self.hco - delta_Q / self.mc
            hc_avg = 0.5 * (self.hco + hci)
            return np.array(
                [delta_Qs, delta_Ql, delta_Q, Tao, Wao, Tso,
                 hci, Ta_avg, Wa_avg, Ts_avg, hc_avg]
            )

        def solve(self) -> Dict[str, Quantity]:
            roots = fsolve(self.equations, self.guess())
            return {
                'delta_Qs': Q_(roots[0], 'W'),
                'delta_Ql': Q_(roots[1], 'W'),
                'delta_Q': Q_(roots[2], 'W'),
                'Tao': Q_(roots[3], 'K'),
                'Wao': Q_(roots[4], 'kg / kg'),
                'Tso': Q_(roots[5], 'K'),
                'hci': Q_(roots[6], 'J/kg'),
                'Ta_avg': Q_(roots[7], 'K'),
                'Wa_avg': Q_(roots[8], 'kg / kg'),
                'Ts_avg': Q_(roots[9], 'K'),
                'hc_avg': Q_(roots[10], 'J / kg')
            }

    def _calculate_surface_temperature(self, air: HumidAir, coolant: FluidState) -> Quantity:

        def eq(unknowns: np.ndarray) -> np.ndarray:
            Ts = Q_(unknowns[0], 'K')
            Tc = coolant.T.to('K')
            he = self._calculate_he(air)
            xi = self._calculate_xi(Ts, Tc)
            eta_surf = self._calculate_eta_surf(he, xi)
            hi = self._calculate_hi(Ts, Tc)
            R = self._calculate_R(he, eta_surf, hi)
            ha = air.h
            hs = HumidAir(Tdb=Ts, RH=Q_(100, 'pct')).h
            out = Tc + R * (ha - hs) - Ts
            return np.array(out.m)

        Ts_ini = (coolant.T.to('K') + air.Tdb.to('K')) / 2
        roots = fsolve(eq, np.array([Ts_ini.m]))
        Ts = Q_(roots[0], 'K')
        return Ts

    def _calculate_he(self, air: HumidAir) -> Quantity:
        return WetSurfaceHeatTransfer(
            air, self.air_stream.vfa, self.coil_geometry.Ac_to_Afa,
            self.coil_geometry.Ao_to_Afa, self.coil_geometry.Ap_to_Afa,
            self.coil_geometry.Do, self.coil_geometry.sf, self.coil_geometry.tf
        ).h_tot

    def _calculate_xi(self, Ts: Quantity, Tc: Quantity) -> Quantity:
        hs = HumidAir(Tdb=Ts, RH=Q_(100, 'pct')).h
        hc_star = HumidAir(Tdb=Tc, RH=Q_(100, 'pct')).h
        xi = (hs - hc_star) / (Ts - Tc)
        return xi

    def _calculate_eta_surf(self, he: Quantity, xi: Optional[Quantity] = None) -> Quantity:
        eta = FinnedSurface(self.coil_geometry, he * xi).efficiency
        return eta

    def _calculate_Uo(self, he: Quantity, eta_surf: Quantity, hi: Quantity, xi: Optional[Quantity] = None) -> Quantity:
        Uo = 1 / (
                1 / (he * eta_surf) +
                xi * self.coil_geometry.Ao_to_Ai / hi
        )
        return Uo

    def size(self) -> None:
        super().size()
        self._calculate_delta_Ao()
        self._calculate_air_out()
        self._set_heat_transfers()

    def _calculate_delta_Ao(self) -> None:
        ha_avg = 0.5 * (self.air_stream.inlet.h + self.air_stream.outlet.h)
        hs_avg = 0.5 * (self.surface.inlet.h + self.surface.outlet.h)
        delta_Ao = self.ht.delta_Q / (self.htp.he * self.htp.eta_surf * (ha_avg - hs_avg))
        self.delta_Ao = delta_Ao

    def _calculate_air_out(self) -> None:
        Ws_avg = 0.5 * (self.surface.inlet.W + self.surface.outlet.W)
        k = 0.5 * self.htp.he * self.delta_Ao
        Wai = self.air_stream.inlet.W
        Wao = (2 * k * Ws_avg + (self.air_stream.m - k) * Wai) / (self.air_stream.m + k)
        air_out = HumidAir(h=self.air_stream.outlet.h, W=Wao)
        self.air_stream.outlet = air_out

    def _set_heat_transfers(self) -> None:
        self.ht.delta_Qs = self.air_stream.m * CP_HUMID_AIR * (self.air_stream.inlet.Tdb - self.air_stream.outlet.Tdb)
        self.ht.delta_Ql = self.ht.delta_Q - self.ht.delta_Qs
        self.ht.delta_mw = self.air_stream.m * (self.air_stream.inlet.W - self.air_stream.outlet.W)

    def _solve_analysis_equations(self) -> Dict[str, Quantity]:
        eqs = WetControlVolume.SystemOfEquations(
            ma=self.air_stream.m,
            mc=self.coolant_stream.m,
            Tai=self.air_stream.inlet.Tdb,
            Wai=self.air_stream.inlet.W,
            Tsi=self.surface.inlet.Tdb,
            Wsi=self.surface.inlet.W,
            hco=self.coolant_stream.outlet.h,
            Pc=self.coolant_stream.outlet.P,
            alpha_e=self.htp.he * CP_HUMID_AIR,
            eta_surf=self.htp.eta_surf,
            R=self.htp.R,
            delta_Ao=self.delta_Ao,
            coolant_type=self.coolant_stream.coolant_type
        )
        sol = eqs.solve()
        return sol

    def analyze(self) -> None:
        super().analyze()
        sol = self._solve_analysis_equations()
        try:
            self.air_stream.outlet = HumidAir(Tdb=sol['Tao'], W=sol['Wao'])
            _ = self.air_stream.outlet.RH
        except ValueError:
            self.air_stream.outlet = HumidAir(Tdb=sol['Tao'], RH=Q_(100, 'pct'))
        self.surface.To = sol['Tso']
        self.coolant_stream.inlet = self.coolant_stream.coolant_type(P=self.coolant_stream.outlet.P, h=sol['hci'])
        self.ht.delta_Q = sol['delta_Q']
        self._set_heat_transfers()

    def get_fanning_friction_factor_air(self) -> float:
        return WetSurfaceFriction(
            self.air_stream.inlet, self.air_stream.vfa, self.coil_geometry.Ac_to_Afa,
            self.coil_geometry.Ao_to_Afa, self.coil_geometry.Ap_to_Afa,
            self.coil_geometry.sv, self.coil_geometry.Do, self.coil_geometry.sf,
            self.coil_geometry.tf
        ).f_fan


class DryControlVolume(AbstractControlVolume):

    class SystemOfEquations:

        def __init__(self, **kwargs):
            self.coolant_type = kwargs.get('coolant_type')
            self.ma = kwargs.get('ma').to('kg / s').m
            self.cpa = CP_HUMID_AIR.to('J / (kg * K)').m
            self.Tai = kwargs.get('Tai').to('K').m
            self.Tsi = kwargs.get('Tsi').to('K').m
            self.mc = kwargs.get('mc').to('kg / s').m
            self.hco = kwargs.get('hco').to('J / kg').m
            self.Pc = kwargs.get('Pc').to('Pa')
            self.alpha_e = kwargs.get('alpha_e').to('W / (m ** 2 * K)').m
            self.eta_surf = kwargs.get('eta_surf').to('frac').m
            self.R = kwargs.get('R').to('frac').m
            self.delta_Ao = kwargs.get('delta_Ao').to('m ** 2').m

        def eq1(self, delta_Q, Tao, Tso):
            Ta_avg = 0.5 * (self.Tai + Tao)
            Ts_avg = 0.5 * (self.Tsi + Tso)
            lhs = self.alpha_e * self.eta_surf * self.delta_Ao * (Ta_avg - Ts_avg)
            rhs = delta_Q
            return rhs - lhs

        def eq2(self, delta_Q, Tao):
            lhs = self.ma * self.cpa * (self.Tai - Tao)
            rhs = delta_Q
            return rhs - lhs

        def eq3(self, delta_Q, hci):
            lhs = self.mc * (self.hco - hci)
            rhs = delta_Q
            return rhs - lhs

        def eq4(self, Tao, hci, Tso):
            Tci = self.coolant_type(h=Q_(hci, 'J / kg'), P=self.Pc).T.to('K').m
            lhs = (self.R * Tao + Tci) / (1 + self.R)
            rhs = Tso
            return rhs - lhs

        def equations(self, unknowns: np.ndarray) -> np.ndarray:
            delta_Q = unknowns[0]
            Tao = unknowns[1]
            Tso = unknowns[2]
            hci = unknowns[3]

            out1 = self.eq1(delta_Q, Tao, Tso)
            out2 = self.eq2(delta_Q, Tao)
            out3 = self.eq3(delta_Q, hci)
            out4 = self.eq4(Tao, hci, Tso)

            return np.array([out1, out2, out3, out4])

        def guess(self) -> np.ndarray:
            delta_Q_ini = 100
            Tao_ini = self.Tai - delta_Q_ini / (self.ma * self.cpa)
            hci_ini = self.hco - delta_Q_ini / self.mc
            Tci_ini = self.coolant_type(h=Q_(hci_ini, 'J / kg'), P=self.Pc).T.to('K').m
            Tso_ini = (self.R * Tao_ini + Tci_ini) / (1 + self.R)
            return np.array([delta_Q_ini, Tao_ini, Tso_ini, hci_ini])

        def solve(self) -> Dict[str, Quantity]:
            sol = fsolve(self.equations, self.guess())
            delta_Q = Q_(sol[0], 'W')
            Tao = Q_(sol[1], 'K')
            Tso = Q_(sol[2], 'K')
            hci = Q_(sol[3], 'J / kg')
            return {'delta_Q': delta_Q, 'Tao': Tao, 'Tso': Tso, 'hci': hci}

    def _calculate_surface_temperature(self, air: HumidAir, coolant: FluidState) -> Quantity:

        def eq(unknown: float) -> float:
            Ta = air.Tdb.to('K')
            Ts = Q_(unknown, 'K')
            Tc = coolant.T.to('K')
            he = self._calculate_he(air)
            eta_surf = self._calculate_eta_surf(he)
            hi = self._calculate_hi(Ts, Tc)
            R = self._calculate_R(he, eta_surf, hi)
            out = (R * Ta + Tc) / (1 + R) - Ts
            return out.m

        # Ts_ini = coolant.T.to('K') + Q_(0.1, 'K')
        # roots = fsolve(eq, np.array([Ts_ini.m]))
        # Ts = Q_(roots[0], 'K')

        sol = root_scalar(eq, method='brentq', bracket=[coolant.T.to('K').m, air.Tdb.to('K').m])
        Ts = Q_(sol.root, 'K')
        return Ts

    def _calculate_he(self, air: HumidAir) -> Quantity:
        return DrySurfaceHeatTransfer(
            air, self.air_stream.vfa, self.coil_geometry.Ac_to_Afa,
            self.coil_geometry.Ao_to_Afa, self.coil_geometry.Ap_to_Afa,
            self.coil_geometry.Do
        ).h

    def _calculate_xi(self, Ts: Quantity, Tc: Quantity) -> Quantity:
        return None

    def _calculate_eta_surf(self, he: Quantity, xi: Optional[Quantity] = None) -> Quantity:
        eta = FinnedSurface(self.coil_geometry, he).efficiency
        return eta

    def _calculate_Uo(self, he: Quantity, eta_surf: Quantity, hi: Quantity, xi: Optional[Quantity] = None) -> Quantity:
        Uo = 1 / (
                1 / (he * eta_surf) +
                self.coil_geometry.Ao_to_Ai / hi
        )
        return Uo

    def size(self) -> None:
        super().size()
        self._calculate_delta_Ao()

    def _calculate_delta_Ao(self) -> None:
        Ta_avg = 0.5 * (self.air_stream.inlet.Tdb + self.air_stream.outlet.Tdb)
        Ts_avg = 0.5 * (self.surface.inlet.Tdb + self.surface.outlet.Tdb)
        delta_Ao = self.ht.delta_Q / (self.htp.he * self.htp.eta_surf * (Ta_avg - Ts_avg))
        self.delta_Ao = delta_Ao

    def _solve_analysis_equations(self) -> Dict[str, Quantity]:
        eqs = DryControlVolume.SystemOfEquations(
            coolant_type=self.coolant_stream.coolant_type,
            ma=self.air_stream.m,
            Tai=self.air_stream.inlet.Tdb,
            Tsi=self.surface.Ti,
            mc=self.coolant_stream.m,
            hco=self.coolant_stream.outlet.h,
            Pc=self.coolant_stream.outlet.P,
            alpha_e=self.htp.he,
            eta_surf=self.htp.eta_surf,
            R=self.htp.R,
            delta_Ao=self.delta_Ao
        )
        sol = eqs.solve()
        return sol

    def analyze(self) -> None:
        super().analyze()
        sol = self._solve_analysis_equations()
        self.air_stream.outlet = HumidAir(Tdb=sol['Tao'], W=self.air_stream.inlet.W)
        self.coolant_stream.inlet = self.coolant_stream.coolant_type(P=self.coolant_stream.outlet.P, h=sol['hci'])
        self.surface.To = sol['Tso']
        self.ht.delta_Q = sol['delta_Q']

    def get_fanning_friction_factor_air(self) -> float:
        return DrySurfaceFriction(
            self.air_stream.inlet, self.air_stream.vfa, self.coil_geometry.Ac_to_Afa,
            self.coil_geometry.Ao_to_Afa, self.coil_geometry.Ap_to_Afa,
            self.coil_geometry.sv, self.coil_geometry.Do, self.coil_geometry.sf,
            self.coil_geometry.tf
        ).f_fan


class CoilSurfaceCondition(Enum):
    WET = 'wet'
    DRY = 'dry'


class ControlVolume:

    @classmethod
    def create(
            cls,
            coil_geometry: CoilGeometry,
            air_stream: AirStream,
            coolant_stream: CoolantStream,
            delta_Ao: Optional[Quantity] = None
    ) -> Union[WetControlVolume, DryControlVolume]:
        cv_wet = WetControlVolume(coil_geometry, air_stream, coolant_stream, delta_Ao)
        cv_dry = DryControlVolume(coil_geometry, air_stream, coolant_stream, delta_Ao)
        surf_cond = cls._get_surf_condition(air_stream.inlet, cv_wet, cv_dry)
        if surf_cond == CoilSurfaceCondition.WET:
            return cv_wet
        else:
            return cv_dry

    @staticmethod
    def _get_surf_condition(
            air_in: HumidAir,
            cv_wet: WetControlVolume,
            cv_dry: DryControlVolume
    ) -> CoilSurfaceCondition:
        try:
            Ts_wet = cv_wet.get_surface_temperature()
        except TypeError:
            return CoilSurfaceCondition.DRY
        Ts_dry = cv_dry.get_surface_temperature()
        Ts = min(Ts_wet, Ts_dry)
        if Ts <= air_in.Tdp:
            return CoilSurfaceCondition.WET
        else:
            return CoilSurfaceCondition.DRY
