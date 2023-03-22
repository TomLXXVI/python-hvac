from typing import List, Optional, Tuple
from hvac import Quantity
from hvac.fluids import HumidAir, FluidState, CP_HUMID_AIR
from .core import CoilGeometry, AirStream, CoolantStream
from .core.control_volumes import DryControlVolume as DCV
from .core.heat_transfer.internal_flow import CondensationHeatTransfer, SinglePhaseHeatTransfer
from .air_cooler import AirCooler

Q_ = Quantity


class DryControlVolume(DCV):

    def _calculate_hi(self, Ts: Quantity, Tc: Quantity):
        if 0 < self.coolant_stream.outlet.x < 1:
            return CondensationHeatTransfer(
                Coolant=self.coolant_stream.coolant_type,
                Tc=Tc,  # self.coolant_stream.outlet.T,
                mc=self.coolant_stream.m,
                Di=self.coil_geometry.Di,
                x=self.coolant_stream.outlet.x,
                P=self.coolant_stream.outlet.P
            ).h
        else:
            return SinglePhaseHeatTransfer(
                coolant=self.coolant_stream.outlet,
                coil_geometry=self.coil_geometry,
                mc_tube=self.coolant_stream.m_tube,
                ma=self.air_stream.m,
                vfa=self.air_stream.vfa,
                air=self.air_stream.inlet,
                e=self.coolant_stream.e
            ).h


class AirCondenser:

    def __init__(
            self,
            coil_geometry: CoilGeometry,
            air_stream: AirStream,
            coolant_stream: CoolantStream,
            n_segments: int
    ) -> None:
        self.coil_geometry = coil_geometry
        self.air_stream = air_stream
        if self.air_stream.m is None:
            self.air_stream.m = self.air_stream.inlet.rho * self.air_stream.vfa * self.coil_geometry.Afa
        self.coolant_stream = coolant_stream
        self.n_segments = n_segments

        self.control_volumes: List[DryControlVolume] = []
        self.air_out: Optional[HumidAir] = None
        self.coolant_in: Optional[FluidState] = None
        self._section_printer = AirCondenser.SectionPrinter(self.n_segments)
        self._segment_printer = AirCondenser.SegmentPrinter()

    def analyze(self, display_info: bool = True) -> Tuple[HumidAir, FluidState]:
        self.control_volumes = []
        if display_info:
            self._section_printer.print_header()
        Ao = (
                self.coil_geometry.Afa *
                self.coil_geometry.Nr *
                self.coil_geometry.Ao_to_Afa
        )
        delta_Ao = Ao / self.n_segments
        air_stream = AirStream(
            inlet=self.air_stream.inlet,
            outlet=None,
            m=self.air_stream.m,
            vfa=self.air_stream.vfa
        )
        coolant_stream = CoolantStream(
            inlet=None,
            outlet=self.coolant_stream.outlet,
            m=self.coolant_stream.m,
            m_tube=self.coolant_stream.m_tube,
            coolant_type=self.coolant_stream.coolant_type,
            e=self.coolant_stream.e
        )
        for i in range(self.n_segments):
            cv = DryControlVolume(
                self.coil_geometry, air_stream,
                coolant_stream, delta_Ao
            )
            cv.analyze()
            self.control_volumes.append(cv)
            if display_info:
                self._section_printer.print_section(i, cv)
            air_stream = AirStream(
                inlet=cv.air_stream.outlet,
                outlet=None,
                m=cv.air_stream.m,
                vfa=cv.air_stream.vfa
            )
            coolant_stream = CoolantStream(
                inlet=None,
                outlet=cv.coolant_stream.inlet,
                m=cv.coolant_stream.m,
                m_tube=cv.coolant_stream.m_tube,
                coolant_type=cv.coolant_stream.coolant_type,
                e=cv.coolant_stream.e
            )
        self.air_out = self.control_volumes[-1].air_stream.outlet
        self.coolant_in = self.control_volumes[-1].coolant_stream.inlet
        return self.air_out, self.coolant_in

    def size(self, display_info: bool = True) -> Tuple[Quantity, Quantity, float]:
        self.control_volumes = []
        if display_info:
            self._segment_printer.print_header()
        hai = self.air_stream.inlet.h
        Tai = self.air_stream.inlet.Tdb.to('K')
        hao = self.air_stream.outlet.h
        Tao = self.air_stream.outlet.Tdb.to('K')
        dha = (hao - hai) / self.n_segments
        dTa = (Tao - Tai) / self.n_segments
        dhc = self.air_stream.m * dha / self.coolant_stream.m
        air_stream = AirStream(
            inlet=self.air_stream.inlet,
            outlet=HumidAir(Tdb=Tai + dTa, W=self.air_stream.inlet.W),
            m=self.air_stream.m,
            vfa=self.air_stream.vfa
        )
        coolant_stream = CoolantStream(
            inlet=self.coolant_stream.coolant_type(
                P=self.coolant_stream.outlet.P,
                h=self.coolant_stream.outlet.h + dhc
            ),
            outlet=self.coolant_stream.outlet,
            m=self.coolant_stream.m,
            m_tube=self.coolant_stream.m_tube,
            coolant_type=self.coolant_stream.coolant_type,
            e=self.coolant_stream.e
        )
        for i in range(self.n_segments):
            cv = DryControlVolume(self.coil_geometry, air_stream, coolant_stream, None)
            cv.size()
            self.control_volumes.append(cv)
            if display_info:
                self._segment_printer.print_segment(i, cv)
            air_stream = AirStream(
                inlet=cv.air_stream.outlet,
                outlet=HumidAir(Tdb=cv.air_stream.outlet.Tdb + dTa, W=cv.air_stream.outlet.W),
                m=cv.air_stream.m,
                vfa=cv.air_stream.vfa
            )
            coolant_stream = CoolantStream(
                inlet=cv.coolant_stream.coolant_type(
                    P=cv.coolant_stream.inlet.P,
                    h=cv.coolant_stream.inlet.h + dhc
                ),
                outlet=cv.coolant_stream.inlet,
                m=cv.coolant_stream.m,
                m_tube=cv.coolant_stream.m_tube,
                coolant_type=cv.coolant_stream.coolant_type,
                e=cv.coolant_stream.e
            )
        Ao = sum(cv.delta_Ao for cv in self.control_volumes)
        self.air_out = self.control_volumes[-1].air_stream.outlet
        self.coolant_in = self.control_volumes[-1].coolant_stream.inlet
        self.coil_geometry.Afa = self.air_stream.m / (self.air_stream.inlet.rho * self.air_stream.vfa)
        self.coil_geometry.Nr = float(Ao / (self.coil_geometry.Ao_to_Afa * self.coil_geometry.Afa))
        return Ao, self.coil_geometry.Afa, self.coil_geometry.Nr

    @property
    def Q(self) -> Quantity:
        Q = sum(cv.ht.delta_Q for cv in self.control_volumes)
        return abs(Q)

    def get_air_states(self) -> List[HumidAir]:
        air_state_list = [self.control_volumes[0].air_stream.inlet]
        air_state_list.extend(cv.air_stream.outlet for cv in self.control_volumes)
        return air_state_list

    def get_surface_temperatures(self) -> List[Quantity]:
        surf_state_list = [self.control_volumes[0].surface.Ti]
        surf_state_list.extend(cv.surface.To for cv in self.control_volumes)
        return surf_state_list

    def get_coolant_states(self) -> List[FluidState]:
        coolant_state_list = [self.control_volumes[0].coolant_stream.outlet]
        coolant_state_list.extend(cv.coolant_stream.inlet for cv in self.control_volumes)
        return coolant_state_list

    @property
    def Uo_avg(self) -> Quantity:
        Uo_avg = sum(cv.htp.Uo * cv.delta_Ao for cv in self.control_volumes) / self.Ao
        return Uo_avg

    @property
    def he_avg(self) -> Quantity:
        he_avg = sum(cv.htp.he * cv.delta_Ao for cv in self.control_volumes) / self.Ao
        return he_avg

    @property
    def eta_s_avg(self) -> Quantity:
        n = len(self.control_volumes)
        eta_s_avg = sum(cv.htp.eta_surf for cv in self.control_volumes) / n
        return eta_s_avg

    @property
    def hi_avg(self) -> Quantity:
        Ai_to_Ao = 1 / self.coil_geometry.Ao_to_Ai
        Ai = Ai_to_Ao * self.Ao
        hi_avg = sum(cv.htp.hi * Ai_to_Ao * cv. delta_Ao for cv in self.control_volumes) / Ai
        return hi_avg

    @property
    def Q_max(self) -> Quantity:
        # only applies to refrigerants; Q_max is based on condensation temperature (saturation temperature)
        # Tc = self.coolant_stream.coolant_type(P=self.coolant_stream.outlet.P, x=Q_(0, 'frac')).T
        delta_hc = self.coolant_in.h - self.coolant_stream.outlet.h
        delta_Tc = self.coolant_in.T - self.coolant_stream.outlet.T
        cpc = delta_hc / delta_Tc
        Cmin = min(self.air_stream.m * CP_HUMID_AIR, self.coolant_stream.m * cpc)
        Q_max = Cmin * (self.coolant_in.T - self.air_stream.inlet.Tdb)
        return Q_max

    @property
    def eps(self) -> Quantity:
        return self.Q / self.Q_max

    @property
    def Ao(self) -> Quantity:
        """Total external heat transfer surface area."""
        return self.coil_geometry.Ao

    @property
    def dPa(self) -> Quantity:
        """Air-side pressure drop."""
        Afa_to_Ac = 1 / self.coil_geometry.Ac_to_Afa
        Ao_to_Afa = self.coil_geometry.Ao_to_Afa
        Ao_to_Ac = Afa_to_Ac * Ao_to_Afa
        rho_in = self.air_stream.inlet.rho
        rho_out = self.air_stream.outlet.rho
        rho_avg = 0.5 * (rho_in + rho_out)
        G = self.air_stream.m / self.coil_geometry.Ac  # mass velocity
        f_fan_avg = sum(cv.get_fanning_friction_factor_air() for cv in self.control_volumes) / len(self.control_volumes)
        f1 = G ** 2 / (2 * rho_in)
        t1 = (1 + self.coil_geometry.Ac_to_Afa ** 2) * (rho_in / rho_out - 1)
        t2 = f_fan_avg * self.coil_geometry.Nr * Ao_to_Ac * rho_in / rho_avg
        dP = f1 * (t1 + t2)
        return dP

    class SectionPrinter(AirCooler.SectionPrinter):
        pass

    class SegmentPrinter(AirCooler.SegmentPrinter):
        pass
