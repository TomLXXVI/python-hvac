#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass, field

from hvac import Quantity
from hvac.fluids.humid_air import HumidAir


Q_ = Quantity


@dataclass
class CondenserParams:
    air_in: HumidAir
    Q_dot: Quantity
    TD_in: Quantity = field(default=Q_(15, 'K'))  # typical 10 to 15 K
    TD_out: Quantity = field(default=Q_(5, 'K'))  # typical 3 to 5 K
    v_fa: Quantity = field(default=Q_(2, 'm / s'))
    aspect_ratio: Quantity = field(default=1 / 2)

    @property
    def T_cnd(self) -> Quantity:
        """Condensing temperature of refrigerant."""
        return self.air_in.Tdb + self.TD_in

    @property
    def air_out(self) -> HumidAir:
        """State of air at the condenser air outlet."""
        air_out_T = self.T_cnd - self.TD_out
        air_out = HumidAir(Tdb = air_out_T, W=self.air_in.W)
        return air_out

    @property
    def air_m_dot(self) -> Quantity:
        """Mass flow rate of air through condenser."""
        if self.Q_dot is None:
            raise ValueError("condenser capacity must be set first")
        return self.Q_dot / (self.air_out.h - self.air_in.h)

    @property
    def air_V_dot(self) -> Quantity:
        """Air volume flow rate referred to air inlet."""
        return self.air_m_dot / self.air_in.rho

    @property
    def A_fa(self) -> Quantity:
        """Air-side face area."""
        return self.air_V_dot / self.v_fa

    @property
    def front_dims(self) -> tuple[Quantity, Quantity]:
        """Width and height of the air-side face area."""
        width = (self.A_fa / self.aspect_ratio) ** 0.5
        height = self.A_fa / width
        return width, height

    @property
    def width(self) -> Quantity:
        """Width of the air-side face area."""
        width, _ = self.front_dims
        return width

    @property
    def height(self) -> Quantity:
        """Height of the air-side face area."""
        _, height = self.front_dims
        return height

    def __str__(self):
        lines = [
            f"air in: {str(self.air_in)}",
            f"air out: {str(self.air_out)}",
            f"TD in: {self.TD_in.to('K'):~P.1f}",
            f"TD out: {self.TD_out.to('K'):~P.1f}",
            f"T condensing: {self.T_cnd.to('degC'):~P.1f}",
            f"aspect ratio air face area: 1:{1 / self.aspect_ratio:.0f}"
        ]
        if self.Q_dot is not None:
            lines_ext = [
                "---------------",
                f"heat transfer rate: {self.Q_dot.to('kW'):~P.3f}",
                f"air m-rate: {self.air_m_dot.to('kg / hr'):~P.3f}",
                f"air v-rate: {self.air_V_dot.to('m ** 3 / hr'):~P.3f}",
                f"air face velocity: {self.v_fa.to('m / s'):~P.1f}",
                f"air face width: {self.width.to('mm'):~P.0f}",
                f"air face height: {self.height.to('mm'):~P.0f}"
            ]
            lines.extend(lines_ext)
        return "\n".join(lines)


@dataclass
class EvaporatorParams:
    air_in: HumidAir
    TD_in: Quantity = field(default=Q_(10, 'K'))  # typical 5 to 10 K
    TD_out: Quantity = field(default=Q_(5, 'K'))  # typical 3 to 5 K
    Q_dot: Quantity | None = field(default=None)
    v_fa: Quantity = field(default=Q_(2, 'm / s'))
    aspect_ratio: Quantity = field(default=1 / 2)
    superheat_offset: Quantity = field(default=Q_(2, 'K'))  # difference between TD_in and dT_sh

    def __post_init__(self):
        self.TD_in.ito('K')
        self.TD_out.ito('K')
        if isinstance(self.Q_dot, Quantity): self.Q_dot.ito('W')
        self.v_fa.ito('m / s')
        self.superheat_offset.ito('K')

    @property
    def T_evp(self) -> Quantity:
        """Evaporating temperature of refrigerant."""
        return self.air_in.Tdb - self.TD_in

    @property
    def air_out(self) -> HumidAir:
        """State of air at the evaporator air outlet."""
        air_out_T = self.T_evp + self.TD_out
        air_out = HumidAir(Tdb=air_out_T, RH=Q_(100, 'pct'))  # guess
        return air_out

    @property
    def air_m_dot(self) -> Quantity:
        """Mass flow rate of air through evaporator."""
        if self.Q_dot is None:
            raise ValueError("evaporator capacity must be set first")
        return self.Q_dot / (self.air_in.h - self.air_out.h)

    @property
    def air_V_dot(self) -> Quantity:
        """Air volume flow rate referred to air inlet."""
        return self.air_m_dot / self.air_in.rho

    @property
    def A_fa(self) -> Quantity:
        """Air-side face area."""
        return self.air_V_dot / self.v_fa

    @property
    def front_dims(self) -> tuple[Quantity, Quantity]:
        """Width and height of the air-side face area."""
        width = (self.A_fa / self.aspect_ratio) ** 0.5
        height = self.A_fa / width
        return width, height

    @property
    def width(self) -> Quantity:
        """Width of the air-side face area."""
        width, _ = self.front_dims
        return width

    @property
    def height(self) -> Quantity:
        """Height of the air-side face area."""
        _, height = self.front_dims
        return height

    @property
    def dT_sh(self) -> Quantity:
        """Amount of refrigerant superheating."""
        return self.TD_in - self.superheat_offset

    def __str__(self):
        lines = [
            f"air in: {str(self.air_in)}",
            f"air out: {str(self.air_out)}",
            f"TD in: {self.TD_in.to('K'):~P.1f}",
            f"TD out: {self.TD_out.to('K'):~P.1f}",
            f"T evaporating: {self.T_evp.to('degC'):~P.1f}",
            f"degree of superheating: {self.dT_sh.to('K'):~P.1f}",
            f"aspect ratio air face area: 1:{1 / self.aspect_ratio:.0f}"
        ]
        if self.Q_dot is not None:
            lines_ext = [
                "---------------",
                f"heat transfer rate: {self.Q_dot.to('kW'):~P.3f}",
                f"air m-rate: {self.air_m_dot.to('kg / hr'):~P.3f}",
                f"air v-rate: {self.air_V_dot.to('m ** 3 / hr'):~P.3f}",
                f"air face velocity: {self.v_fa.to('m / s'):~P.1f}",
                f"air face width: {self.width.to('mm'):~P.0f}",
                f"air face height: {self.height.to('mm'):~P.0f}"
            ]
            lines.extend(lines_ext)
        return "\n".join(lines)
