from __future__ import annotations

import math
import shelve
import typing
from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
import numpy as np
from scipy.optimize import minimize
from hvac import Quantity

Q_ = Quantity


@dataclass
class Material:
    """
    Dataclass that groups material properties of a construction layer.

    Attributes
    ----------
    k: Quantity
        Conductivity
    rho: Quantity
        Mass density
    c: Quantity
        Specific heat capacity
    R: Quantity or None
        Unit thermal resistance
    """
    k: Quantity = Q_(0.0, 'W / (m * K)')
    rho: Quantity = Q_(0.0, 'kg / m**3')
    c: Quantity = Q_(0.0, 'J / (kg * K)')
    R: Quantity | None = None


@dataclass
class Geometry:
    """
    Dataclass that groups geometrical properties of a construction layer.

    Attributes
    ----------
    A: Quantity
        Area
    t: Quantity
        Thickness
    w: Quantity | None
        Width
    h: Quantity | None
        Height
    
    Notes
    -----
    If both width `w` and height `h` are specified (not None), area `A` will
    be the product of `w` and `h`. 
    """
    A: Quantity = Q_(1.0, 'm ** 2')
    t: Quantity = Q_(0.0, 'm')
    w: Quantity | None = None  # Q_(1.0, 'm')
    h: Quantity | None = None  # Q_(1.0, 'm')

    def __post_init__(self):
        if self.w is not None and self.h is not None:
            self.A = self.w * self.h


class HeatFlowDirection(Enum):
    """
    Enum class that defines the direction of heat flow through a building
    component. In case of a vertical oriented building component (a wall), the 
    heat flow will always be horizontal, and it doesn't matter if the sense of
    heat flow is from left to right or from right to left. In case of a 
    horizontal oriented building component (a ceiling or a roof) the heat flow 
    will always be vertical, but it needs to be specified if the heat flow has 
    an upwards or downwards sense. 
    """
    HORIZONTAL = 'horizontal'
    UPWARDS = 'upwards'
    DOWNWARDS = 'downwards'


class MechanicalFastening:
    """
    Class for determining the correction factor for the thermal transmittance
    of an insulation layer due to mechanical fasteners penetrating this 
    insulation layer according to standard ISO 6946:2017, Annex F.3.
    """
    def __init__(self):
        self._t_ins: Quantity = Q_(10, 'cm')
        self._alpha: Quantity = Q_(0.8, 'frac')
        self._kf: Quantity = Q_(52, 'W / (m * K)')  # steel
        self._Af: Quantity = Q_(0, 'mm ** 2')
        self._nf: Quantity = Q_(3, '1 / m ** 2')

    @classmethod
    def create(
        cls,
        diameter: Quantity = Q_(3, 'mm'),
        number_per_unit_area: Quantity = Q_(3, '1 / m ** 2'),
        insulation_thickness: Quantity = Q_(10, 'cm'),
        length: Quantity = Q_(10, 'cm'),
        conductivity: Quantity = Q_(52, 'W / (m * K)')
    ) -> MechanicalFastening:
        """
        Creates a `MechanicalFastening` object.

        Parameters
        ----------
        diameter: Quantity
            Diameter of a single fastener.
        number_per_unit_area: int
            Number of fasteners per unit area.
        insulation_thickness: Quantity
            Thickness of the insulation layer.
        length: Quantity
            Length of the fasteners that penetrate the insulation layer.
        conductivity: Quantity, default 52 W / (m.K) for steel
            Thermal conductivity of the fastener.

        Returns
        -------
        `MechanicalFastening` object. The correction factor can be accessed
        through property `correction_factor`.
        """
        fastener = cls()
        fastener._Af = math.pi * (diameter ** 2) / 4
        fastener._nf = number_per_unit_area
        fastener._alpha = 0.8 * length / insulation_thickness
        fastener._t_ins = insulation_thickness
        fastener._kf = conductivity
        return fastener

    @property
    def correction_factor(self) -> Quantity:
        """
        Returns the correction factor for thermal transmittance due to the 
        presence of mechanical fasteners in the insulation layer of the building
        element.
        """
        cf = self._alpha * self._kf * self._Af * self._nf / self._t_ins
        return cf.to('W / (m ** 2 * K)')


class ConstructionLayer:
    """
    Common base class that represents a flat layer in a construction assembly.

    An ID can identify a layer in a construction assembly. A construction layer
    has geometry (see class `Geometry`), thermal resistance (or thermal
    transmittance) and possibly also thermal capacity.

    To create a linear thermal network of an exterior building element, a layer
    can be further subdivided into a number of slices for better modeling of the
    true thermal inertia of the exterior building element.
    """
    def __init__(self):
        self.ID: str = ''
        self.geometry: Geometry = Geometry()
        self.R: Quantity | None = None
        self.C: Quantity = Q_(0.0, 'J / (K * m ** 2)')
        self.num_slices: int = 1

    @property
    def U(self) -> Quantity:
        """
        Returns the unit transmittance of the construction layer, which is
        the inverse of the unit thermal resistance.
        """
        return 1 / self.R

    @U.setter
    def U(self, v: Quantity) -> None:
        """Sets the unit transmittance of the construction layer."""
        self.R = 1 / v

    def __add__(self, other: ConstructionLayer) -> ConstructionLayer:
        """
        Creates a new `ConstructionLayer` object by adding two existing
        `ConstructionLayer` objects (`self` and `other`) in series.

        The thermal resistance, thermal capacity and thickness of both layers
        are added together.
        The resulting layer will have the same surface area as the construction
        layer on the left of the + operator (i.e., `self`).
        """
        layer = ConstructionLayer()
        layer.ID = f"{self.ID} + {other.ID}"
        A = self.geometry.A
        t = self.geometry.t + other.geometry.t
        layer.geometry = Geometry(A, t)
        layer.R = self.R + other.R
        layer.C = self.C + other.C
        layer.num_slices = self.num_slices + other.num_slices
        return layer

    def __radd__(self, other: ConstructionLayer | int) -> ConstructionLayer:
        """
        Creates a new `ConstructionLayer` object by adding two existing,
        stacked `ConstructionLayer` objects (`other` and `self`) in series.

        The thermal resistance, thermal capacity and thickness of both layers
        are added together.
        The resulting layer will have the same surface area as the construction
        layer on the right of the + operator (i.e., `self`).
        """
        if isinstance(other, int):
            return self
        else:
            return other.__add__(self)

    def __floordiv__(self, other: ConstructionLayer) -> ConstructionLayer:
        """
        Creates a new `ConstructionLayer` object by adding two existing,
        adjacent `ConstructionLayer` objects (`other` and `self`) in parallel.

        The unit thermal transmittance and the unit thermal capacity of the new
        layer are determined by the area-weighted sum of the unit thermal
        transmittance/unit thermal capacity of the two layers.
        """
        layer = ConstructionLayer()
        layer.ID = f"{self.ID} // {other.ID}"
        A = self.geometry.A + other.geometry.A
        t = self.geometry.t
        geometry = Geometry(A, t)
        layer.geometry = geometry
        U = (self.geometry.A / A) * self.U + (other.geometry.A / A) * other.U
        layer.R = 1 / U
        layer.C = (self.geometry.A / A) * self.C + (other.geometry.A / A) * other.C
        layer.num_slices = min(self.num_slices, other.num_slices)
        return layer

    def __str__(self):
        return (
            f"{str(self.__class__.__name__)} '{self.ID}': "
            f"t = {self.geometry.t.to('m'):~P.3f}, "
            f"R = {self.R.to('m ** 2 * K / W'):~P.2f}, "
            f"C = {self.C.to('J / (m ** 2 * K)'):~P.2f}, "
            f"nodes = {self.num_slices}"
        )


class SolidLayer(ConstructionLayer):
    """
    Represents a flat layer of solid building material in a construction
    assembly.
    """
    def __init__(self):
        super().__init__()
        self.material: Material | None = None

    @classmethod
    def create(
        cls,
        ID: str,
        geometry: Geometry,
        material: Material,
        num_slices: int = 1
    ) -> SolidLayer:
        """
        Creates a `SolidLayer` object.

        Parameters
        ----------
        ID:
            Name to identify the layer in the construction assembly.
        geometry:
            `Geometry` object that holds the surface area (`A`) and thickness
            (`t`) of the layer.
        material:
            `Material` object that holds the material properties of the layer
            (`k`, `rho`, `c`, `R`) where `k` is the thermal conductivity of
            the material, `rho` is the mass density, `c` is the specific heat,
            and `R` is the unit thermal resistance (default None).
        num_slices:
            A solid layer can be subdivided into a number of slices. Each slice
            will correspond with a temperature node in the linear thermal 
            network model of an exterior building element. 
        """
        layer = SolidLayer()
        layer.ID = ID
        layer.geometry = geometry
        layer.material = material
        layer.num_slices = num_slices
        layer.R, layer.C = layer._calc_R_and_C()
        return layer

    def _calc_R_and_C(self) -> tuple[Quantity, Quantity]:
        if self.material.R is None:
            try:
                self.material.R = self.geometry.t / self.material.k
            except ZeroDivisionError:
                self.material.R = Q_(float('inf'), 'K * m**2 / W')
        C = self.material.rho * self.material.c * self.geometry.t
        return self.material.R, C


class AirLayer(ConstructionLayer):
    """
    Represents an airspace between two building material layers in a
    construction assembly.

    The thermal capacity of the airspace is ignored (C = 0).

    The thermal resistance is calculated according to standard ISO 6946 (2017),
    annex D.2 for unventilated air spaces with length and width both more than
    10 times the thickness or annex D.4 for small or divided unventilated air
    spaces (air voids) with a width less than 10 times the thickness.
    """
    @classmethod
    def create(
        cls,
        ID: str,
        geometry: Geometry,
        dT: Quantity,
        heat_flow_dir: HeatFlowDirection,
        T_mn: Quantity,
        surf_emissivities: tuple[Quantity, Quantity] = (Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        angle: Quantity | None = None
    ) -> AirLayer:
        """
        Creates an `AirLayer` object.

        Parameters
        ----------
        ID: str
            Name to identify the airspace in a building component.
        geometry: Geometry
            Note: the maximum allowable thickness of the airspace is 0.3 m.
            The thickness is the dimension in the direction of heat flow.
        dT: Quantity
            The temperature difference across the airspace.
        heat_flow_dir: HeatFlowDirection
            Enum-type object indicating the direction of the heat flow
            perpendicular to the surface. Can be HORIZONTAL in case of a vertical
            surface, and UPWARDS or DOWNWARDS in case of a horizontal surface.
        T_mn: Quantity
            The mean thermodynamic temperature of the surfaces and the airspace.
        surf_emissivities: Tuple[Quantity, Quantity], (Q_(0.9, 'frac'), Q_(0.9, 'frac'))
            Emissivity of surfaces on both sides of the airspace.
        angle: Quantity, optional
            Inclination angle of the airspace with respect to the horizontal.
            Value can be between 0° and 90° included.
        """
        air_space = cls()
        air_space.ID = ID
        air_space.geometry = geometry
        if angle is not None:
            a = angle.to('deg').m
            ha90 = cls._determine_ha(
                geometry,
                dT,
                HeatFlowDirection.HORIZONTAL
            )
            ha0 = cls._determine_ha(
                geometry,
                dT,
                HeatFlowDirection.UPWARDS
            )
            if 0 < a < 90:
                ha = ha90 + (ha90 - ha0) * (a - 90) / 90
            elif a == 0.0:
                ha = ha0
            else:
                # a = 90.0
                ha = ha90
        else:
            ha = cls._determine_ha(
                geometry,
                dT,
                heat_flow_dir
            )
        hr = air_space._determine_hr(T_mn, surf_emissivities)
        air_space.R = 1 / (ha + hr)
        return air_space

    @staticmethod
    def _determine_ha(
        geometry: Geometry,
        dT: Quantity,
        heat_flow_dir: HeatFlowDirection
    ) -> Quantity:
        """Calculates the convection coefficient of the airspace."""
        t = geometry.t.to('m').m
        dT = dT.to('K').m
        ha = 0.0
        if dT <= 5:
            # ISO 6946 Table D.1
            if heat_flow_dir == HeatFlowDirection.HORIZONTAL:
                ha = max(1.25, 0.025 / t)
            elif heat_flow_dir == HeatFlowDirection.UPWARDS:
                ha = max(1.95, 0.025 / t)
            elif heat_flow_dir == HeatFlowDirection.DOWNWARDS:
                ha = max(0.12 * (t ** -0.44), 0.025 / t)
        if dT > 5:
            # ISO 6946 Table D.2
            if heat_flow_dir == HeatFlowDirection.HORIZONTAL:
                ha = max(0.73 * (dT ** (1 / 3)), 0.025 / t)
            elif heat_flow_dir == HeatFlowDirection.UPWARDS:
                ha = max(1.14 * (dT ** (1 / 3)), 0.025 / t)
            elif heat_flow_dir == HeatFlowDirection.DOWNWARDS:
                ha = max(0.09 * (dT ** 0.187) * (t ** -0.44), 0.025 / t)
        return Q_(ha, 'W / (m ** 2 * K)')

    def _determine_hr(
        self,
        T_mn: Quantity,
        surf_emissivities: tuple[Quantity, Quantity]
    ) -> Quantity:
        """Calculates the radiative coefficient of the airspace."""
        SIGMA = Q_(5.67e-8, 'W / (m ** 2 * K ** 4)')
        hro = 4 * SIGMA * T_mn.to('K') ** 3
        e1, e2 = surf_emissivities
        if self.geometry.w is not None and self.geometry.w <= 10.0 * self.geometry.t:
            # air void
            r = self.geometry.t / self.geometry.w
            d = 1 + math.sqrt(1 + r ** 2) - r
            hr = hro / ((1 / e1) + (1 / e2) - 2 + 2 / d)
        else:
            hr = hro / ((1 / e1) + (1 / e2) - 1)
        return hr


class AirLayerTemperatureSolver:
    """Solves for the temperatures on both sides of an air layer, when the
    following conditions are known:
    - the exterior temperature
    - the unit thermal resistance between the exterior and the outer air layer
      side
    - the unit thermal resistance between the inner air layer side and the
      interior zone air
    - the zone air temperature.

    The solving method searches for the temperatures on the outer and the inner
    side of the air layer, so that the heat flow from the exterior to the
    outer side of the air layer balances with the heat flow from the inner side
    of the air layer to the interior zone.
    """

    def __init__(
        self,
        T_ext: Quantity,
        T_int: Quantity,
        R_ea: Quantity,
        R_ai: Quantity
    ) -> None:
        """Creates an `AirLayerTemperatureSolver` instance.

        Parameters
        ----------
        T_ext:
            Exterior temperature.
        T_int:
            Interior temperature.
        R_ea:
            Unit thermal resistance between exterior and outer air layer side.
        R_ai:
            Unit thermal resistance between inner air layer side and interior.
        """
        self.T_ext = T_ext.to('K').m
        self.T_int = T_int.to('K').m
        self.R_ea = R_ea.to('K * m**2 / W').m
        self.R_ai = R_ai.to('K * m**2 / W').m

    def _eq1(self, T_ae: float | np.ndarray) -> float:
        """Heat flow from the exterior to the outer side of the air layer."""
        return (self.T_ext - T_ae) / self.R_ea

    def _eq2(self, T_ai: float | np.ndarray) -> float:
        """Heat flow from the inner side of the air layer to the interior zone."""
        return (T_ai - self.T_int) / self.R_ai

    def _sys_eqs(self, unknowns: np.ndarray) -> float:
        T_ae = unknowns[0]
        T_ai = unknowns[1]
        Q_dot_ea = self._eq1(T_ae)
        Q_dot_az = self._eq2(T_ai)
        dev_Q_dot = abs(Q_dot_ea - Q_dot_az)
        return dev_Q_dot

    def solve(
        self,
        T_ae_guess: Quantity,
        T_ai_guess: Quantity
    ) -> tuple[Quantity, ...]:
        """Solves for the temperatures on the outer and the inner side of the
        air layer.

        Parameters
        ----------
        T_ae_guess:
            Initial guess for the temperature on the outer side of the air layer.
        T_ai_guess:
            Initial guess for the temperature on the inner side of the air layer.

        Returns
        -------
        1. temperature on the outer side of the air layer
        2. temperature on the inner side of the air layer
        3. temperature difference between outer and inner side of the air layer
        4. average temperature in the air layer
        """
        ini_T_ae = T_ae_guess.to('K').m
        ini_T_ai = T_ai_guess.to('K').m
        sol = minimize(self._sys_eqs, np.array([ini_T_ae, ini_T_ai]))
        T_ae = Q_(sol.x[0], 'K')
        T_ai = Q_(sol.x[1], 'K')
        dT = abs(T_ae - T_ai)
        T_avg = (T_ae + T_ai) / 2
        return T_ae, T_ai, dT, T_avg


class SurfaceFilm(ConstructionLayer):
    """
    Represents the internal or external surface film of a construction 
    assembly. 
    
    The thermal resistance is calculated according to standard ISO 6946 (2017), 
    annex C.
    """
    @classmethod
    def create(
        cls,
        ID: str,
        geometry: Geometry,
        heat_flow_dir: HeatFlowDirection,
        T_mn: Quantity = Q_(20, 'degC'),
        surf_emissivity: Quantity = Q_(0.9, 'frac'),
        is_internal_surf: bool = True,
        wind_speed: Quantity = Q_(4.0, 'm / s')
    ) -> SurfaceFilm:
        """
        Creates a `SurfaceFilm` object.

        Parameters
        ----------
        ID: str
            Name to identify the surface film.
        geometry: Geometry
            `Geometry` object with the area (A) of the surface film.
            Note: in case of a surface film, the film thickness is not 
            considered.
        heat_flow_dir: HeatFlowDirection
            `HeatFlowDirection` object indicating the direction of the heat flow
            perpendicular to the surface. Can be HORIZONTAL in case of a vertical
            surface, and UPWARDS or DOWNWARDS in case of a horizontal surface.
        T_mn: Quantity
            The mean thermodynamic temperature of the surface and of its
            surroundings (e.g. operative temperature for internal surfaces,
            external air temperature or sol-air temperature for external 
            surfaces).
        surf_emissivity: Quantity, default 0.9
            The hemispherical emissivity of the surface.
        is_internal_surf: bool, default True
            Indicates whether the surface is internal or external.
        wind_speed: Quantity, default 4 m/s
            Wind speed. Only used with external surfaces.
        """
        surf_layer = cls()
        surf_layer.ID = ID
        surf_layer.geometry = geometry
        hr = cls._determine_hr(T_mn, surf_emissivity)
        hc = cls._determine_hc(heat_flow_dir, is_internal_surf, wind_speed)
        surf_layer.R = 1 / (hr + hc)
        return surf_layer

    @staticmethod
    def _determine_hr(T_mn: Quantity, e: Quantity) -> Quantity:
        """Calculates the radiative coefficient of the surface film."""
        SIGMA = Q_(5.67e-8, 'W / (m ** 2 * K ** 4)')
        hro = 4 * SIGMA * T_mn.to('K') ** 3
        hr = e * hro
        return hr

    @staticmethod
    def _determine_hc(
        heat_flow_dir: HeatFlowDirection, 
        is_internal_surf: bool, 
        v_wind: Quantity
    ) -> Quantity:
        """Calculates the convection coefficient of the surface film."""
        if is_internal_surf:
            if heat_flow_dir == HeatFlowDirection.UPWARDS:
                hc = Q_(5.0, 'W / (m ** 2 * K)')
            elif heat_flow_dir == HeatFlowDirection.HORIZONTAL:
                hc = Q_(2.5, 'W / (m ** 2 * K)')
            else:
                hc = Q_(0.7, 'W / (m ** 2 * K)')
        else:
            hc = Q_(4 + 4 * v_wind.to('m / s').m, 'W / (m ** 2 * K)')
        return hc


def apply_insulation_correction(
    insulation_layer: SolidLayer,
    building_component: ConstructionLayer,
    installation_quality_level: int = 0,
    mechanical_fastening: MechanicalFastening | None = None
) -> Quantity:
    """
    Apply correction to thermal resistance/transmittance to allow for the effects
    of air voids in insulation and mechanical fasteners penetrating the insulation
    layer according to ISO 6946:2017, Annex F.

    Parameters
    ----------
    insulation_layer: ConstructionLayer
        The insulation layer within the building component or element.
    building_component: ConstructionLayer
        The complete building component or element that contains the insulation
        layer.
    installation_quality_level: int, [0 (default), 1, or 2]
        Installation quality level of the insulation according to ISO 6946:2017,
        table F.1:
            - level 0: no air voids
            - level 1: air voids, but no air circulation
            - level 2: air voids with free air circulation between the warm and 
                       cold side of the insulation.
    mechanical_fastening: MechanicalFastening, optional
        Object that describes the mechanical fastening of the insulation layer.
        See `MechanicalFastening` class.

    Returns
    -------
    The correction term `delta_U` to be added to the thermal transmittance U of
    the building component or element.
    """
    correction_factor = Q_(0.0, 'W / (m ** 2 * K)')
    # correction for air voids
    if installation_quality_level == 1:
        correction_factor += Q_(0.01, 'W / (m ** 2 * K)')
    if installation_quality_level == 2:
        correction_factor += Q_(0.04, 'W / (m ** 2 * K)')
    # correction for mechanical fasteners
    if isinstance(mechanical_fastening, MechanicalFastening):
        correction_factor += mechanical_fastening.correction_factor
    # calculate delta_U
    R1 = insulation_layer.R.to('m ** 2 * K / W')
    Rtot = building_component.R.to('m ** 2 * K / W')
    delta_U = correction_factor * (R1 / Rtot) ** 2
    return delta_U


class ConstructionAssembly:
    """
    Represents a flat construction assembly composed of `ConstructionLayer` 
    objects.
    """
    db_path: str

    def __init__(self):
        self.ID: str = ''
        self.layers: dict[str, ConstructionLayer] | None = None
        self._building_component: ConstructionLayer = ConstructionLayer()
        # `_building_component` is a `ConstructionLayer` object that contains
        # the total thermal resistance of the construction assembly.

    @classmethod
    def create(
        cls,
        ID: str,
        layers: list[ConstructionLayer] | None = None,
        U: Quantity | None = None,
        R: Quantity | None = None,
        geometry: Geometry | None = None
    ) -> ConstructionAssembly:
        """
        Creates a `ConstructionAssembly` object.

        Parameters
        ----------
        ID: str
            Name to identify the construction assembly.
        layers: list of `ConstructionLayer` objects, optional, default None
            A construction assembly is thought to be composed of successive
            construction layers, arranged from the outdoor to the indoor
            environment.
            Note that it is expected that a construction assembly always has an
            exterior and an interior surface film.
        U: Quantity, default None
            U-value of the construction assembly per unit area. Can be used if
            the layered composition of the construction assembly is unknown.
        R: Quantity, default None
            R-value of the construction assembly per unit area. Can be used if
            the layered composition of the construction assembly is unknown.
        geometry: Geometry, default None
            The geometrical properties of the construction assembly (area and
            thickness). Use this if `layers` is None, and `U` or `R` are used
            instead.
        """
        constr_assem = cls()
        constr_assem.ID = ID
        if layers is not None:
            constr_assem.layers = {layer.ID: layer for layer in layers}
            # noinspection PyTypeChecker
            constr_assem._building_component = sum(layers)
        else:
            if geometry is not None:
                constr_assem._building_component.geometry = geometry
            if U is not None:
                constr_assem._building_component.U = U
            elif R is not None:
                constr_assem._building_component.R = R
            else:
                raise TypeError(
                    "Either `U`or `R` must be specified when instantiating"
                    "a construction assembly without layers."
                )
        return constr_assem

    def apply_insulation_correction(
        self,
        insulation_layer_ID: str,
        insulation_level: int = 0,
        mechanical_fastening: MechanicalFastening | None = None
    ) -> ConstructionAssembly:
        """
        Applies a correction to the thermal resistance and transmittance of
        the construction assembly taking air voids and/or mechanical fasteners
        in the insulation layer into account.

        Parameters
        ----------
        insulation_layer_ID: str
            The ID of the insulation layer in the construction assembly.
        insulation_level: int, [0 (default), 1, 2]
            Insulation level according to ISO 6946:2017, Annex F, table F.1:
            - Level 0: No air voids within the insulation, or where only minor
            air voids are present that have no significant effect on the thermal
            transmittance.
            - Level 1: Air gaps bridging between the hot and cold side of the
            insulation, but not causing air circulation between the warm and
            cold side of the insulation.
            - Level 2: Air gaps bridging between the hot and cold side of the
            insulation, combined with cavities resulting in free air circulation
            between the warm and cold sides of the insulation.
        mechanical_fastening: MechanicalFastening, optional
            Object that describes the mechanical fastening of the insulation layer.
            See `MechanicalFastening` class.

        Returns
        -------
        A new `ConstructionAssembly` object with the applied correction.
        """
        new_assembly = deepcopy(self)
        insulation_layer = new_assembly.layers.get(insulation_layer_ID)
        if insulation_layer:
            # noinspection PyTypeChecker
            delta_U = apply_insulation_correction(
                insulation_layer, new_assembly._building_component,
                insulation_level, mechanical_fastening
            )
            insulation_layer.U += delta_U
            # Also update total R of the construction assembly:
            # noinspection PyTypeChecker
            new_assembly._building_component = sum(new_assembly.layers.values())
        return new_assembly

    def apply_ventilated_layer_correction(
        self,
        ventilated_layer_ID: str,
        area_of_openings: Quantity,
        heat_flow_dir: HeatFlowDirection,
        T_mn: Quantity,
        surf_emissivity: Quantity = Q_(0.9, 'frac'),
    ) -> ConstructionAssembly:
        """
        Applies a correction to the thermal resistance and transmittance of
        the construction assembly due to a slightly or well-ventilated air layer
        according to ISO 6946:2017 6.9.3 and 6.9.4.

        Parameters
        ----------
        ventilated_layer_ID: str
            ID of the ventilated air layer.
        area_of_openings: Quantity
            A slightly ventilated air layer is one in which there is provision
            for limited air flow through it from the external environment by
            openings of area, A_ve, within the following ranges:
            — >500 mm2 but <1 500 mm2 per metre of length (in the horizontal
            direction) for vertical air layers;
            — >500 mm2 but <1 500 mm2 per square metre of surface area for
            horizontal air layers.
            A well-ventilated air layer is one for which the openings between
            the air layer and the external environment are equal to or exceed
            — 1 500 mm2 per metre of length (in the horizontal direction) for
            vertical air layers, and
            — 1 500 mm2 per square of metre of surface area for horizontal air
            layers.
        heat_flow_dir: HeatFlowDirection
            `HeatFlowDirection` object indicating the direction of the heat flow
            perpendicular to the surface. Can be HORIZONTAL in case of a vertical
            surface, and UPWARDS or DOWNWARDS in case of a horizontal surface.
        T_mn: Quantity
            Mean thermodynamic temperature of airspace.
        surf_emissivity: Quantity
            The hemispherical emissivity of the surface adjacent to the
            ventilated air layer.

        Returns
        -------
        A new `ConstructionAssembly` object with the correction for a slightly
        or well ventilated air layer (if needed, else the original
        construction assembly is returned).
        """
        A_ve = area_of_openings.to('mm ** 2').m
        if A_ve > 500:
            # Get layers between the interior surface and the ventilated air layer:
            i = list(self.layers.keys()).index(ventilated_layer_ID) + 1
            interior_layers = list(self.layers.values())[i:]

            # Sum the interior layers (returns a `ConstructionLayer` object):
            # noinspection PyTypeChecker
            interior_part: ConstructionLayer = sum(interior_layers)

            # Add external surface layer for still air to the internal part:
            ext_surf_layer = SurfaceFilm.create(
                ID='ext_surf_film',
                geometry=Geometry(),
                heat_flow_dir=heat_flow_dir,
                T_mn=T_mn,
                surf_emissivity=surf_emissivity,
                is_internal_surf=False,
                wind_speed=Q_(0, 'm / s')
            )

            # Get the total thermal resistance of building component with
            # well-ventilated air layer:
            R_ve = (interior_part.R + ext_surf_layer.R).to('m ** 2 * K / W')

            # Get the total thermal resistance of building component with unventilated
            # air layer
            R_nve = self._building_component.R.to('m ** 2 * K / W')

            # Slightly ventilated air layer:
            if 500 < A_ve < 1500:
                Rtot = (1500 - A_ve) / 1000 * R_nve + (A_ve - 500) / 1000 * R_ve
                new_layers = deepcopy(list(self.layers.values()))
                new_assembly = ConstructionAssembly.create(
                    ID=self.ID,
                    layers=new_layers
                )
                air_layer = new_assembly.layers[ventilated_layer_ID]
                air_layer.R = self.R - Rtot
                new_assembly._building_component.R = Rtot
                return new_assembly

            # Well-ventilated air layer:
            if A_ve >= 1500:
                new_layers = [ext_surf_layer]
                new_layers.extend(interior_layers)
                new_assembly = ConstructionAssembly.create(
                    ID=self.ID,
                    layers=new_layers
                )
                new_assembly._building_component.R = R_ve
                return new_assembly
        return self

    @property
    def R(self) -> Quantity:
        """Returns the unit thermal resistance of the construction assembly."""
        return self._building_component.R

    @property
    def U(self) -> Quantity:
        """
        Returns the unit thermal transmittance of the construction
        assembly.
        """
        return 1 / self._building_component.R

    @property
    def R_surf_ext(self) -> Quantity | None:
        """
        Returns the surface resistance at the exterior side of the
        construction assembly.
        """
        if self.layers is not None:
            exterior_surface_layer = list(self.layers.values())[0]
            return exterior_surface_layer.R
        return None

    @property
    def thickness(self) -> Quantity:
        """Returns the thickness of the construction assembly."""
        if self.layers is not None:
            t = sum(layer.geometry.t for layer in self.layers.values())
            return t
        else: 
            return self._building_component.geometry.t 

    @property
    def area(self) -> Quantity:
        """Returns the surface area of the construction assembly."""
        return self._building_component.geometry.A

    @area.setter
    def area(self, v: Quantity) -> None:
        """Sets the surface area of the construction assembly."""
        self._building_component.geometry.A = v

    def __str__(self):
        _str = f'construction assembly: {self.ID}\n'
        for layer_str in (str(layer) for layer in self.layers.values()):
            _str += layer_str + '\n'
        return _str

    def save(self) -> None:
        """
        Saves the construction assembly to a shelf on disk.

        Class variable `db_path` must be assigned the path string to this shelf,
        otherwise a `ValueError` will be raised.
        """
        if self.db_path:
            with shelve.open(self.db_path) as shelf:
                shelf[self.ID] = self
        else:
            raise ValueError('Path to shelf is not specified.')

    @classmethod
    def load(cls, ID: str) -> ConstructionAssembly:
        """
        Loads the construction assembly with the given ID from a shelf on
        disk.

        Class variable `db_path` must be assigned the path string to this shelf,
        otherwise a `ValueError` will be raised.
        """
        if cls.db_path:
            with shelve.open(cls.db_path) as shelf:
                constr_assem = typing.cast('ConstructionAssembly', shelf[ID])
                return constr_assem
        else:
            raise ValueError('Path to shelf is not specified.')

    @classmethod
    def overview(cls) -> list[str]:
        """
        Returns a list with all the ID's of the construction assemblies
        stored in the shelf.

        Class variable `db_path` must be assigned the path string to this shelf,
        otherwise a `ValueError` will be raised.
        """
        if cls.db_path:
            with shelve.open(cls.db_path) as shelf:
                return list(shelf.keys())
        else:
            raise ValueError('Path to shelf is not specified.')
    
    @classmethod
    def remove(cls, ID: str) -> None:
        """
        Removes the construction assembly identified by `ID` from the shelf.
        
        Parameters
        ----------
        ID:
            ID of the construction assembly to be removed.

        Returns
        -------
        None
        
        Raises
        ------
        ValueError:
            If the class variable `db_path` is not set.
        KeyError:
            If the construction assembly with the given ID is not on the shelf.
        """
        if cls.db_path:
            with shelve.open(cls.db_path) as shelf:
                try:
                    del shelf[ID]
                except KeyError:
                    raise KeyError(
                        f"Construction assembly {ID} is not on the shelf."
                    )
        else:
            raise ValueError('Path to shelf is not specified.')
    
    def add_layers(self, new_layers: list[tuple[int, ConstructionLayer]]) -> None:
        """
        Adds one or more new layers to the construction assembly.

        Parameters
        ----------
        new_layers:
            List of tuples with 2 elements. The first element is the index
            of the new layer in the layered composition of the construction
            assembly. The second element is the new `ConstructionLayer` object.
            Keep in mind that layers must be ordered from the exterior side
            towards the interior side.
        """
        layers = list(self.layers.values())
        for pos, new_layer in new_layers:
            layers.insert(pos, new_layer)
        self.layers = {layer.ID: layer for layer in layers}
        self._building_component += sum(layer for _, layer in new_layers)
