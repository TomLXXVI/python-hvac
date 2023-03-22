from typing import Callable
import pandas as pd
from hvac import Quantity
from hvac.climate import ClimateData
from hvac.climate.sun.solar_time import time_from_decimal_hour
from .construction_assembly import (
    ThermalComponent,
    ConstructionAssembly
)
from .thermal_network import (
    AbstractNode,
    ExteriorSurfaceNode,
    BuildingMassNode,
    InteriorSurfaceNode,
    ThermalNetwork,
    ThermalNetworkSolver
)
from .exterior_surface import ExteriorSurface
from .fenestration import (
    WindowThermalProperties,
    ExteriorShadingDevice,
    InteriorShadingDevice,
    Window
)


Q_ = Quantity


class ThermalNetworkBuilder:
    """Builds the linear thermal network of a construction assembly."""

    @classmethod
    def build(cls, construction_assembly: ConstructionAssembly) -> ThermalNetwork | None:
        """
        Get the linear thermal network of the construction assembly.

        Returns
        -------
        `ThermalNetwork` object.

        Notes
        -----
        The layers of the construction assembly must be arranged from
        the exterior surface towards the interior surface.
        """
        thermal_network = cls._compose(list(construction_assembly.layers.values()))
        reduced_thermal_network = cls._reduce(thermal_network)
        nodes = cls._transform(reduced_thermal_network)
        return ThermalNetwork(nodes)

    @staticmethod
    def _compose(layers: list[ThermalComponent]) -> list[Quantity]:
        # create a list of resistors and capacitors
        thermal_network = []
        for layer in layers:
            n = layer.slices
            R_slice = layer.R / (2 * n)
            C_slice = layer.C / n
            slices = [R_slice, C_slice, R_slice] * n
            thermal_network.extend(slices)
        return thermal_network

    @staticmethod
    def _reduce(thermal_network: list[Quantity]) -> list[Quantity]:
        # sum adjacent resistors between capacitors
        R_dummy = Q_(0, 'm ** 2 * K / W')
        reduced_thermal_network = []
        R = thermal_network[0]
        for i in range(1, len(thermal_network)):
            if R_dummy.check(thermal_network[i].dimensionality):
                # thermal_network[i] is a resistance
                R += thermal_network[i]
            else:
                # thermal_network[i] is a capacitance: only keep C
                # if it is > 0, unless for the last C (to set the interior
                # surface node, see _transform)
                if thermal_network[i].m > 0 or i == len(thermal_network) - 2:
                    reduced_thermal_network.append(R)
                    reduced_thermal_network.append(thermal_network[i])
                    R = Q_(0.0, 'm ** 2 * K / W')
        if R.m > 0:
            reduced_thermal_network.append(R)
        return reduced_thermal_network

    @staticmethod
    def _transform(reduced_thermal_network: list[Quantity]) -> list[AbstractNode] | None:
        # create list of nodes, starting at the exterior surface towards the interior surface
        if len(reduced_thermal_network) >= 5:
            i = 1
            node_index = 1
            nodes = []
            while True:
                if i == 1:
                    node = ExteriorSurfaceNode(
                        ID=f'N{node_index}',
                        R_unit_lst=[
                            reduced_thermal_network[i - 1],
                            reduced_thermal_network[i + 1]
                        ],
                        C=reduced_thermal_network[i]
                    )
                    nodes.append(node)
                    i += 2
                    node_index += 1
                elif i == len(reduced_thermal_network) - 2:
                    node = InteriorSurfaceNode(
                        ID=f'N{node_index}',
                        R_unit_lst=[
                            reduced_thermal_network[i - 1],
                            reduced_thermal_network[i + 1]
                        ]
                    )
                    nodes.append(node)
                    break
                else:
                    node = BuildingMassNode(
                        ID=f'N{node_index}',
                        R_unit_lst=[
                            reduced_thermal_network[i - 1],
                            reduced_thermal_network[i + 1]
                        ],
                        C=reduced_thermal_network[i]
                    )
                    nodes.append(node)
                    i += 2
                    node_index += 1
            return nodes
        return None


class ExteriorBuildingElement:
    """
    An `ExteriorBuildingElement` is composed of:
    - an `ExteriorSurface` that defines the orientation and size of the
    building element and that is responsible for calculating the sol-air
    temperature under clear-sky conditions (based on `ClimateData`).
    - a `ConstructionAssembly` that defines the layered material composition of
    the building element and knows about the thermal resistance and thermal
    capacity of the building element.
    - a `ThermalNetwork` model of the construction assembly that is used to
    calculate the dynamic and steady heat transfer through the building element.
    """

    def __init__(self):
        self.ID: str = ''
        self._exterior_surface: ExteriorSurface | None = None
        self.construction_assembly: ConstructionAssembly | None = None
        self._thermal_network: ThermalNetwork | None = None
        self.T_int_fun: Callable[[float], float] | None = None
        self.F_rad: float | None = None
        self._heat_transfer: dict[str, list[Quantity]] | None = None
        self.windows: dict[str, Window] = {}
        self.doors: dict[str, 'ExteriorBuildingElement'] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        azimuth: Quantity,
        tilt: Quantity,
        width: Quantity,
        height: Quantity,
        climate_data: ClimateData,
        construction_assembly: ConstructionAssembly,
        T_int_fun: Callable[[float], float],
        F_rad: float = 0.46,
        surface_absorptance: Quantity | None = None,
        surface_color: str = 'dark-colored'
    ) -> 'ExteriorBuildingElement':
        """
        Create an `ExteriorBuildingElement` object.

        Parameters
        ----------
        ID: str
            Identifier for the building element.
        azimuth: Quantity
            The azimuth angle of the exterior building element measured clockwise
            from North (East = 90°, South = 180°, West = 270 °, North = 0°)
        tilt: Quantity
            The tilt angle of the exterior surface with respect to the horizontal
            surface. A vertical wall has a tilt angle of 90°.
        width: Quantity
            The width of the exterior surface.
        height: Quantity
            The height of the exterior surface.
        climate_data: ClimateData
            Object that determines the outdoor dry- and wet-bulb temperature and
            solar radiation for each hour of the design day. See class
            ClimateData.
        construction_assembly: ConstructionAssembly
            Object that defines the construction of the building element. See
            class ConstructionAssembly.
        T_int_fun: Callable
            Function that returns the indoor air temperature at time t expressed
            in seconds from 00:00:00. May also be an instance of `TemperatureSchedule`.
        F_rad: float, default 0.46
            Fraction of conductive heat flow that is transferred by radiation
            to the internal mass of the space (ASHRAE Fundamentals 2017, Ch. 18,
            Table 14).
        surface_absorptance: Quantity | None (default)
            The absorptance of the exterior surface of the building element.
        surface_color: ['dark-colored' (default), 'light-colored']
            Indicate if the surface is either dark-colored or light-colored. You
            can use this instead of specifying `surface_absorptance`.

        Returns
        -------
        ExteriorBuildingElement
        """
        ext_building_element = cls()
        ext_building_element.ID = ID
        ext_building_element.construction_assembly = construction_assembly
        ext_building_element._exterior_surface = ExteriorSurface(
            azimuth=azimuth,
            tilt=tilt,
            width=width,
            height=height,
            climate_data=climate_data,
            surface_resistance=construction_assembly.R_surf_ext,
            surface_absorptance=surface_absorptance,
            surface_color=surface_color
        )
        ext_building_element.T_int_fun = T_int_fun
        ext_building_element.F_rad = F_rad
        return ext_building_element

    @property
    def thermal_network(self) -> ThermalNetwork:
        """
        Get thermal network.
        """
        if self._thermal_network is None:
            self._thermal_network = ThermalNetworkBuilder.build(self.construction_assembly)
            self._thermal_network.T_ext = self._exterior_surface.T_sol
            self._thermal_network.T_int = self.T_int_fun
            for node in self._thermal_network.nodes: node.A = self.area_net
        return self._thermal_network

    @property
    def area_net(self) -> Quantity:
        A_net = self._exterior_surface.area
        if self.windows:
            A_wnd = sum(window.area for window in self.windows.values())
            A_net -= A_wnd
        if self.doors:
            A_drs = sum(door.area_gross for door in self.doors.values())
            A_net -= A_drs
        return A_net

    @property
    def area_gross(self) -> Quantity:
        return self._exterior_surface.area

    def irr_profile(self, unit_theta_i: str = 'deg', unit_I: str = 'W / m ** 2') -> pd.DataFrame:
        d = self._exterior_surface.irr_profile
        time = [dt.time() for dt in d['t']]
        theta_i = [theta_i.to(unit_theta_i).m for theta_i in d['theta_i']]
        glo_sur = [glo_sur.to(unit_I).m for glo_sur in d['glo_sur']]
        dir_sur = [dir_sur.to(unit_I).m for dir_sur in d['dir_sur']]
        dif_sur = [dif_sur.to(unit_I).m for dif_sur in d['dif_sur']]
        d = {
            'theta_i': theta_i,
            'glo_sur': glo_sur,
            'dir_sur': dir_sur,
            'dif_sur': dif_sur
        }
        df = pd.DataFrame(data=d, index=time)
        return df

    def temp_profile(self, unit_T: str = 'degC') -> pd.DataFrame:
        d_T_sol = self._exterior_surface.T_sol_profile
        d_T_db = self._exterior_surface.climate_data.Tdb_profile
        d_T_wb = self._exterior_surface.climate_data.Twb_profile
        time = [dt.time() for dt in d_T_sol['t']]
        T_sol = [T_sol.to(unit_T).m for T_sol in d_T_sol['T']]
        T_db = [T_db.to(unit_T).m for T_db in d_T_db['T']]
        T_wb = [T_wb.to(unit_T).m for T_wb in d_T_wb['T']]
        d = {
            'T_db': T_db,
            'T_wb': T_wb,
            'T_sol': T_sol
        }
        df = pd.DataFrame(data=d, index=time)
        return df

    def add_window(
        self,
        ID: str,
        width: Quantity,
        height: Quantity,
        therm_props: WindowThermalProperties,
        F_rad: float = 0.46,
        ext_shading_dev: ExteriorShadingDevice | None = None,
        int_shading_dev: InteriorShadingDevice | None = None
    ) -> None:
        """
        Add a window to the exterior building element.

        Parameters
        ----------
        ID:
            Identifier for the window.
        width:
            The width of the window.
        height:
            The height (or length) of the window.
        therm_props:
            Instance of class `WindowThermalProperties` containing the U-value
            of the entire window, the center-of-glass SHGCs for direct solar
            radiation at different solar incidence angles, the center-of-glass
            SHGC for diffuse solar radiation, and an overall SHGC for the
            entire window at normal incidence.
        F_rad: default 0.46
            The fraction of conductive and diffuse solar heat gain that is
            transferred by radiation to the interior thermal mass of the space.
            (see ASHRAE Fundamentals 2017, chapter 18, table 14).
        ext_shading_dev: default None
            See class `ExternalShadingDevice`. E.g. overhang or recessed window.
            In case a window is equipped with an external shading device, or in
            case of a recessed window, part of the window may be shaded depending
            on the position of the sun during the course of day.
        int_shading_dev: default None
            See class `InteriorShadingDevice`. E.g. louvered shades, roller
            shades, draperies, insect screens.
        """
        window = Window.create(
            ID=ID,
            azimuth=self._exterior_surface.azimuth,
            tilt=self._exterior_surface.tilt,
            width=width,
            height=height,
            climate_data=self._exterior_surface.climate_data,
            therm_props=therm_props,
            F_rad=F_rad,
            ext_shading_dev=ext_shading_dev,
            int_shading_dev=int_shading_dev
        )
        self.windows[window.ID] = window

    def add_door(
        self,
        ID: str,
        width: Quantity,
        height: Quantity,
        construction_assembly: ConstructionAssembly,
        F_rad: float = 0.46,
        surface_absorptance: Quantity | None = None,
        surface_color: str = 'dark-colored'
    ) -> None:
        """
        Add door to exterior building element. An exterior door is also
        regarded as an exterior building element. See class method `create(...)`
        for more explanation about the parameters.
        """
        door = ExteriorBuildingElement.create(
            ID=ID,
            azimuth=self._exterior_surface.azimuth,
            tilt=self._exterior_surface.tilt,
            width=width,
            height=height,
            climate_data=self._exterior_surface.climate_data,
            construction_assembly=construction_assembly,
            T_int_fun=self.T_int_fun,
            F_rad=F_rad,
            surface_absorptance=surface_absorptance,
            surface_color=surface_color
        )
        self.doors[door.ID] = door

    @property
    def UA(self) -> float:
        # for internal use only
        U = self.construction_assembly.U.to('W / (m ** 2 * K)').m
        A = self.area_net.to('m ** 2').m
        return U * A

    def T_sol(self, t: float) -> float:
        # for internal use only
        return self._exterior_surface.T_sol(t)

    def get_conductive_heat_gain(self, t: float, T_int: float) -> dict[str, float]:
        # for internal use only
        # note: the **steady-state** conductive heat gain is calculated.
        U = self.construction_assembly.U.to('W / (m ** 2 * K)').m
        A = self.area_net.to('m ** 2').m
        T_sol = self._exterior_surface.T_sol(t)
        Q = U * A * (T_sol - T_int)
        Q_rad = self.F_rad * Q
        Q_conv = Q - Q_rad
        return {'rad': Q_rad, 'conv': Q_conv}

    def get_heat_transfer(
        self,
        dt_hr: float = 1.0,
        n_cycles: int = 6,
        unit: str = 'W'
    ) -> pd.DataFrame:
        """
        Get the diurnal heat transfer cycle through the exterior building
        element.

        Parameters
        ----------
        dt_hr:
            The time step of the calculations in decimal hours.
        n_cycles:
            The number of repeated cycles before returning the result.
        unit: default 'W'
            The desired unit in which thermal power is to be expressed.

        Returns
        -------
        A Pandas DataFrame with the dynamic heat flows at the exterior side (key
        `Q_ext`) and interior side (key `Q_int`) of the building element, and
        also the calculated steady-state value (based on the temperature
        difference between exterior, sol-air temperature and indoor temperature
        and the thermal resistance of the building element) (key `Q_steady`).
        """
        if self._heat_transfer is None:
            tnw_solved = ThermalNetworkSolver.solve(
                (self.thermal_network, None),
                None,
                dt_hr,
                n_cycles
            )
            dt = dt_hr * 3600
            self._heat_transfer = {
                'Q_ext': [],
                'Q_int': [],
                'Q_steady': []
            }
            t_ax = []
            for k, T_node_list in enumerate(tnw_solved.T_node_table):
                t_ax.append(time_from_decimal_hour(k * dt_hr))
                T_sol = Q_(self._exterior_surface.T_sol(k * dt), 'degC').to('K')
                T_esn = T_node_list[0].to('K')
                R_ext = tnw_solved.R_ext
                Q_ext = (T_sol - T_esn) / R_ext
                T_int = Q_(self.T_int_fun(k * dt), 'degC').to('K')
                T_node_int = T_node_list[-1].to('K')
                R_int = tnw_solved.R_int
                Q_int = (T_node_int - T_int) / R_int
                R_tot = self.construction_assembly.R.to('m ** 2 * K / W')
                A = self.area_net.to('m ** 2')
                q_steady = (T_sol - T_int) / R_tot
                Q_steady = q_steady * A
                self._heat_transfer['Q_ext'].append(Q_ext.to(unit).m)
                self._heat_transfer['Q_int'].append(Q_int.to(unit).m)
                self._heat_transfer['Q_steady'].append(Q_steady.to(unit).m)
            self._heat_transfer = pd.DataFrame(data=self._heat_transfer, index=t_ax)
        return self._heat_transfer


class InteriorBuildingElement:

    def __init__(self):
        self.ID: str = ''
        self.width: Quantity | None = None
        self.height: Quantity | None = None
        self.construction_assembly: Quantity | None = None
        self.T_adj_fun: Callable[[float], float] | None = None
        self.F_rad: float = 0.46
        self.doors: dict[str, 'InteriorBuildingElement'] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        width: Quantity,
        height: Quantity,
        construction_assembly: Quantity,
        T_adj_fun: Callable[[float], float],
        F_rad: float = 0.46
    ) -> 'InteriorBuildingElement':
        int_build_elem = cls()
        int_build_elem.ID = ID
        int_build_elem.width = width
        int_build_elem.height = height
        int_build_elem.construction_assembly = construction_assembly
        int_build_elem.T_adj_fun = T_adj_fun
        int_build_elem.F_rad = F_rad
        return int_build_elem

    @property
    def area_net(self) -> Quantity:
        A_net = self.area_gross
        if self.doors:
            A_drs = sum(door.area_gross for door in self.doors.values())
            A_net -= A_drs
        return A_net

    @property
    def area_gross(self) -> Quantity:
        return self.width * self.height

    @property
    def UA(self) -> float:
        U = self.construction_assembly.U.to('W / (m ** 2 * K)').m
        A = self.area_net.to('m ** 2').m
        return U * A

    def T_adj(self, t: float) -> float:
        return self.T_adj_fun(t)

    def get_conductive_heat_gain(self, t: float, T_int: float) -> dict[str, float]:
        U = self.construction_assembly.U.to('W / (m ** 2 * K)').m
        A = self.area_net.to('m ** 2').m
        T_adj = self.T_adj_fun(t)
        Q = U * A * (T_adj - T_int)
        Q_rad = self.F_rad * Q
        Q_conv = Q - Q_rad
        return {'rad': Q_rad, 'conv': Q_conv}

    def add_door(
        self,
        ID: str,
        width: Quantity,
        height: Quantity,
        construction_assembly: ConstructionAssembly,
        F_rad: float = 0.46
    ) -> None:
        door = InteriorBuildingElement.create(
            ID=ID,
            width=width,
            height=height,
            construction_assembly=construction_assembly,
            F_rad=F_rad,
            T_adj_fun=self.T_adj_fun
        )
        self.doors[door.ID] = door
