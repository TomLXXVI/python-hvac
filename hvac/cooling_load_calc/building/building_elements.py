from __future__ import annotations

import warnings
from typing import Type
from hvac import Quantity
from .weather_data import WeatherData, ExteriorSurface
from .construction_assembly import ConstructionAssembly
from .fenestration import (
    WindowThermalProperties, 
    Window, 
    ExteriorShadingDevice, 
    InteriorShadingDevice
)
from ..thermal_models import (
    ExteriorBuildingElementModel, 
    Units
)

Q_ = Quantity


class ExteriorBuildingElement:
    """Represents an opaque exterior building element.

    An exterior building element is on one side exposed to the outdoor
    environment. The temperature on the exterior surface of the building element
    is the sol-air temperature. The sol-air-temperature will be calculated using
    the given climatic design information valid for the specified design day on 
    the specified geographic location. 

    The location, the design day and all climatic design information are
    encapsulated in an instance of class `WeatherData`.

    The material construction of an exterior building element is defined by its
    construction assembly, an object of class `ConstructionAssembly`.

    The conductive heat transfer through an exterior building element depends on
    its thermal mass (thermal inertia). To determine this time-dependent heat 
    transfer, the exterior building element is modeled by a linear thermal 
    network consisting of temperature nodes, which have a thermal capacitor 
    (a "heat reservoir") and which are interconnected by thermal resistors (see 
    module `thermal_models.ext_build_elem.py`).
    """
    def __init__(self):
        self.name: str = ''
        self.gross_area: Quantity | None = None
        self.net_area: Quantity | None = None
        self.constr_assem: ConstructionAssembly | None = None
        self.ext_surf: ExteriorSurface | None = None
        self.F_rad: float = 0.0
        self.windows: dict[str, Window] = {}
        self.doors: dict[str, ExteriorBuildingElement] = {}
        self._model_factory: Type[ExteriorBuildingElementModel] = ExteriorBuildingElementModel
        self.ext_input_label: str = 'T_sa'  # refers to sol-air temperature
        self.int_input_label: str = 'T_is'  # refers to interior surface temperature
        self.parent: ExteriorBuildingElement | None = None  # used with exterior doors
        
    @classmethod
    def create(
        cls,
        name: str,
        gross_area: Quantity,
        constr_assem: ConstructionAssembly,
        weather_data: WeatherData,
        azimuth_angle: Quantity,
        slope_angle: Quantity,
        surface_color: str = 'dark',
        R_surf: Quantity | None = None,
        a_surf: Quantity | None = None,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = '',
        F_rad: float = 0.46
    ) -> ExteriorBuildingElement:
        """Creates an `ExteriorBuildingElement` object.

        Parameters
        ----------
        name:
            Identifies the exterior building element.
        gross_area:
            Gross surface area of the exterior building element, i.e. also
            including any large openings, such as doors and windows.
        constr_assem:
            Construction assembly the exterior building element is made of.
        weather_data:
            Instance of class `WeatherData` encapsulating the climatic design
            information needed to determine the solar radiation incident on the
            exterior surface of the building element and its sol-air temperature
            during the design day, which was specified on instantiation of the
            `WeatherData` object.
        azimuth_angle:
            The azimuth angle of the exterior building element. South = 0°,
            West = +90° and East = -90°.
        slope_angle:
            Slope angle of the exterior building element. E.g., a vertical
            exterior wall has a slope angle of 90°.
        surface_color: optional
            Either a 'dark' (default) or 'light' colored surface.
        R_surf: optional
            Thermal resistance of the exterior surface film.
        a_surf: optional
            Absorption factor of the exterior surface.
        rho_g: Quantity, default 0.2 frac
            Ground reflectance.
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            Sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.
        F_rad: default 0.46
            Radiative fraction of the conduction heat gain through the exterior
            building element, taken up by the interior thermal mass of the
            space (see ASHRAE Fundamentals 2017, chapter 18, table 14).

        Notes
        -----
        If `R_surf` and `a_surf` are specified, parameter `surface_color` is
        ignored.

        Returns
        -------
        `ExteriorBuildingElement` object.
        """
        ext_build_elem = cls()
        ext_build_elem.name = name
        ext_build_elem.gross_area = gross_area
        ext_build_elem.net_area = gross_area
        ext_build_elem.constr_assem = constr_assem
        ext_build_elem.ext_surf = ExteriorSurface(
            weather_data=weather_data,
            gamma=azimuth_angle,
            beta=slope_angle,
            surface_color=surface_color,
            R_surf=R_surf,
            a_surf=a_surf,
            rho_g=rho_g,
            sky_model=sky_model
        )
        ext_build_elem.F_rad = F_rad
        return ext_build_elem

    @property
    def UA(self) -> Quantity:
        """Returns the steady-state total transmittance of the exterior building
        element.
        """
        U = self.constr_assem.U.to('W / (m**2 * K)')
        A = self.net_area.to('m**2')
        return U * A
    
    def get_input_name(self) -> str:
        """Returns the input label that will be used when creating the zone 
        system. For internal use only.
        """
        return f"{self.ext_input_label}@{self.name}"
    
    def create_model(
        self, 
        ext_input_label: str | None = None, 
        int_input_label: str | None = None, 
        closed: bool = False
    ) -> ExteriorBuildingElementModel:
        """Creates and then returns the linear thermal network model of the 
        exterior building element.
        
        Parameters
        ----------
        ext_input_label: optional
            Name for the temperature node on the exterior side of the building
            element, where the sol-air temperature is present. Always an 
            external temperature node (i.e. an input of the linear thermal 
            network system). If `None`, the default name 'T_sa' is given. 
        int_input_label: optional
            Name for the temperature node on the interior surface side of the 
            building element, where the interior surface temperature is present. 
            Can be either an internal temperature node (part of the linear 
            thermal network) or an external temperature node (i.e. an input of 
            the linear thermal network system). If `None`, the default name 
            `T_is` is given.
        closed:
            Indicates if the linear thermal network model should be closed or
            not. 
            If `closed` is true, the linear thermal network is terminated on its
            interior side by an "interior surface temperature node". In that 
            case the linear thermal network system is a single-input/
            single-output system (SISO) with the input being the sol-air 
            temperature on the exterior side of the building element, and the 
            output of the system is the temperature of the interior surface.
            The interior surface node has no resistor that connects it with an
            external environment (the indoor environment).
            Note that in this case it is assumed that no external heat flow can 
            leave or enter the interior surface node (it is isolated from the 
            external environment); only internal heat transfer between the 
            interior surface node and its preceding node in the linear thermal 
            network is possible.
            If `closed` is false (default), the thermal resistor of the last 
            (building mass) temperature node in the network is connected to an 
            external temperature node on the interior side of the building 
            element. In that case the linear thermal network system has now two 
            inputs (the sol-air temperature on the exterior side and the indoor 
            air temperature on the interior side), and the temperature of the 
            last (building mass) node is now the output of the system.
        
        Notes
        -----
        To determine the conduction heat gain at the interior side of the
        exterior building element for a given indoor air temperature, the
        linear thermal network should be open-ended (`closed=False`).
        
        To connect the linear thermal network model of the exterior building
        element with a zone model, the linear thermal network should be
        closed (`closed=True`). 
        The connection with a zone model is actually done with two separate 
        connections. 
        One connection connects the "interior surface node" of the exterior 
        building element with the "zone-air mass node" of the zone model. This 
        connection has a "convective heat-transfer" resistor (`R_conv`), that 
        represents convective heat transfer between the interior surface and the
        zone air. 
        The other connection connects the "interior surface node" of the 
        exterior building element with the "interior mass node" of the zone 
        model, which models thermal storage of radiative heat inside the 
        interior mass of a zone. This connection has a "radiative heat-transfer"
        resistor (`R_rad`), that represents radiative heat transfer between the 
        interior surface and the interior mass of the zone. 
        The convective resistor (`R_conv`) is the thermal resistance of the 
        interior surface film. The value of the radiative resistor is determined
        by the specified parameter `F_rad`. If the given construction assembly
        has been configured with an interior surface film layer, the attributes 
        `R_conv` and `R_rad` of the `ExteriorBuildingElementModel` object are 
        instances of class `Resistor`. Otherwise, these attributes will remain 
        `None`. 
        """
        model = self._model_factory(
            self.name,
            self.constr_assem,
            self.net_area,
            self.F_rad,
            ext_input_label or self.ext_input_label,
            int_input_label or self.int_input_label,
            closed
        )
        return model
   
    def conduction_heat_gain(self, T_za: Quantity) -> tuple[Quantity, ...]:
        """Returns the hourly values of conductive heat gain at the interior 
        side of the exterior building element on the selected design day, 
        together with its convective and radiative components.
        
        Parameters
        ----------
        T_za:
            Zone-air temperature, i.e. the indoor air temperature at the 
            interior side of the exterior building element.
        
        Returns
        -------
        Q_dot_cond:
            `Quantity` array containing the hourly average values of conduction 
            heat gain on the selected design day and with the specified zone-air
            temperature.
        Q_dot_cond_conv:
            Part of the conduction heat gain which is transferred by convection 
            to the zone air.
        Q_dot_cond_rad:
            Part of the conduction heat gain which is transferred by radiation 
            to the interior thermal mass of the zone.
        """
        T_sa = Quantity.from_list([
            self.ext_surf.T_sa(hr * 3600) for hr in range(0, 24)
        ])
        T_za = Quantity.from_list([T_za] * len(T_sa))
        model = self._model_factory(
            self.name, 
            self.constr_assem, 
            self.net_area, 
            self.F_rad,
            input_label='sa',
            output_label='za',
            closed=False
        )
        time, T_n = model.node_temperatures(T_sa, T_za, num_cycles=5)
        T_nf = T_n[-1, -1]  # temperature of last node at final hour
        i, i_max = 0, 10
        # repeat calculations until tolerance condition is met or maximum number
        # of iterations is reached
        while i < i_max:
            time, T_n = model.node_temperatures(
                T_sa, T_za, 1, 
                T_n0=T_n[:, -1].flatten()
            )
            T_nf_new = T_n[-1, -1]
            dev = T_nf_new - T_nf
            if abs(dev) < Q_(0.1, 'K'):
                break
            T_nf = T_nf_new
        else:
            warnings.warn(
                'No accurate solution could be determined.',
                category=UserWarning
            )
        _, R = model.nodes[-1].resistors_out[-1]  # thermal resistor between last node and zone air
        T_nf = T_n[-1, :]  # hourly temperatures of the last node
        Q_dot_cnd = (T_nf - T_za.to(Units.unit_T)) / R.value.to(Units.unit_R)
        Q_dot_cnd_rad = self.F_rad * Q_dot_cnd
        Q_dot_cnd_conv = Q_dot_cnd - Q_dot_cnd_rad
        return Q_dot_cnd, Q_dot_cnd_conv, Q_dot_cnd_rad
    
    def add_window(
        self,
        name: str,
        width: Quantity,
        height: Quantity,
        props: WindowThermalProperties,
        F_rad: float = 0.46,
        ext_shading: ExteriorShadingDevice | None = None,
        int_shading: InteriorShadingDevice | None = None
    ) -> Window:
        """Adds a window to the exterior building element.

        Parameters
        ----------
        name:
            Identifies the window. `Window` objects are hold in a dictionary 
            `windows`. The name of the window is used as the key in this 
            dictionary.
        width:
            The width of the window.
        height:
            The height of the window.
        props:
            A `WindowsThermalProperties` object that holds the thermal
            and solar properties of the window.
        F_rad: default 0.46
            Radiative fraction of solar heat gain through window to the interior
            thermal mass of the space. (see ASHRAE Fundamentals 2017, chapter 18,
            table 14).
        ext_shading: default None
            E.g. an overhang or recessed window. See: class `ExternalShadingDevice`.
            If a window has an external shading device, or if the window is
            recessed, part of the window may be shaded depending on the position
            of the sun during the day.
        int_shading: default None
            E.g. a louvered shade, roller shade, drapery, or insect screen.
            See: class `InteriorShadingDevice`.

        Returns
        -------
        The `Window` object that was added to the `ExteriorBuildingElement`
        object.
        """
        # Create `Window` object.
        window = Window.create(
            name=name,
            azimuth_angle=self.ext_surf.gamma,
            slope_angle=self.ext_surf.beta,
            width=width,
            height=height,
            weather_data=self.ext_surf.weather_data,
            props=props,
            F_rad=F_rad,
            ext_shading=ext_shading,
            int_shading=int_shading
        )
        # Add a reference to the exterior building element the window is part of.
        window.parent = self
        # Hold the windows that are part of the exterior building element inside
        # a dictionary.
        self.windows[name] = window
        # When a window is added to the exterior building element, the opaque
        # surface area of the exterior building element is reduced. Therefore, 
        # it will be needed to recreate the linear thermal network model of the 
        # exterior building element each time a window has been added.
        self.net_area -= window.area
        return window
    
    def add_door(
        self,
        name: str,
        width: Quantity,
        height: Quantity,
        constr_assem: ConstructionAssembly,
        surface_color: str = 'dark',
        R_surf: Quantity | None = None,
        a_surf: Quantity | None = None,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = '',
        F_rad: float = 0.46
    ) -> ExteriorBuildingElement:
        """Adds a door to the exterior building element.
        
        An exterior door is also represented by the `ExteriorBuildingElement` 
        class.
        
        Parameters
        ----------
        name:
            Name to identify the door.
        width:
            Width of the door.
        height:
            Height of the door. 
        constr_assem:
            Construction assembly the door is made of.
        surface_color: optional
            Either a 'dark' (default) or 'light' colored surface.
        R_surf: optional
            Thermal resistance of the exterior surface film.
        a_surf: optional
            Absorption factor of the exterior surface.
        rho_g: Quantity, default 0.2 frac
            Ground reflectance.
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            Sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.
        F_rad: default 0.46
            Radiative fraction of the conduction heat gain through the door
            building element, taken up by the interior thermal mass of the
            space (see ASHRAE Fundamentals 2017, chapter 18, table 14).

        Notes
        -----
        If `R_surf` and `a_surf` are specified, parameter `surface_color` is
        ignored.

        Returns
        -------
        `ExteriorBuildingElement` object.
        """
        # Create the `ExteriorBuildingElement` object to represent the exterior 
        # door.
        door = ExteriorBuildingElement.create(
            name=name,
            constr_assem=constr_assem,
            gross_area=width * height,
            weather_data=self.ext_surf.weather_data,
            azimuth_angle=self.ext_surf.gamma,
            slope_angle=self.ext_surf.beta,
            surface_color=surface_color,
            R_surf=R_surf,
            a_surf=a_surf,
            rho_g=rho_g,
            sky_model=sky_model,
            F_rad=F_rad
        )
        # Add a reference to the exterior building element the door is part of.
        door.parent = self
        # Hold the exterior doors that are part of the exterior building element
        # inside a dictionary.
        self.doors[door.name] = door
        # As a door is added to the exterior building element, the opaque
        # surface area of the exterior building element, which is exposed to
        # conduction heat transfer, is reduced. Therefore, it is needed to
        # recreate the linear thermal network model of the exterior building
        # element with the net surface area of the exterior building element.
        self.net_area -= door.gross_area
        return door


class InteriorBuildingElement:
    """Represents a building element inside the building (e.g., interior walls,
    ceilings).
    """
    def __init__(self):
        self.name: str = ''
        self.adj_zone_name: str = ''
        self.constr_assem: ConstructionAssembly | None = None
        self.gross_area: Quantity | None = None
        self.net_area: Quantity | None = None
        self.F_rad: float = 0.46
        self.doors: dict[str, InteriorBuildingElement] = {}
        self.T_za_adj: Quantity | None = None
        self.V_dot_trf: Quantity | None = None
        self.parent: InteriorBuildingElement | None = None
        self.input_label: str = 'T_za'  # refers to the zone-air temperature in the adjacent zone.
        
    @classmethod
    def create(
        cls,
        name: str,
        adjacent_zone_name: str,
        constr_assem: ConstructionAssembly,
        gross_area: Quantity,
        adjacent_zone_temperature: Quantity,
        V_dot_trf: Quantity | None = None,
        F_rad: float = 0.46
    ) -> InteriorBuildingElement:
        """
        Creates an `InteriorBuildingElement` object.

        Parameters
        ----------
        name:
            Name to identify the interior building element.
        adjacent_zone_name:
            Name given to the adjacent zone.
        constr_assem:
            The construction assembly the interior building element is made of.
        gross_area:
            The gross surface area of the interior building element, i.e.
            including any large openings such as doors and windows.
        adjacent_zone_temperature:
            The zone-air temperature on the other side of the interior building
            element (i.e. the zone-air temperature of the adjacent zone).
        V_dot_trf: optional
            Volume flow rate of air transferred from the adjacent zone. If 
            `None`, no air is transferred. 
        F_rad: default 0.46
            Radiative fraction of the conduction heat gain through the interior
            building element, taken up by the interior thermal mass of the
            space. (see ASHRAE Fundamentals 2017, chapter 18, table 14).

        Returns
        -------
        `InteriorBuildingElement` object.
        """
        int_build_elem = cls()
        int_build_elem.name = name
        int_build_elem.adj_zone_name = adjacent_zone_name
        int_build_elem.constr_assem = constr_assem
        int_build_elem.gross_area = gross_area
        int_build_elem.net_area = gross_area
        int_build_elem.T_za_adj = adjacent_zone_temperature
        int_build_elem.V_dot_trf = V_dot_trf
        int_build_elem.F_rad = F_rad
        return int_build_elem
    
    def get_input_name(self) -> str:
        """Returns the input label that will be used when creating the zone 
        system. For internal use only.
        """
        return f"{self.input_label}@{self.adj_zone_name}"
    
    def add_door(
        self,
        name: str,
        width: Quantity,
        height: Quantity,
        constr_assem: ConstructionAssembly,
        F_rad: float = 0.46
    ) -> InteriorBuildingElement:
        """Adds a door into the interior building element.

        An interior door is also represented by an `InteriorBuildingElement` 
        object.

        Parameters
        ----------
        name:
            Identifies the door.
            Doors are hold in a dictionary `self.doors`. The ID of the door is
            used as the key in the dictionary.
        width:
            The width of the door.
        height:
            The height of the door.
        constr_assem:
            The construction assembly the door is made off.
        F_rad: default 0.46
            Radiative fraction of the conduction heat gain through the door,
            taken up by the interior thermal mass of the space. (see ASHRAE
            Fundamentals 2017, chapter 18, table 14).
        """
        door = InteriorBuildingElement.create(
            name=name,
            adjacent_zone_name=self.adj_zone_name,
            constr_assem=constr_assem,
            gross_area=width * height,
            adjacent_zone_temperature=self.T_za_adj,
            F_rad=F_rad
        )
        # Add a reference to the interior building element the door is part of.
        door.parent = self
        # Hold the interior doors that are part of the interior building element
        # inside a dictionary.
        self.doors[door.name] = door
        # As a door is added to the interior building element, the surface area 
        # of the interior building element is reduced.
        self.net_area -= door.gross_area
        return door

    @property
    def UA(self) -> Quantity:
        """Returns the total transmittance of the interior building element."""
        U = self.constr_assem.U.to('W / (m**2 * K)')
        A = self.net_area.to('m**2')
        return U * A
