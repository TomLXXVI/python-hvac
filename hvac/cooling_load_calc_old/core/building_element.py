from __future__ import annotations

from typing import Callable

import pandas as pd

from hvac import Quantity
from .construction_assembly import ConstructionAssembly
from .fenestration import (
    WindowThermalProperties,
    Window,
    ExteriorShadingDevice,
    InteriorShadingDevice
)
from .thermal_models import (
    ExteriorBuildingElementLTN,
    ExteriorBuildingElementLTNBuilder
)
from .weather_data import WeatherData, ExteriorSurface

Q_ = Quantity


class ExteriorBuildingElement:
    """\
    Models an opaque exterior building element.

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
    (a "heat reservoir") and are interconnected by thermal resistors (see module 
    `thermal_models.exterior_building_element.py).
    """
    def __init__(self):
        self.ID: str = ''
        self.T_zone: Callable[[float], Quantity] | None = None
        self.ext_surf: ExteriorSurface | None = None
        self.constr_assem: ConstructionAssembly | None = None
        self.model: ExteriorBuildingElementLTN | None = None
        self.gross_area: Quantity | None = None
        self.net_area: Quantity | None = None
        self.F_rad: float = 0.46
        self.windows: dict[str, Window] = {}
        self.doors: dict[str, ExteriorBuildingElement] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        T_zone: Quantity | Callable[[float], Quantity],
        constr_assem: ConstructionAssembly,
        gross_area: Quantity,
        weather_data: WeatherData,
        gamma: Quantity,
        beta: Quantity,
        surface_color: str = 'dark',
        R_surf: Quantity | None = None,
        a_surf: Quantity | None = None,
        rho_g: Quantity = Q_(0.2, 'frac'),
        sky_model: str = '',
        F_rad: float = 0.46
    ) -> ExteriorBuildingElement:
        """\
        Creates an `ExteriorBuildingElement` object.

        Parameters
        ----------
        ID:
            Identifies the exterior building element.
        T_zone:
            The zone air temperature. A function can be passed with signature
            `f(t_sol_sec: float) -> Quantity` which takes the solar time in
            seconds and returns the temperature in the zone as a `Quantity`
            object. This may allow a time-variable temperature in the zone.
        constr_assem:
            The construction assembly the exterior building element is made off.
        gross_area:
            The gross surface area of the exterior building element, i.e.
            including any large openings such as doors and windows.
        weather_data:
            `WeatherData` instance that encapsulates the climatic design
            information needed to determine the solar radiation incident on the
            exterior surface of the building element and its sol-air temperature
            during the design day, which was specified on instantiation of the
            `WeatherData` object.
        gamma:
            The azimuth angle of the exterior building element. South = 0°,
            West = +90° and East = -90°.
        beta:
            The slope angle of the exterior building element. E.g., a vertical
            exterior wall has a slope angle of 90°.
        surface_color: optional
            Either a 'dark' (default) or 'light' colored surface.
        R_surf: optional
            Thermal resistance of the exterior surface film.
        a_surf: optional
            Absorption factor of the exterior surface.
        rho_g: Quantity, default 0.2 frac
            ground reflectance
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.
        F_rad: default 0.46
            Radiative fraction of the conduction heat gain through the exterior
            building element, taken up by the interior thermal mass of the
            space. (see ASHRAE Fundamentals 2017, chapter 18, table 14).

        Notes
        -----
        If `R_surf` and `a_surf` are specified, parameter `surface_color` is
        ignored.

        Returns
        -------
        `ExteriorBuildingElement` object.
        """
        ebe = cls()
        ebe.ID = ID
        if isinstance(T_zone, Quantity):
            ebe.T_zone = lambda t: T_zone
        else:
            ebe.T_zone = T_zone
        ebe.gross_area = gross_area
        ebe.net_area = gross_area
        ebe.constr_assem = constr_assem
        ebe.F_rad = F_rad
        ebe.ext_surf = ExteriorSurface(
            weather_data=weather_data,
            gamma=gamma,
            beta=beta,
            surface_color=surface_color,
            R_surf=R_surf,
            a_surf=a_surf,
            rho_g=rho_g,
            sky_model=sky_model
        )
        # Create the linear thermal network model of the exterior building element:
        ebe.model = ExteriorBuildingElementLTNBuilder.build(constr_assem, gross_area)
        return ebe

    def solve(
        self,
        dt_hr: float = 1.0,
        num_cycles: int | None = 5,
        units: dict[str, str] | None = None
    ) -> tuple[pd.DataFrame, pd.DataFrame] | None:
        """\
        Solves for the node temperatures and the conductive heat transfer
        through the exterior building element on the specified design day using
        the linear thermal network model of the exterior building element.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step for the calculations is one quarter of an hour.
            The default value is 1 hour.
        num_cycles:
            Number of diurnal calculation cycles before results are returned.
            The default number is 5 cycles.
        units: optional
            Dictionary with the units in which node temperatures and heat flows
            should be expressed. If None, default units are:
            {'T': 'degC', 'Q_dot': 'W'}

        Returns
        -------
        Either 2-tuple or None.

        The 2-tuple has following items:
        - A table (Pandas `DataFrame` object) with the node temperatures at each
        time step of the last calculation cycle.
        - A table (Pandas `DataFrame` object) with the heat flows between nodes
        at each time step of the last calculation cycle.

        Returns None if the exterior building element has no linear thermal
        network.

        Notes
        -----
        To determine the heat transfer through the exterior building element,
        the linear thermal network model of the exterior building element is
        solved.
        At the exterior side of the building element, the known input function
        is the sol-air temperature, being a function of time.
        At the interior side of the building element, the known input function
        is the controlled zone air temperature. If a single `Quantity` was 
        passed to `T_zone` when instantiating the `ExteriorBuildingElement` 
        object, it is assumed that the temperature control of the zone keeps 
        the zone air temperature constant in time.
        
        The calculations using the linear thermal network model, being a linear
        system of first order differential equations, start at time t = 0 s, 
        which is at midnight of the design day. The second time moment is at 
        `dt_hr` * 3600 s, the third time moment is at 2 * `dt_hr` * 3600 s, 
        and so on. The total number of steps is determined as 24 hr / `dt_hr`,
        rounded to the nearest integer. So, the duration of the design day is
        divided in (24 hr / `dt_hr`) equal time intervals. When the last
        calculation step at t = (n - 1) * `dt_hr` * 3600 s is finished,
        one diurnal calculation cycle is finished.
        
        To solve the linear thermal network for the node temperatures, it is
        necessary to know the initial values of all the node temperatures at 
        t = 0 s. The initial values of the node temperatures are internally set
        all equal to the sol-air temperature at the exterior surface of the
        building element at t = 0 s. As the values of the node temperatures at
        t > 0 s will depend on these first initial values, the diurnal cycle 
        needs to be repeated multiple times to reduce the influence of these
        arbitrarily chosen initial values. The number of repetitions can be set 
        with parameter `num_cycles`. 
        The results from the final time step in the previous cycle are then used
        as initial values for the next cycle. It is as if the same daily weather
        cycle of the design day repeats itself for several days in a row.
        """
        if self.model is not None and num_cycles is not None:
            self.model.solve(
                num_steps=int(round(24 / dt_hr)),
                dt_hr=dt_hr,
                T_ext=self.ext_surf.T_sa,
                T_zone=self.T_zone,
                init_values=[
                    [self.ext_surf.T_sa(0)]
                    * len(self.model) for _ in range(2)
                ],
                num_cycles=num_cycles
            )
            if units is None:
                units = {'T': 'degC', 'Q_dot': 'W'}
            T_node_table = self.model.get_node_temperatures(units['T'])
            Q_dot_table = self.model.get_heat_flows(units['Q_dot'])
            return T_node_table, Q_dot_table
        return None

    def conductive_heat_gain(
        self,
        dt_hr: float = 1.0,
        num_cycles: int | None = 5
    ) -> tuple[Quantity, Quantity, Quantity]:
        """\
        Returns the conductive heat rate output at the interior side of the
        exterior building element, and also its convective and radiative
        component, for each time index k of the design day.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step for the calculations is one quarter of an hour.
            The default value is 1 hour.
        num_cycles:
            Number of diurnal calculation cycles before the results of the last
            diurnal cycle are returned. If `None`, the static heat conduction
            through the building element is calculated. The default number is
            5 cycles.

        Returns
        -------
        3-tuple with:
        - a `Quantity`-array with the conductive heat gains at each time index
        - a `Quantity`-array with the convective components
        - a `Quantity`-array with the radiative components

        Notes
        -----
        The conductive heat gain to the zone is the heat that flows from the
        interior surface node to the zone air. It is the last column, titled
        'ISN', in the heat flow table returned by method `solve()`.

        The time in decimal hours from midnight that corresponds with a given
        time index depends on the time step `dt_hr` (`t_hr_dec = k * `dt_hr`).
        """
        res = self.solve(dt_hr, num_cycles)
        if res is not None:
            _, Q_dot_table = res
            Q_dot_cond = Q_(Q_dot_table['ISN'].values, 'W')
            Q_dot_rad = self.F_rad * Q_dot_cond
            Q_dot_conv = (1.0 - self.F_rad) * Q_dot_cond
        else:
            # The exterior building element has no thermal linear network.
            # Fall back on a steady-state -only thermal resistance- calculation:
            num_steps = int(round(24 / dt_hr))
            dt_sec = dt_hr * 3600
            U = self.constr_assem.U.to('W / (m**2 * K)')
            A = self.net_area.to('m**2')
            dT = Quantity.from_list([
                self.ext_surf.T_sa(k * dt_sec).to('K') - self.T_zone(k * dt_sec).to('K')
                for k in range(num_steps)
            ])
            Q_dot_cond = U * A * dT
            Q_dot_rad = self.F_rad * Q_dot_cond
            Q_dot_conv = (1.0 - self.F_rad) * Q_dot_cond
        return Q_dot_cond, Q_dot_conv, Q_dot_rad

    def add_window(
        self,
        ID: str,
        width: Quantity,
        height: Quantity,
        props: WindowThermalProperties,
        F_rad: float = 0.46,
        ext_shading: ExteriorShadingDevice | None = None,
        int_shading: InteriorShadingDevice | None = None
    ) -> Window:
        """
        Adds a window to the exterior building element.

        Parameters
        ----------
        ID:
            Identifies the window.
            `Window` objects are hold in a dictionary `self.windows`. The ID
            of the window is used as the key in the dictionary.
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
        window = Window.create(
            ID=ID,
            T_zone=self.T_zone,
            gamma=self.ext_surf.gamma,
            beta=self.ext_surf.beta,
            width=width,
            height=height,
            weather_data=self.ext_surf.weather_data,
            props=props,
            F_rad=F_rad,
            ext_shading=ext_shading,
            int_shading=int_shading
        )
        self.windows[ID] = window
        # As a window is added to the exterior building element, the opaque
        # surface area of the exterior building element, which is exposed to
        # conduction heat transfer, is reduced. Therefore, it is needed to
        # recreate the linear thermal network model of the exterior building
        # element with the net surface area of the exterior building element.
        self.net_area -= window.area
        self.model = ExteriorBuildingElementLTNBuilder.build(self.constr_assem, self.net_area)
        return window

    def add_door(
        self,
        ID: str,
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
        """
        Adds a door to the exterior building element.

        A door is also represented as an `ExteriorBuildingElement` object.

        Parameters
        ----------
        ID:
            Identifies the door.
            Doors are hold in a dictionary `self.doors`. The ID of the door is
            used as the key in the dictionary.
        width:
            The width of the door.
        height:
            The height of the door.
        constr_assem:
            The construction assembly the door is made off.
        surface_color: optional
            Either a 'dark' (default) or 'light' colored surface.
        R_surf: optional
            Thermal resistance of the exterior surface film.
        a_surf: optional
            Absorption factor of the exterior surface.
        rho_g: Quantity, default 0.2 frac
            ground reflectance
        sky_model: str, {'anisotropic.hdkr, 'isotropic'}, optional
            sky model to be used for calculating the irradiation on the tilted
            surface. By default, the anisotropic sky model according to Perez
            is used.
        F_rad: default 0.46
            Radiative fraction of the conduction heat gain through the exterior
            building element, taken up by the interior thermal mass of the
            space. (see ASHRAE Fundamentals 2017, chapter 18, table 14).

        Notes
        -----
        If `R_surf` and `a_surf` are specified, parameter `surface_color` is
        ignored.
        """
        door = ExteriorBuildingElement.create(
            ID=ID,
            T_zone=self.T_zone,
            constr_assem=constr_assem,
            gross_area=width * height,
            weather_data=self.ext_surf.weather_data,
            gamma=self.ext_surf.gamma,
            beta=self.ext_surf.beta,
            surface_color=surface_color,
            R_surf=R_surf,
            a_surf=a_surf,
            rho_g=rho_g,
            sky_model=sky_model,
            F_rad=F_rad
        )
        self.doors[door.ID] = door
        # As a door is added to the exterior building element, the opaque
        # surface area of the exterior building element, which is exposed to
        # conduction heat transfer, is reduced. Therefore, it is needed to
        # recreate the linear thermal network model of the exterior building
        # element with the net surface area of the exterior building element.
        self.net_area -= door.gross_area
        self.model = ExteriorBuildingElementLTNBuilder.build(self.constr_assem, self.net_area)
        return door

    @property
    def UA(self) -> Quantity:
        """Returns the total transmittance of the exterior building element."""
        U = self.constr_assem.U.to('W / (m**2 * K)')
        A = self.net_area.to('m**2')
        return U * A


class InteriorBuildingElement:
    """
    Represents a building element inside the building (e.g., interior walls,
    ceilings).
    """
    def __init__(self):
        self.ID: str = ''
        self.T_zone: Callable[[float], Quantity] | None = None
        self.T_adj: Callable[[float], Quantity] | None = None
        self.constr_assem: ConstructionAssembly | None = None
        self.gross_area: Quantity | None = None
        self.net_area: Quantity | None = None
        self.F_rad: float = 0.46
        self.doors: dict[str, InteriorBuildingElement] = {}

    @classmethod
    def create(
        cls,
        ID: str,
        T_zone: Callable[[float], Quantity],
        T_adj: Callable[[float], Quantity],
        constr_assem: ConstructionAssembly,
        gross_area: Quantity,
        F_rad: float = 0.46
    ) -> InteriorBuildingElement:
        """
        Creates an `InteriorBuildingElement` object.

        Parameters
        ----------
        ID:
            Identifies the interior building element.
        T_zone:
            The zone air temperature, being a function with signature
            `f(t_sol_sec: float) -> Quantity` which takes the solar time in
            seconds and returns the temperature in the zone as a `Quantity`
            object. This may allow a time-variable temperature in the zone.
        T_adj:
            The temperature in the adjacent space, being a function with
            signature `f(t_sol_sec: float) -> Quantity` which takes the solar
            time in seconds and returns the temperature in the adjacent space
            as a `Quantity` object. This allows a time-variable temperature in
            the adjacent space that may be an unconditioned space.
        constr_assem:
            The construction assembly the interior building element is made of.
        gross_area:
            The gross surface area of the interior building element, i.e.
            including any large openings such as doors and windows.
        F_rad: default 0.46
            Radiative fraction of the conduction heat gain through the interior
            building element, taken up by the interior thermal mass of the
            space. (see ASHRAE Fundamentals 2017, chapter 18, table 14).

        Returns
        -------
        `InteriorBuildingElement` object.
        """
        ibe = cls()
        ibe.ID = ID
        ibe.T_zone = T_zone
        ibe.T_adj = T_adj
        ibe.constr_assem = constr_assem
        ibe.gross_area = gross_area
        ibe.net_area = gross_area
        ibe.F_rad = F_rad
        return ibe

    def conductive_heat_gain(
        self,
        dt_hr: float = 1.0
    ) -> tuple[Quantity, Quantity, Quantity]:
        """
        Returns the conductive heat output at the interior side of the
        interior building element, and also its convective and radiative
        component, at each time index k of the design day.

        Parameters
        ----------
        dt_hr:
            Time step expressed as a fraction of 1 hour, e.g., `dt_hr` = 1/4
            means the time step for the calculations is one quarter of an hour.
            The default value is 1 hour.

        Returns
        -------
        3-tuple with:
        - a `Quantity`-array with the conductive heat gains at each time index
        - a `Quantity`-array with the convective components
        - a `Quantity`-array with the radiative components

        Notes
        -----
        The time in decimal hours from midnight of the design day that
        corresponds with a given time index depends on the time step `dt_hr`
        (`t_hr_dec = k * `dt_hr`).
        """
        num_steps = int(round(24 / dt_hr))
        dt_sec = dt_hr * 3600
        U = self.constr_assem.U.to('W / (m**2 * K)')
        A = self.net_area.to('m**2')
        dT = Quantity.from_list([
            self.T_adj(k * dt_sec).to('K') - self.T_zone(k * dt_sec).to('K')
            for k in range(num_steps)
        ])
        Q_dot_cond = U * A * dT
        Q_dot_rad = self.F_rad * Q_dot_cond
        Q_dot_conv = (1.0 - self.F_rad) * Q_dot_cond
        return Q_dot_cond, Q_dot_conv, Q_dot_rad

    def add_door(
        self,
        ID: str,
        width: Quantity,
        height: Quantity,
        constr_assem: ConstructionAssembly,
        F_rad: float = 0.46
    ) -> InteriorBuildingElement:
        """
        Adds a door to the interior building element.

        A door is also represented as an `InteriorBuildingElement` object.

        Parameters
        ----------
        ID:
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
            ID=ID,
            T_zone=self.T_zone,
            T_adj=self.T_adj,
            constr_assem=constr_assem,
            gross_area=width * height,
            F_rad=F_rad
        )
        self.doors[door.ID] = door
        # As a door is added to the exterior building element, the opaque
        # surface area of the interior building element, which is exposed to
        # conduction heat transfer, is reduced.
        self.net_area -= door.gross_area
        return door

    @property
    def UA(self) -> Quantity:
        """Returns the total transmittance of the interior building element."""
        U = self.constr_assem.U.to('W / (m**2 * K)')
        A = self.net_area.to('m**2')
        return U * A
