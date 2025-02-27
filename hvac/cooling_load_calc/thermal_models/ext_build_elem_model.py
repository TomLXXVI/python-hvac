import numpy as np
import control as ct
from hvac import Quantity
from ..building.construction_assembly import (
    ConstructionAssembly,
    ConstructionLayer,
    SurfaceFilm
)
from .thermal_network import (
    LinearThermalNetwork,
    TemperatureNode,
    Resistor,
    Capacitor,
    ConnectionSide
)
from .units import Units

Q_ = Quantity


class ExteriorBuildingElementModel(LinearThermalNetwork):
    """Models a flat, opaque exterior building element as a linear thermal 
    network that can also be connected to the model of a temperature zone in 
    the building (see class `ZoneModel`).
    """
    def __init__(
        self,
        name: str,
        constr_assem: ConstructionAssembly,
        area: Quantity,
        F_rad: float,
        input_label: str,
        output_label: str,
        closed: bool = True
    ) -> None:
        """Creates an `ExteriorBuildingElementModel` instance.
        
        Parameters
        ----------
        name:
            Name that identifies the model.
        constr_assem:
            Construction assembly the exterior building element is made of.
        area:
            The (exterior) surface area of the exterior building element.
        F_rad: default 0.46
            Radiative fraction of the conduction heat gain through the exterior
            building element, taken up by the interior thermal mass of the
            space (see ASHRAE Fundamentals 2017, chapter 18, table 14).
        input_label:
            Name for the temperature node on the exterior side of the building
            element, where the sol-air temperature is present. Always an 
            external temperature node (i.e. an input of the linear thermal 
            network).
        output_label:
            Name for the temperature node on the interior side of the building
            element, where the interior surface temperature is present. Can be
            either an internal temperature node (part of the linear thermal 
            network) or an external temperature node (i.e. an input of the 
            linear thermal network).
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
            If `closed` is false, the thermal resistor of the last (building 
            mass) temperature node in the network is connected to an external 
            temperature node on the interior side of the building element. In
            that case the linear thermal network system has now two inputs (the 
            sol-air temperature on the exterior side and the indoor air 
            temperature on the interior side), and the temperature of the last 
            (building mass) node is now the output of the system.
        
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
        self.name = name
        self.constr_assem = constr_assem
        self.area = area
        self.F_rad = F_rad
        self.ext_input_label = input_label
        self.int_input_label = output_label
        self.R_rad: Resistor | None = None
        self.R_conv: Resistor | None = None
        layers = self._get_construction_layers(closed)
        nodes = self._create_nodes(layers, input_label, output_label, closed)
        super().__init__(name, tuple(nodes))
    
    def _get_construction_layers(self, closed: bool) -> list[ConstructionLayer]:
        """Returns the construction layers from which the temperature nodes in
        the linear thermal network model will be derived. If the linear thermal 
        network should be closed, also removes the interior surface film from 
        the linear network model and the linear network is closed with an 
        "interior surface node".
        Furthermore, based on parameter `F_rad`, determines the convective 
        resistor `R_conv` and radiative resistor `R_rad`, which are intended to 
        connect the model of the exterior building element with a zone model 
        (for this the linear thermal network of the exterior building element
        must be closed). 
        """
        if self.constr_assem.layers is not None:
            layers = list(self.constr_assem.layers.values())
            # Remove the interior surface film when making the model of the
            # exterior building element. It will be replaced by a convective
            # and radiative thermal resistor.
            if isinstance(layers[-1], SurfaceFilm):
                if closed:
                    int_surf_film = layers.pop()
                else:
                    int_surf_film = layers[-1]
                R_isf = int_surf_film.R / self.area
                # self.R_conv = Resistor(R_isf)
                # self.R_rad = Resistor((1 - self.F_rad) / self.F_rad * R_isf)
                self.R_conv = Resistor(R_isf / (1 - self.F_rad))
                self.R_rad = Resistor(R_isf / self.F_rad)
            return layers
        else:
            raise ValueError(
                f"Construction assembly '{self.constr_assem.ID}' has no layers "
                f"to build the linear thermal network from."
            )
    
    def _create_nodes(
        self,
        layers: list[ConstructionLayer],
        ext_input_label: str,
        int_input_label: str,
        closed: bool = True
    ) -> list[TemperatureNode]:
        """Creates a list of temperature nodes from the construction layers
        in the construction assembly.
        """         
        # A construction layer with thermal capacity has a number of slices (at 
        # least 1 or more). Each slice represents a temperature node with a 
        # resistor on its left and a resistor on its right side. A layer without
        # thermal capacity has only a resistor. A list of resistors and 
        # temperature nodes results, e.g. `[R, (R, N, R), R, (R, N, R), ... R]`.
        nodes_a = []
        for layer in layers:
            if layer.C.magnitude > 0.0:
                n = layer.num_slices
                R_slice = layer.R / (2 * n) * (1 / self.area)
                C_slice = layer.C / n * self.area
                for i in range(n):
                    resistor_left = Resistor(R_slice)
                    node = TemperatureNode(
                        name=f"{layer.ID}[{i + 1}]",
                        network_name=self.name,
                        capacitor=Capacitor(C_slice)
                    )
                    resistor_right = Resistor(R_slice)
                    nodes_a.extend([resistor_left, node, resistor_right])
            else:
                nodes_a.append(Resistor(layer.R * (1 / self.area)))
        
        # Multiple resistors in series between two nodes are added and replaced
        # by a single resistor. A new list results with alternating one resistor
        # and one node, e.g. `[R, N, R, N, R, N, R, ... R]`
        nodes_b: list[Resistor | TemperatureNode | int] = []
        resistors: list[Resistor] = []
        for item in nodes_a:
            if isinstance(item, Resistor):
                resistors.append(item)
            if isinstance(item, TemperatureNode):
                if resistors:
                    resistor = sum(resistors)
                    nodes_b.extend([resistor, item])
                    resistors = []
                else:
                    nodes_b.append(item)
        if resistors:
            resistor = sum(resistors)
            nodes_b.append(resistor)
        
        # Get indices of the nodes in the list `[R, N, R, N, R, N, R, ... R]`
        node_indices = [
            i 
            for i, item in enumerate(nodes_b) 
            if isinstance(item, TemperatureNode)
        ]
        # A linear thermal network can be open-ended on the interior side, or it
        # can be terminated by an "interior surface node" that must have a small 
        # thermal capacity to be recognized as a temperature node (state 
        # variable) in the state-space representation of the linear thermal 
        # network.
        is_node = TemperatureNode(
            name=int_input_label,
            network_name=self.name,
            capacitor=Capacitor(Q_(1.0, 'J / K')),
        )
        # Configure the temperature nodes by connecting the resistor on the left
        # side and on the right side of each node to this node.
        for i in node_indices:
            resistor_left = nodes_b[i - 1]
            resistor_right = nodes_b[i + 1]
            node = nodes_b[i]
            if i - 2 < 0:
                pre_node = ext_input_label  # preceding node
            else:
                pre_node = nodes_b[i - 2]
            try:
                next_node = nodes_b[i + 2]  # next node
            except IndexError:
                if closed:
                    # if the lineair network is closed, the next and also last 
                    # node of the network will be the "interior surface node" 
                    next_node = is_node
                else:
                    # `int_node` is a string that refers to an external node at
                    # the interior side of the linear thermal network
                    next_node = int_input_label
            node.connect_resistor(
                other_terminal=pre_node,
                resistor=resistor_left,
                this_side=ConnectionSide.IN
            )
            node.connect_resistor(
                other_terminal=next_node,
                resistor=resistor_right,
                this_side=ConnectionSide.OUT
            )
        # Configure the "interior surface node".
        is_node.connect_resistor(
            other_terminal=nodes_b[-2],
            resistor=nodes_b[-1],
            this_side=ConnectionSide.IN
        )
        # Move the configured temperatures nodes to a separate list.
        nodes_c = [nodes_b[i] for i in node_indices]
        if closed:
            # If the linear thermal network is closed, also append the interior 
            # surface node to this list.
            nodes_c.append(is_node)
        return nodes_c
    
    def node_temperatures(
        self, 
        T_sa: Quantity, 
        T_za: Quantity, 
        num_cycles: int = 1,
        T_n0: Quantity | None = None
    ) -> Quantity:
        """Returns the node temperatures at each hour of the selected design 
        day.

        Parameters
        ----------
        T_sa:
            `Quantity` array with the hourly average sol-air temperature at 
            each hour of the design day.
        T_za:
            `Quantity` array with the zone air temperature at each hour of the
            design day.
        num_cycles:
            Number of daily cycles to repeat before the results of the final 
            daily cycle will be returned. 
        T_n0: optional
            `Quantity` array with the initial temperature of each node in the
            network. By default, the initial temperature of each node is zero. 

        Returns
        -------
        time:
            Numpy array with the hours of the day (from 0 h until 23 h).
        T_n:
            2D-`Quantity` array with, for each node in the network, its 
            temperature at each hour of the day. The number of rows equals the
            number of nodes in the network, and the array has 24 columns.
        """
        system = self.create_system()
        time, T_sa, T_za = self._get_input_data(T_sa, T_za, num_cycles)
        if T_n0 is None:
            resp = ct.forced_response(
                system,
                time,
                [T_sa, T_za]
            )
        else:
            resp = ct.forced_response(
                system,
                time,
                [T_sa, T_za],
                X0=T_n0.to(Units.unit_T).m
            )
        time = np.array([Q_(h, 'hr').to(Units.unit_t).m for h in range(24)])
        T_n = resp.states[:, -24:]
        return time, Q_(T_n, 'K')
    
    @staticmethod
    def _get_input_data(
        T_sa: Quantity,
        T_za: Quantity,
        num_cycles: int
    ) -> tuple[np.ndarray, ...]:
        """Internal function, helper function of `node_temperatures()`.
        Takes the daily hourly values of sol-air and zone-air temperature and
        creates arrays in which these daily values are repeated `num_cycles`
        times.
        """
        T_sa = T_sa.to(Units.unit_T).m
        T_za = T_za.to(Units.unit_T).m
        time_list, T_sa_list, T_za_list = [], [], []
        for i in range(num_cycles):
            time = np.array([
                Q_(h, 'hr').to(Units.unit_t).m
                for h in range(i * 24, (i + 1) * 24)
            ])
            time_list.append(time)
            T_sa_list.append(T_sa)
            T_za_list.append(T_za)
        time = np.concatenate(time_list)
        T_sa = np.concatenate(T_sa_list)
        T_za = np.concatenate(T_za_list)
        return time, T_sa, T_za
 