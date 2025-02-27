from __future__ import annotations

import typing
import control as ct
from hvac import Quantity
from hvac.fluids.fluid_experimental import PureFluid
from ..building.fenestration import Window
from .thermal_network import (
    LinearThermalNetwork, 
    TemperatureNode, 
    Capacitor,
    Resistor,
    ConnectionSide,
    Connection
)
from .ext_build_elem_model import ExteriorBuildingElementModel
from .units import Units

if typing.TYPE_CHECKING:
    from ..building.building_elements import (
        ExteriorBuildingElement,
        InteriorBuildingElement
    )
    from ..ventilation import ThermalZoneVentilation


Q_ = Quantity
Air = PureFluid('Air')


class ThermalZoneModel(LinearThermalNetwork):
    """Represents a temperature zone modeled as a linear thermal network.

    The linear thermal network has two temperature nodes: the interior thermal
    mass node and the zone-air node. The interior thermal mass node represents
    the interior thermal mass in the zone which acts as a thermal storage 
    buffer. The buffer is loaded by radiative heat gains in the zone.
    Between the interior thermal mass node and the zone air node is a coupling
    thermal resistor through which heat can be transferred between these two
    nodes. Convective heat gains in the zone are transferred to the zone-air
    node directly. Also, heat release to, or heat extraction from the zone air 
    by a HVAC-system acts on the zone-air node.
    """
    # The input variables of the zone system need a name to uniquely identify 
    # them, so that we can assign the right values to the right inputs when 
    # simulating the zone system. `_input_names` defines the names of the
    # input variables that exist at the level of the thermal zone only.
    input_names = {
        'T_ext_db': 'T_ext_db',  # refers to outdoor air dry-bulb temperature
        'Q_gain_rad': 'Q_gain_rad',  # refers to radiative heat gain in the zone
        'Q_gain_conv': 'Q_gain_conv',  # refers to conductive heat gain in the zone
        'Q_hvac': 'Q_hvac',  # refers to heat input from the HVAC-system
        'T_sp': 'T_sp',  # refers to the setpoint temperature for the zone
        'T_sup': 'T_sup',  # refers to temperature of ventilation supply air
    }
    
    def __init__(
        self,
        name: str,
        floor_area: Quantity,
        ceiling_height: Quantity,
        unit_R_im: Quantity,
        unit_C_im: Quantity,
        factor_surf_im: float,
        **kwargs
    ) -> None:
        """Creates a `ZoneModel` instance.
        
        Parameters
        ----------
        name:
            Name of the temperature zone.
        floor_area:
            Floor area of the zone.
        ceiling_height:
            Ceiling height in the zone.
        unit_R_im:
            Unit thermal resistance of the interior mass surface. 
        unit_C_im:
            Unit thermal capacity of the interior mass, i.e. per unit area.
        factor_surf_im:
            Multiplication factor for the interior mass surface area. 
        
        Other Parameters
        ----------------
        air_temperature:
            Reference zone air temperature used for calculating the thermal
            capacity of the zone air. Default value is 20 °C.
        air_pressure:
            Reference zone air pressure used for calculating the thermal 
            capacity of the zone air. Default value is 101,325 Pa.
        input_names:
            Identifiers for the inputs of the zone system.
        """
        self.name = name
        self.floor_area = floor_area
        self.ceiling_height = ceiling_height
        self.unit_R_m = unit_R_im
        self.unit_C_m = unit_C_im
        self.factor_surf_im = factor_surf_im
        self.air_temperature = kwargs.get('air_temperature', Q_(20, 'degC'))
        self.air_pressure = kwargs.get('air_pressure', Q_(101_325, 'Pa'))
        # Create the "bare" zone model with an interior mass node and a zone-air
        # node.
        self.int_mass_node, self.zone_air_node = self._create_nodes()
        super().__init__(self.name, (self.int_mass_node, self.zone_air_node))
   
    def _create_nodes(self) -> tuple[TemperatureNode, ...]:
        """Creates and returns the linear thermal network model of the 
        temperature zone which consists of a zone-air node and a thermal storage
        node.
        """
        int_mass_node = self._create_interior_mass_node()
        zone_air_node = self._create_zone_air_node()
        R_im = self.unit_R_m / (self.factor_surf_im * self.floor_area)
        R_im = Resistor(R_im)
        int_mass_node.connect_resistor(
            other_terminal=zone_air_node,
            resistor=R_im,
            this_side=ConnectionSide.OUT
        )
        int_mass_node.add_heat_flow(
            name=self.input_names['Q_gain_rad'],
            side=ConnectionSide.IN
        )
        zone_air_node.connect_resistor(
            other_terminal=int_mass_node,
            resistor=R_im,
            this_side=ConnectionSide.IN
        )
        zone_air_node.add_heat_flow(
            name=self.input_names['Q_gain_conv'],
            side=ConnectionSide.IN
        )
        zone_air_node.add_heat_flow(
            name=self.input_names['Q_hvac'],
            side=ConnectionSide.IN
        )
        return int_mass_node, zone_air_node

    def _create_zone_air_node(self) -> TemperatureNode:
        air = Air(T=self.air_temperature, P=self.air_pressure)
        volume = self.floor_area * self.ceiling_height
        C_za = air.rho * air.cp * volume
        capacitor_za = Capacitor(C_za)
        zone_air_node = TemperatureNode(
            name='T_za',
            network_name=self.name,
            capacitor=capacitor_za
        )
        return zone_air_node

    def _create_interior_mass_node(self) -> TemperatureNode:
        C_im = self.factor_surf_im * self.floor_area * self.unit_C_m
        capacitor_im = Capacitor(C_im)
        int_mass_node = TemperatureNode(
            name='T_im',
            network_name=self.name,
            capacitor=capacitor_im
        )
        return int_mass_node
    
    def assemble(
        self, 
        ext_build_elem_models: tuple[ExteriorBuildingElementModel, ...],
        int_build_elems: tuple[InteriorBuildingElement, ...] | None = None,
        windows: tuple[Window, ...] | None = None,
        exterior_doors: tuple[ExteriorBuildingElement, ...] | None = None,
        interior_doors: tuple[InteriorBuildingElement, ...] | None = None,
        ventilation: ThermalZoneVentilation | None = None
    ) -> ThermalZoneModel:
        """Connects the linear thermal network model of the zone with the
        linear thermal network models of the exterior building elements that
        surround the zone.
        
        Parameters
        ----------
        ext_build_elem_models:
            Linear thermal network models of the exterior building elements
            surrounding the zone.
        int_build_elems:
            Interior building elements surrounding the zone.
        windows:
            Windows surrounding the zone.
        exterior_doors:
            Exterior doors surrounding the zone.
        interior_doors:
            Interior doors surrounding the zone.
        ventilation:
            The ventilation of the zone.
        
        Returns
        -------
        ThermalZoneModel
            The resulting, new linear thermal network of the "assembled zone".
        """
        # Connect the thermal networks of exterior building elements with the 
        # thermal network of the zone. `connect()` returns a new object of 
        # class `LinearThermalNetwork`.
        connections = self._create_connections(ext_build_elem_models)
        zone_network = self.connect(connections, self.name)
        # Re-create the `ZoneModel` object and assign it the connected nodes of 
        # `zone_network` (which is only a basic `LinearThermalNetwork` object).
        zone_model = ThermalZoneModel(
            name=self.name,
            floor_area=self.floor_area,
            ceiling_height=self.ceiling_height,
            unit_R_im=self.unit_R_m,
            unit_C_im=self.unit_C_m,
            factor_surf_im=self.factor_surf_im,
            air_temperature=self.air_temperature,
            air_pressure=self.air_pressure,
            input_names=self.input_names
        )
        zone_model.nodes = zone_network.nodes
        zone_model.int_mass_node = self.int_mass_node
        zone_model.zone_air_node = self.zone_air_node
        # Connect windows to the zone model.
        if windows is not None:
            # Convective conduction heat gain through windows.
            UA_wnds_conv = sum((1 - window.F_rad) * window.UA for window in windows)
            R_wnds_conv = Resistor(1 / UA_wnds_conv)
            zone_model.zone_air_node.connect_resistor(
                other_terminal=self.input_names['T_ext_db'],
                resistor=R_wnds_conv,
                this_side=ConnectionSide.IN
            )
            # Radiative conduction heat gain through windows.
            UA_wnds_rad = sum(window.F_rad * window.UA for window in windows)
            R_wnds_rad = Resistor(1 / UA_wnds_rad)
            zone_model.int_mass_node.connect_resistor(
                other_terminal=self.input_names['T_ext_db'],
                resistor=R_wnds_rad,
                this_side=ConnectionSide.IN
            )
        # Connect exterior doors to the zone model.
        if exterior_doors is not None:
            for ext_door in exterior_doors:
                R_rad = Resistor(1 / (ext_door.F_rad * ext_door.UA))
                R_conv = Resistor(1 / ((1 - ext_door.F_rad) * ext_door.UA))
                # Radiative conduction heat gain through exterior door.
                zone_model.int_mass_node.connect_resistor(
                    other_terminal=ext_door.get_input_name(),
                    resistor=R_rad,
                    this_side=ConnectionSide.IN
                )
                # Convective conduction heat gain through exterior door.
                zone_model.zone_air_node.connect_resistor(
                    other_terminal=ext_door.get_input_name(),
                    resistor=R_conv,
                    this_side=ConnectionSide.IN
                )
        # Connect interior building elements to the zone model.
        if int_build_elems is not None:
            for int_build_elem in int_build_elems:
                R_rad = Resistor(1 / (int_build_elem.F_rad * int_build_elem.UA))
                R_conv = Resistor(1 / ((1 - int_build_elem.F_rad) * int_build_elem.UA))
                # Radiative conduction heat gain through interior building element
                zone_model.int_mass_node.connect_resistor(
                    other_terminal=int_build_elem.get_input_name(),
                    resistor=R_rad,
                    this_side=ConnectionSide.IN
                )
                # Convective conduction heat gain through interior building element
                zone_model.zone_air_node.connect_resistor(
                    other_terminal=int_build_elem.get_input_name(),
                    resistor=R_conv,
                    this_side=ConnectionSide.IN
                )
        # Connect interior doors to the zone model.
        if interior_doors is not None:
            for int_door in interior_doors:
                R_rad = Resistor(1 / (int_door.F_rad * int_door.UA))
                R_conv = Resistor(1 / ((1 - int_door.F_rad) * int_door.UA))
                # Radiative conduction heat gain through interior door
                zone_model.int_mass_node.connect_resistor(
                    other_terminal=int_door.get_input_name(),
                    resistor=R_rad,
                    this_side=ConnectionSide.IN
                )
                # Convective conduction heat gain through interior door
                zone_model.zone_air_node.connect_resistor(
                    other_terminal=int_door.get_input_name(),
                    resistor=R_conv,
                    this_side=ConnectionSide.IN
                )
        # Connect ventilation to the zone model.
        if ventilation:
            air = Air(T=self.air_temperature, P=self.air_pressure)
            # Heat gain associated with outdoor air infiltration. 
            if ventilation.V_dot_ext > 0.0:
                R_inf = Resistor(1 / (air.rho * air.cp * ventilation.V_dot_ext))
                zone_model.zone_air_node.connect_resistor(
                    other_terminal=self.input_names['T_ext_db'],
                    resistor=R_inf,
                    this_side=ConnectionSide.IN
                )
            # Heat gain associated with ventilation supply air.
            if ventilation.V_dot_sup > 0.0:
                R_sup = Resistor(1 / (air.rho * air.cp * ventilation.V_dot_sup))
                zone_model.zone_air_node.connect_resistor(
                    other_terminal=self.input_names['T_sup'],
                    resistor=R_sup,
                    this_side=ConnectionSide.IN
                )
            # Heat gain associated with air transfer from adjacent zones.
            if ventilation.V_dot_trf_dict is not None:
                int_build_elem_dict = {
                    ibe.adj_zone_name: ibe.get_input_name() 
                    for ibe in int_build_elems
                }
                for adj_zone_name, V_dot_trf in ventilation.V_dot_trf_dict.items():
                    R_trf = Resistor(1 / (air.rho * air.cp * V_dot_trf))
                    zone_model.zone_air_node.connect_resistor(
                        other_terminal=int_build_elem_dict[adj_zone_name],
                        resistor=R_trf,
                        this_side=ConnectionSide.IN
                    )
        return zone_model

    def _create_connections(
        self,
        ext_build_elem_models: tuple[ExteriorBuildingElementModel, ...]
    ) -> list[Connection]:
        connections = []
        for ext_build_elem_model in ext_build_elem_models:
            conn_im = Connection(
                node=self.int_mass_node,
                other_network=ext_build_elem_model,
                other_node=ext_build_elem_model.nodes[-1],  # interior surface node
                resistor=ext_build_elem_model.R_rad,
                other_side=ConnectionSide.OUT
            )
            connections.append(conn_im)
            conn_za = Connection(
                node=self.zone_air_node,
                other_network=ext_build_elem_model,
                other_node=ext_build_elem_model.nodes[-1],  # interior surface node
                resistor=ext_build_elem_model.R_conv,
                other_side=ConnectionSide.OUT
            )
            connections.append(conn_za)
        return connections
        
    def create_feedback_system(
        self, 
        K_p: Quantity,
        K_i: Quantity,
        reduced_order: int | None = None
    ) -> ct.StateSpace:
        """Creates and returns a control feedback system incorporating the zone.
        The feedback system is used to control the zone-air temperature to the 
        desired setpoint temperature. The temperature controller (a P- or 
        PI-controller) controls the heat output of the HVAC-system to the zone.
        
        Parameters
        ----------
        K_p:
            Proportional gain of the temperature controller.
        K_i:
            Integral gain of the temperature controller.
        reduced_order:
            Reduces the order of the feedback system to the indicated order.
            If left `None`, the order of the feedback system is not reduced.

        Returns
        -------
        control.StateSpace
            State-space representation of the feedback system incorporating the
            zone.
        """
        # Create the zone system from its linear thermal network model. Outputs
        # of the zone system are the temperatures of the zone-air node, the
        # interior mass node, and the interior surface nodes of exterior 
        # building elements surrounding the zone. 
        output_nodes = [self.zone_air_node, self.int_mass_node]
        output_nodes.extend([node for node in self.nodes if '_is@' in node.name])
        zone_system = super().create_system(tuple(output_nodes))
        # Create summing junction of the feedback system.
        sum_block = ct.summing_junction(
            name='sum',
            inputs=['T_sp', 'T_za'],
            output='error'  # error signal at input of feedback controller
        )
        # Create feedback controller.
        controller = self._create_PI_controller(K_p, K_i)
        # Collect the input names of the zone system. The controlled input of 
        # the zone system is the heat released into/extracted from the zone by 
        # the HVAC-system. The controlled input of the zone is connected with 
        # the output of the controller. The other inputs of the zone system are 
        # the sol-air temperatures at exterior building elements surrounding the
        # zone, the outdoor air dry-bulb temperature, and the convective and 
        # radiative heat gains in the zone. 
        controlled_input = ''
        other_inputs = []
        for input_label in zone_system.input_labels:
            if input_label.startswith(self.input_names['Q_hvac']):
                controlled_input = input_label
            else:
                other_inputs.append(input_label)
        # Define the inputs of the feedback system. These are the setpoint 
        # temperature 'T_sp' at the summing junction and the `other_inputs` of 
        # the zone system.
        inputs = [f"T_sp@{self.name}"] + other_inputs
        inplist = ["sum.T_sp"] + [f"{self.name}.{inp}" for inp in other_inputs]
        # Define the outputs of the feedback system. 
        # Note that we can also monitor internal signals of the feedback system 
        # by adding them to the outputs.
        outputs = [f"{out}" for out in zone_system.output_labels]
        outputs.extend(["Q_hvac", "error"])
        outlist = [f"{self.name}.{out}" for out in zone_system.output_labels]
        outlist.extend([f"{controller.name}.command", "sum.error"])
        # Create the feedback system by interconnecting the summing junction,
        # the controller, and the zone system.
        feedback_system = ct.interconnect(
            [zone_system, controller, sum_block],
            connections=[
                [f'{self.name}.{controlled_input}', f'{controller.name}.command'],
                ['sum.T_za', f'-{outlist[0]}'],  # note negative feedback
                [f'{controller.name}.error', 'sum.error'],
            ],
            inplist=inplist,
            outlist=outlist,
            name=self.name,
            inputs=inputs,
            outputs=outputs,
            debug=False
        )
        # Reduce the order of the feedback system. The order cannot be reduced
        # below the number of system outputs. Also, reduction is only possible 
        # if the number of system state variabels is greater than the given
        # reduced order.
        if reduced_order is not None: 
            reduced_order = max(reduced_order, len(outputs))
            if feedback_system.nstates > reduced_order:
                red_system = ct.balanced_reduction(feedback_system, reduced_order, method='matchdc')
                red_system.name = self.name
                red_system.set_inputs(inputs)
                red_system.set_outputs(outputs)
                return red_system
        return feedback_system

    @staticmethod
    def _create_PI_controller(
        K_p: Quantity,
        K_i: Quantity
    ) -> ct.StateSpace:
        K_p = K_p.to(f"{Units.unit_Q} / {Units.unit_T}").m
        K_i = K_i.to(f"{Units.unit_Q} / ({Units.unit_T} * {Units.unit_t})").m
        if K_i != 0.0:
            num = [K_p, K_i]
            den = [1, 0]
            PI_controller = ct.tf2ss(
                num, den,
                inputs='error',
                outputs='command',
                name='PI-control'
            )
            return PI_controller
        else:
            num = [K_p]
            den = [1]
            P_controller = ct.tf2ss(
                num, den,
                inputs='error',
                outputs='command',
                name='P-control'
            )
            return P_controller
