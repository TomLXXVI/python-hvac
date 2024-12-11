"""Implementation of the detailed method according to standard EN 16282-1 for
the calculation of the design extraction air volume flow rate in commercial
kitchens.
"""
from enum import Enum
from hvac import Quantity

Q_ = Quantity


class KitchenAppliance:

    def __init__(
        self,
        name: str,
        P: Quantity,
        q_dot_sen: Quantity,
        m_dot_w: Quantity,
        h_d: Quantity,
        width: Quantity | None = None,
        length: Quantity | None = None,
        f_sen_conv: Quantity = Q_(0.5, 'frac')
    ) -> None:
        """Creates a `KitchenAppliance` object.

        Parameters
        ----------
        name:
            For identifying the kitchen appliance.
        P:
            The rated power consumption of the kitchen appliance.
        q_dot_sen:
            Specific sensible heat emission of the kitchen appliance (see
            EN 16282-1, annex A.1).
        m_dot_w:
            Specific moisture emission of the kitchen appliance (see EN 16282-1,
            annex A.1).
        h_d:
            Free height between the top of the cooking appliance and the bottom
            plane of the extraction hood. In case of a kitchen appliance with a
            vertical opening, e.g. salamander or oven, the free height is
            measured from the middle of the opening height.
            In the case of a kitchen appliance that is not underneath a hood
            (i.e. unhooded), or in the case of a kitchen with an air extraction
            ceiling, the free height is the distance between the top of the
            appliance and a height of 2.5 m.
        width:
            Width of the kitchen appliance. If the appliance belongs to a
            block with other appliances close together, the width of the
            appliance can be ignored as the total width of the block will be
            used to calculate the required extraction airflow of the block.
        length:
            Length of the kitchen appliance. If the appliance belongs to a
            block with other appliances close together, the length of the
            appliance can be ignored as the total length of the block will be
            used to calculate the required extraction airflow of the block.
        f_sen_conv:
            The convective fraction of the sensible heat emission.
        """
        self.name = name
        self.P = P.to('kW')
        self.q_dot_sen = q_dot_sen.to('W / kW')
        self.m_dot_w = m_dot_w.to('g / (h * kW)')
        self.h_d = h_d.to('m')
        self.width = width.to('m') if isinstance(width, Quantity) else None
        self.length = length.to('m') if isinstance(length, Quantity) else None
        self.f_sen_conv = f_sen_conv.to('frac')
        self.Q_dot_sen_conv = self._calc_convective_sensible_load()
        self.M_dot_w = self._calc_moisture_emission()

    def _calc_convective_sensible_load(self) -> Quantity:
        """The convectively transmitted proportion of the sensible heat load
        of the kitchen appliance.
        """
        Q_dot_sen_conv = self.f_sen_conv * self.q_dot_sen * self.P
        return Q_dot_sen_conv.to('W')

    def _calc_moisture_emission(self) -> Quantity:
        """Moisture released from the appliance to the kitchen space."""
        M_dot_w = self.m_dot_w * self.P
        return M_dot_w.to('g / h')

    def _calc_thermally_induced_airflow(
        self,
        r: float,
        f_simul: float
    ) -> Quantity:
        """The non-isothermal free stream above the kitchen appliance, which
        induces air from the environment.
        """
        k = 18  # m ** (4 / 3) / (W ** 3 * h)
        d_hyd = 2 * self.length * self.width / (self.length + self.width)
        V_dot_therm = k * r * (self.Q_dot_sen_conv.m * f_simul) ** (1 / 3)
        V_dot_therm *= (self.h_d.m + 1.7 * d_hyd.m) ** (5 / 3)
        return Q_(V_dot_therm, 'm ** 3 / h')


class KitchenBlockArrangement(Enum):
    """Depending on the arrangement of the kitchen block, a certain reduction
    factor may be applied in the calculation of the thermally induced air flow.
    """
    ANYWHERE = 1.0
    AGAINST_WALL = 0.63


class SupplyAirFlowType(Enum):
    """Depending on the type of supply air diffusion in the kitchen a correction
    factor is used to increase the required extraction airflow at the hood of a
    cooking block in order to prevent flush-out of cooking plumes from the hood
    (improve containment). This correction factor integrates the influence of
    various air velocities and turbulences in the kitchen.
    """
    MIXED_FLOW_RADIAL = 1.25
    MIXED_FLOW_PLANE = 1.20
    DISPLACEMENT_FLOW_CEILING = 1.10
    DISPLACEMENT_FLOW_WALL = 1.05


class KitchenBlock:

    def __init__(
        self,
        name: str,
        width: Quantity | None = None,
        length: Quantity | None = None,
        arrangement: KitchenBlockArrangement = KitchenBlockArrangement.ANYWHERE,
        hooded: bool = True,
        V_dot_sup_dir: Quantity = Q_(0.0, 'm ** 3 / h')
    ) -> None:
        """Creates a `KitchenBlock` object, which can contain a single kitchen
        appliance or a group of multiple kitchen appliances.

        Parameters
        ----------
        name:
            For identifying the kitchen block.
        width:
            Width of the kitchen block surface if the appliances within the
            block are close together. Otherwise, the width of each individual
            appliance in the block will be considered in the calculation of the
            extraction airflow; in that case the width of each appliance in the
            block must have been specified and the width of the block will be
            ignored, i.e. the width of the block can be left to the default
            `None`.
        length:
            Length of the kitchen block surface if the appliances within the
            block are close together. Otherwise, the length of each individual
            appliance in the block will be considered in the calculation of the
            extraction airflow; in that case the length of each appliance in the
            block must have been specified and the length of the block will be
            ignored, i.e. the length of the block can be left to the default
            `None`.
        arrangement:
            See docstring of `KitchenBlockArrangement`. The kitchen block is
            either installed anywhere, or installed against a wall.
        hooded:
            Indicates whether the kitchen block is underneath an extraction hood
            or not.
        V_dot_sup_dir:
            Supply air blown directly into the hood in case of an extraction
            system with integrated supply air.

        Notes
        -----
        Supply air blown directly in the hood shall be temperature-controlled
        to avoid additional condensation and shall not exceed 15 % to 20 %.
        """
        self.name = name
        self.width = width.to('m') if isinstance(width, Quantity) else None
        self.length = length.to('m') if isinstance(length, Quantity) else None
        self.arrangement = arrangement
        self.hooded = hooded
        self.V_dot_sup_dir = V_dot_sup_dir.to('m ** 3 / h')
        self.kitchen_appliances: list[KitchenAppliance] = []
        self.f_simul: float = 1.0
        self.V_dot_ext: Quantity | None = None
        self.V_dot_therm_ext: Quantity | None = None
        self.V_dot_w_ext: Quantity | None = None

    def add_kitchen_appliances(self, *kitchen_appliances) -> None:
        """Adds a single or multiple kitchen appliances to a kitchen block."""
        self.kitchen_appliances.extend(kitchen_appliances)

    def _calc_thermally_induced_airflow(self) -> Quantity:
        """The non-isothermal free stream above the kitchen appliance, which
        induces air from the environment.
        """
        k = 18  # m ** (4 / 3) / (W ** 3 * h)
        if self.width is not None and self.length is not None:
            # If the width and length of the kitchen block are specified, it is
            # assumed that the kitchen appliances are close together and a
            # single thermal airflow for the block is calculated.
            d_hyd = 2 * self.length * self.width / (self.length + self.width)
            h_d_sum = sum(app.h_d.m for app in self.kitchen_appliances)
            num_cooking_apps = len(self.kitchen_appliances)
            h_d_avg = h_d_sum / num_cooking_apps
            Q_dot_sen_conv = sum(ka.Q_dot_sen_conv.m for ka in self.kitchen_appliances)
            r = self.arrangement.value
            V_dot_therm = k * r * (Q_dot_sen_conv * self.f_simul) ** (1 / 3)
            V_dot_therm *= (h_d_avg + 1.7 * d_hyd.m) ** (5 / 3)
            return Q_(V_dot_therm, 'm ** 3 / h')
        else:
            V_dot_therm = sum(
                ka._calc_thermally_induced_airflow(
                    r=self.arrangement.value,
                    f_simul=self.f_simul
                ) for ka in self.kitchen_appliances
            )
            return V_dot_therm

    def _calc_thermal_extraction_airflow(
        self,
        a: SupplyAirFlowType,
    ) -> Quantity:
        """The required extraction airflow at the extraction hood above the
        kitchen block, but the method can also be applied to a unhooded kitchen
        block.

        Parameters
        ----------
        a:
            See docstring of `SupplyAirFlowType`.
        """
        V_dot_therm = self._calc_thermally_induced_airflow()
        V_dot_therm_ext = V_dot_therm * a.value + self.V_dot_sup_dir
        return V_dot_therm_ext

    def _calc_moisture_extraction_airflow(
        self,
        w_ext: Quantity,
        w_sup: Quantity,
        rho: Quantity
    ) -> Quantity:
        """The required extraction airflow to compensate for moisture release
        in the kitchen space and to protect against excessive condensation.

        Parameters
        ----------
        w_ext:
            Final absolute moisture content of the extraction air.
        w_sup:
            Absolute moisture content of supply air to the kitchen space.
        rho:
            Mass density of air.
        """
        w_ext.ito('g / kg')
        w_sup.ito('g / kg')
        rho.ito('kg / m ** 3')
        M_dot_w = sum(ka.M_dot_w for ka in self.kitchen_appliances)
        M_dot_w *= self.f_simul
        V_dot_w_ext = M_dot_w / ((w_ext - w_sup) * rho)
        return V_dot_w_ext.to('m ** 3 / h')

    def _calc_extraction_airflow(
        self,
        a: SupplyAirFlowType,
        w_ext: Quantity = Q_(16.5, 'g / kg'),
        w_sup: Quantity = Q_(10.5, 'g / kg'),
        rho: Quantity = Q_(1.2, 'kg / m ** 3')
    ) -> None:
        """Calculates:
        1. the required extraction airflow to compensate for the thermal load of
        the kitchen block, accessible through instance attribute `V_dot_therm_ext`
        after the calculation.
        2. the required extraction airflow to compensate for the moisture load
        of the kitchen block, accessible through instance attribute `V_dot_w_ext`
        after the calculation.
        3. the actual required extraction airflow for the kitchen block, which
        is the largest of the required extraction airflow to compensate for the
        thermal load and the required extraction airflow to compensate for the
        moisture load, accessible through instance attribute `V_dot_ext` after
        the calculation.

        Parameters
        ----------
        a:
            See docstring of `SupplyAirFlowType`.
        w_ext:
            Final absolute moisture content of the extraction air.
        w_sup:
            Absolute moisture content of supply air to the kitchen space.
        rho:
            Mass density of air.

        Notes
        -----
        If the required extraction airflow for the moisture load is greater than
        the required extraction airflow for the thermal load, this difference
        in airflow volume can be supplied directly into the hood and is added to
        instance attribute `V_dot_sup_dir`.
        """
        self.V_dot_therm_ext = self._calc_thermal_extraction_airflow(a)
        self.V_dot_w_ext = self._calc_moisture_extraction_airflow(w_ext, w_sup, rho)
        self.V_dot_ext = max(self.V_dot_therm_ext, self.V_dot_w_ext)
        if self.V_dot_w_ext > self.V_dot_therm_ext:
            self.V_dot_sup_dir += (self.V_dot_w_ext - self.V_dot_therm_ext)


class Kitchen:

    def __init__(
        self,
        a: SupplyAirFlowType,
        f_simul: float
    ) -> None:
        """Creates a `Kitchen` object, which can contain a collection of hooded
        and unhooded kitchen blocks.

        Parameters
        ----------
        a:
            See `SupplyAirFlowType`.
        f_simul:
            Simultaneity factor, i.e. the ratio of the actual power consumption
            of the block to the total rated power of the kitchen block
            appliances. This factor depends on the type and size of the kitchen
            (see EN 16282-1, annex A.2).
        """
        self.a = a
        self.f_simul = f_simul
        self.kitchen_blocks: list[KitchenBlock] = []

    def add_kitchen_blocks(self, *kitchen_blocks) -> None:
        """Adds a single or multiple kitchen blocks to the kitchen."""
        self.kitchen_blocks.extend(kitchen_blocks)
        for kb in self.kitchen_blocks: kb.f_simul = self.f_simul

    def calc_extract_airflow(self) -> tuple[Quantity, ...]:
        """The total required extraction airflow for the kitchen, taking the
        thermal load as well as the moisture load of the kitchen blocks into
        account.

        Returns
        -------
        V_dot_ext_tot:
            Total required extraction airflow for all hooded and unhooded
            kitchen blocks in the kitchen.
        V_dot_ext_unhooded:
            Required extraction airflow for unhooded kitchen blocks in the
            kitchen.
        V_dot_comp:
            Compensation airflow to be added to the required extraction airflow
            for unhooded kitchen blocks in the kitchen such that the total
            extraction airflow required for unhooded kitchen blocks in the
            kitchen would be equal to 10 % of the required extraction airflow
            for hooded kitchen blocks in the kitchen.

        Notes
        -----
        If the total required extraction airflow for unhooded kitchen blocks
        in the kitchen `V_dot_ext_unhooded` is smaller than 10 % of the total
        required extraction airflow for hooded kitchen blocks in the kitchen
        (in that case `V_dot_comp` will be greater than zero), the total
        required extraction airflow for the kitchen `V_dot_ext_tot` can be
        exhausted completely via the extraction hoods in the kitchen.
        Otherwise, the total required extraction airflow for unhooded kitchen
        blocks in the kitchen needs to be added to the general kitchen
        ventilation system.
        """
        V_dot_ext_tot = Q_(0.0, 'm ** 3 / h')
        V_dot_ext_unhooded = Q_(0.0, 'm ** 3 / h')
        for kb in self.kitchen_blocks:
            kb._calc_extraction_airflow(self.a)
            V_dot_ext_tot += kb.V_dot_ext
            if not kb.hooded:
                V_dot_ext_unhooded += kb.V_dot_ext
        V_dot_ext_hooded = V_dot_ext_tot - V_dot_ext_unhooded
        V_dot_comp = Q_(0.0, 'm ** 3 / h')
        if V_dot_ext_unhooded < 0.1 * V_dot_ext_hooded:
            # When the required extraction airflow for unhooded kitchen
            # blocks is less than 10 % of the required extraction airflow
            # for hooded kitchen blocks, a compensation airflow is assumed
            # for extracting the airflow of unhooded kitchen blocks via the
            # extraction hoods.
            V_dot_comp = 0.1 * V_dot_ext_hooded - V_dot_ext_unhooded
        V_dot_ext_tot += V_dot_comp
        return V_dot_ext_tot, V_dot_ext_unhooded, V_dot_comp
