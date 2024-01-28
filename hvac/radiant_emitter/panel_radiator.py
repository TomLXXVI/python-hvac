import numpy as np
from scipy.optimize import root_scalar
from hvac import Quantity
from hvac.fluids import Fluid

Q_ = Quantity
Water = Fluid('Water')


class PanelRadiator:
    """
    Models a panel radiator.
    """
    SYSTEM_WATER_PRESSURE: Quantity = Q_(1.5, 'bar')
    # the average pressure in the hydronic water distribution system

    def __init__(
        self,
        Qe_dot_nom: Quantity,
        Tw_sup_nom: Quantity,
        Tw_ret_nom: Quantity,
        Ti_nom: Quantity,
        n_exp: float
    ) -> None:
        """Creates a `PanelRadiator` instance from manufacturer's specifications.

        Parameters
        ----------
        Qe_dot_nom:
            Emitted heat power at nominal conditions.
        Tw_sup_nom:
            Nominal water supply temperature.
        Tw_ret_nom:
            Nominal water return temperature.
        Ti_nom:
            Nominal indoor air temperature.
        n_exp:
            PanelRadiator exponent.
        """
        Tw_avg_nom = (Tw_sup_nom.to('K') + Tw_ret_nom.to('K'))/2
        self._water = Water(
            T=Tw_avg_nom,
            P=self.SYSTEM_WATER_PRESSURE
        )
        self.Qe_dot_nom = Qe_dot_nom
        self.Tw_sup_nom = Tw_sup_nom
        self.Tw_ret_nom = Tw_ret_nom
        self.Ti_nom = Ti_nom
        self.n_exp = n_exp
        self.Vw_dot_nom = self._calculate_Vw_nom()

        self._Qe_dot_nom = Qe_dot_nom.to('W').m
        self._Tw_sup_nom = Tw_sup_nom.to('K').m
        self._Tw_ret_nom = Tw_ret_nom.to('K').m
        self._Ti_nom = Ti_nom.to('K').m
        self._UA = self._calculate_thermal_transmittance()

        self._Ti: float | None = None
        self._Vw_dot: float | None = None
        self._Tw_sup: float | None = None
        self._Qe_dot: float | None = None

    def _calculate_thermal_transmittance(self) -> float:
        """Calculates the radiator's thermal transmittance from the nominal
        specifications.
        """
        DT_max = self._Tw_sup_nom - self._Ti_nom
        DT_min = self._Tw_ret_nom - self._Ti_nom
        lmtd = (DT_max - DT_min) / np.log(DT_max / DT_min)
        UA = self._Qe_dot_nom / (lmtd ** self.n_exp)
        return UA

    def _calculate_Vw_nom(self) -> Quantity:
        """Calculates the water volume flow rate through the radiator at nominal
        conditions.
        """
        DTw_nom = self.Tw_sup_nom.to('K') - self.Tw_ret_nom.to('K')
        Vw_nom = self.Qe_dot_nom / (self._water.rho * self._water.cp * DTw_nom)
        return Vw_nom.to('m ** 3 / s')

    def _calculate_Qe_dot(self, Tw_sup: float, Vw_dot: float, Ti: float) -> float:
        """Calculates the emitted heat power of the radiator for the given
        conditions.

        Parameters
        ----------
        Tw_sup:
            Water supply temperature.
        Vw_dot:
            Water volume flow rate.
        Ti:
            Indoor air temperature.
        """
        supply_water = Water(T=Q_(Tw_sup, 'K'), P=self.SYSTEM_WATER_PRESSURE)
        rho_sw = supply_water.rho.to('kg / m ** 3').m
        cp_sw = supply_water.cp.to('J / (kg * K)').m
        DT_max = Tw_sup - Ti
        Qe_dot_max = rho_sw * cp_sw * Vw_dot * DT_max

        def eq(Qe_dot: float) -> float:
            try:
                DTw = Qe_dot / (rho_sw * cp_sw * Vw_dot)
            except ZeroDivisionError:
                Qe_dot_new = 0.0
            else:
                Tw_ret = Tw_sup - DTw
                DT_min = Tw_ret - Ti
                lmtd = (DT_max - DT_min) / np.log(DT_max / DT_min)
                Qe_dot_new = self._UA * lmtd ** self.n_exp
            return Qe_dot - Qe_dot_new

        sol = root_scalar(eq, bracket=(1e-9, Qe_dot_max - 1e-9))
        Qe = sol.root
        return Qe

    def _calculate_Vw_dot(self, Qe_dot: float, Tw_sup: float, Ti: float) -> float:
        """Calculates the needed water volume flow rate through the radiator to
        get the desired heat output when the water supply temperature and
        indoor air temperature are given.

        Parameters
        ----------
        Qe_dot:
            Desired heat power from the radiator.
        Tw_sup:
            Water supply temperature.
        Ti:
            Indoor air temperature.
        """
        supply_water = Water(T=Q_(Tw_sup, 'K'), P=self.SYSTEM_WATER_PRESSURE)
        rho_sw = supply_water.rho.to('kg / m ** 3').m
        cp_sw = supply_water.cp.to('J / (kg * K)').m
        DT_max = Tw_sup - Ti
        Vw_dot_min = Qe_dot / (rho_sw * cp_sw * DT_max)

        def eq(Vw_dot):
            DTw = Qe_dot / (rho_sw * cp_sw * Vw_dot)
            Tw_ret = Tw_sup - DTw
            DT_min = Tw_ret - Ti
            lmtd = (DT_max - DT_min) / np.log(DT_max / DT_min)
            Qe_new = self._UA * lmtd ** self.n_exp
            return Qe_dot - Qe_new

        Vw_dot_max = 1000  # m ** 3 / s
        sol = root_scalar(eq, bracket=(Vw_dot_min + 1e-9, Vw_dot_max))
        Vw_dot = sol.root
        return Vw_dot

    def _calculate_Tw_sup(self, Qe_dot: float, Vw_dot: float, Ti: float) -> float:
        """Calculates the required water supply temperature to get the desired
        heat output from the radiator when water volume flow rate and indoor
        air temperature are given.

        Parameters
        ----------
        Qe_dot:
            Desired heat output from the radiator.
        Vw_dot:
            Water volume flow rate through the radiator.
        Ti:
            Indoor air temperature.
        """
        def eq(Tw_sup):
            supply_water = Water(T=Q_(Tw_sup, 'K'), P=self.SYSTEM_WATER_PRESSURE)
            rho_sw = supply_water.rho.to('kg / m ** 3').m
            cp_sw = supply_water.cp.to('J / (kg * K)').m
            DTw = Qe_dot / (rho_sw * cp_sw * Vw_dot)
            Tw_ret = Tw_sup - DTw
            DT_max = Tw_sup - Ti
            DT_min = Tw_ret - Ti
            lmtd = (DT_max - DT_min) / np.log(DT_max / DT_min)
            Qe_dot_new = self._UA * lmtd ** self.n_exp
            return Qe_dot - Qe_dot_new

        Tw_sup_max = 363.15  # K, 90 °C
        sol = root_scalar(eq, bracket=(Ti, Tw_sup_max))
        Tw_sup = sol.root
        return Tw_sup

    def _calculate_Vw_dot_prim(
        self,
        Qe_dot: float,
        Vw_dot: float,
        Ti: float,
        Tw_sup_prim: float
    ) -> float:
        """Mixing case. Through the radiator flows a constant water volume flow
        rate `Vw_dot`. Calculates the required primary water volume flow rate to
        the mixing valve to get at the required water supply temperature for
        which the given, desired heat output `Qe_dot` from the radiator will be
        established.

        Parameters
        ----------
        Qe_dot:
            Desired heat output from the radiator.
        Vw_dot:
            Constant water volume flow rate through the radiator.
        Ti:
            Indoor air temperature.
        Tw_sup_prim:
            The constant water supply temperature in the primary circuit.
        """
        Tw_sup = self._calculate_Tw_sup(Qe_dot, Vw_dot, Ti)
        supply_water = Water(T=Q_(Tw_sup, 'K'), P=self.SYSTEM_WATER_PRESSURE)
        rho_sw = supply_water.rho.to('kg / m ** 3').m
        cp_sw = supply_water.cp.to('J / (kg * K)').m
        DTw = Qe_dot / (rho_sw * cp_sw * Vw_dot)
        Tw_ret = Tw_sup - DTw
        Vw_dot_prim = DTw / (Tw_sup_prim - Tw_ret) * Vw_dot
        return Vw_dot_prim

    def inputs(
        self,
        Ti: Quantity,
        Vw_dot: Quantity | None = None,
        Tw_sup: Quantity | None = None,
        Qe: Quantity | None = None
    ) -> None:
        """Sets the known input parameters of the `PanelRadiator` instance. Only one
        input parameter can be left to `None`. This parameter will be calculated
        and returned by method `output()`.

        Parameters
        ----------
        Ti: always required
            Indoor air temperature.
        Vw_dot: optional, default None
            Water volume flow rate through the radiator.
        Tw_sup: optional, default None
            Water supply temperature at the radiator's inlet.
        Qe: optional, default None
            The emitted heat power from the radiator.
        """
        self._Ti = Ti.to('K').m
        self._Vw_dot = Vw_dot.to('m ** 3 / s').m if Vw_dot is not None else None
        self._Tw_sup = Tw_sup.to('K').m if Tw_sup is not None else None
        self._Qe_dot = Qe.to('W').m if Qe is not None else None

    def output(self) -> Quantity:
        """Returns the output from the `PanelRadiator` instance. The output will
        depend on the inputs given in method `inputs()`. E.g. if the water
        volume flow rate `Vw_dot` was not set, the output will be the water
        volume flow rate.
        """
        if self._Vw_dot is None:
            Vw_dot = self._calculate_Vw_dot(self._Qe_dot, self._Tw_sup, self._Ti)
            return Q_(Vw_dot, 'm ** 3 / s')
        if self._Tw_sup is None:
            Tw_sup = self._calculate_Tw_sup(self._Qe_dot, self._Vw_dot, self._Ti)
            return Q_(Tw_sup, 'K')
        if self._Qe_dot is None:
            Qe = self._calculate_Qe_dot(self._Tw_sup, self._Vw_dot, self._Ti)
            return Q_(Qe, 'W')

    def __call__(
        self,
        Ti: Quantity,
        Vw_dot: Quantity | None = None,
        Tw_sup: Quantity | None = None,
        Qe_dot: Quantity | None = None
    ) -> Quantity:
        """Returns the output of the `PanelRadiator` instance depending on the
        input parameters that are set.

        Parameters
        ----------
        Ti: always required
            Indoor air temperature.
        Vw_dot: optional, default None
            The water volume flow rate through the radiator.
        Tw_sup: optional, default None.
            The water supply temperature at the radiator's inlet.
        Qe_dot: optional, default None.
            The thermal power emitted from the radiator.

        Returns
        -------
        - If input parameter `Vw_dot` is None, the required water volume flow rate
        through the radiator.
        - If input parameter `Tw_sup` is None, the required water supply
        temperature at the radiator's inlet.
        - If input parameter `Qe_dot` is None, the resulting thermal power
        emitted from the radiator.
        """
        self.inputs(Ti, Vw_dot, Tw_sup, Qe_dot)
        return self.output()

    def QV_characteristic(
        self,
        Tw_sup: Quantity,
        Ti: Quantity,
        Vw_dot_max: Quantity,
        percent: bool = False
    ) -> tuple[Quantity, Quantity]:
        """Returns the QV-characteristic of the radiator, i.e. the emitted
        thermal power from the radiator as a function of the water volume flow
        rate through the radiator when water supply temperature and indoor air
        temperature remain constant.

        Parameters
        ----------
        Tw_sup:
            Water supply temperature.
        Ti:
            Indoor air temperature.
        Vw_dot_max:
            The maximum water volume flow rate through the radiator when the
            radiator valve is fully open.
        percent: optional, default False
            Indicates if the characteristic should contain the real volume
            flow rate and heat output (this is the default), or if the
            flow rates should be expressed as a percentage of the maximum
            flow rate and the heat outputs as a percentage of the maximum
            heat output that corresponds with the maximum flow rate.

        Returns
        -------
        A tuple of which the first element is the range of water volume flow
        rates at which the heat output of the radiator is calculated. The
        second element is the range of calculated heat outputs that corresponds
        with the range of water volume flow rates.
        """
        Qe_dot_max = self(Ti=Ti, Tw_sup=Tw_sup, Vw_dot=Vw_dot_max).to('W').m
        Vw_dot_frac = np.linspace(0, 1, endpoint=True)
        Vw_dot_max = Vw_dot_max.to('m ** 3 / s').m
        Vw_dot_rng = Vw_dot_frac * Vw_dot_max
        Tw_sup = Tw_sup.to('K').m
        Ti = Ti.to('K').m
        Qe_dot_rng = [self._calculate_Qe_dot(Tw_sup, Vw_dot, Ti) for Vw_dot in Vw_dot_rng]
        if not percent:
            return Q_(Vw_dot_rng, 'm ** 3 / s'), Q_(Qe_dot_rng, 'W')
        else:
            Qe_dot_rng = [Qe / Qe_dot_max for Qe in Qe_dot_rng]
            return Q_(Vw_dot_frac, 'frac').to('pct'), Q_(Qe_dot_rng, 'frac').to('pct')

    def QT_characteristic(
        self,
        Vw_dot: Quantity,
        Ti: Quantity,
        Tw_sup_max: Quantity
    ) -> tuple[Quantity, Quantity]:
        """Returns the QT-characteristic of the radiator, i.e. the emitted
        thermal power from the radiator as a function of the water supply
        temperature when water volume flow rate and indoor air temperature
        remain constant.

        Parameters
        ----------
        Vw_dot:
            Water volume flow rate through the radiator.
        Ti:
            Indoor air temperature.
        Tw_sup_max:
            The maximum possible water supply temperature at the radiator's
            inlet.

        Returns
        -------
        A tuple of which the first element is the range of water supply
        temperatures starting at the given indoor air temperature. The
        second element is the range of calculated heat outputs that corresponds
        with the range of water supply temperatures.
        """
        Vw_dot = Vw_dot.to('m ** 3 / s').m
        Ti = Ti.to('K').m
        Tw_sup_max = Tw_sup_max.to('K').m
        Tw_sup_rng = np.linspace(Ti, Tw_sup_max, endpoint=True)
        Qe_dot_rng = [
            self._calculate_Qe_dot(Tw_sup, Vw_dot, Ti)
            for Tw_sup in Tw_sup_rng
        ]
        return Q_(Tw_sup_rng, 'K'), Q_(Qe_dot_rng, 'W')

    def Qe_dot_percent_range(
        self,
        Tw_sup: Quantity,
        Ti: Quantity,
        Vw_dot_max: Quantity,
        Vw_dot_per_rng: list[Quantity]
    ) -> tuple[Quantity, Quantity]:
        """Returns the heat output of the radiator for a range of water volume
        flow rates.

        Parameters
        ----------
        Tw_sup:
            Water supply temperature at the radiator's inlet.
        Ti:
            Indoor air temperature of space.
        Vw_dot_max:
            The maximum water volume flow rate through the radiator.
        Vw_dot_per_rng:
            A range of water volume flow rates expressed as a fraction or
            percentage of `Vw_dot_max`.

        Returns
        -------
        A tuple of which the first element is the range of heat outputs expressed
        as a percentage that correspond with the range of water volume flow rates.
        The second element is the maximum heat output that corresponds with the
        given maximum water volume flow rate `Vw_dot_max`.
        """
        Vw_dot_max = Vw_dot_max.to('m ** 3 / s').m
        Vw_dot_rng = [
            Vw_dot_frac.to('frac').m * Vw_dot_max
            for Vw_dot_frac in Vw_dot_per_rng
        ]
        Tw_sup = Tw_sup.to('K').m
        Ti = Ti.to('K').m
        Qe_dot_rng = [
            self._calculate_Qe_dot(Tw_sup, Vw_dot, Ti)
            for Vw_dot in Vw_dot_rng
        ]
        Qe_dot_max = self._calculate_Qe_dot(Tw_sup, Vw_dot_max, Ti)
        Qe_dot_per_rng = [
            Q_(Qe / Qe_dot_max, 'pct')
            for Qe in Qe_dot_rng
        ]
        return Qe_dot_per_rng, Q_(Qe_dot_max, 'W')
