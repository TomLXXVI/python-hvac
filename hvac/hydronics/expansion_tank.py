from hvac import Quantity
from hvac.fluids import Fluid

Q_ = Quantity
Water = Fluid("water")
P_atm = Q_(101_325, 'Pa')


class ExpansionTank:
    """
    Class for sizing a diaphragm-type expansion tank in a closed-loop hydronic
    system.
    """
    def __init__(
        self,
        T_cold: Quantity,
        T_hot: Quantity,
        V_sys: Quantity,
        h_sys: Quantity,
        h_rv: Quantity,
        P_max: Quantity,
        **kwargs
    ) -> None:
        """
        Creates an `ExpansionTank` object.

        Parameters
        ----------
        T_cold:
            Initial cold water temperature when the system is filled.
        T_hot:
            Maximum operating temperature of the system.
        V_sys:
            Fluid volume in the system.
        h_sys:
            Height from the inlet of the expansion tank to top of the system.
        h_rv:
            Height from the inlet of the pressure-relief valve to the inlet of
            the expansion tank. If the inlet of the expansion tank is below the
            inlet of the pressure-relief valve, `h_rv` has a negative value.
        P_max:
            Maximum allowable system pressure, i.e. the rated pressure of the
            system's pressure-relief valve (gage pressure).
        kwargs:
            g:
                Gravity constant. Default value is 9.81 m/s2.
            P_stat_margin:
                Additional safety margin for static system pressure. Default
                value is 0.3 bar.
            V_fill_factor:
                The volume of water that the expansion tank must absorb when
                filling the system expressed as a fraction of the system volume.
                Default value is 1 %.
            P_max_margin:
                Additional safety margin applied to the maximum system pressure.
                Default value is 0.5 bar.
        """
        self.T_cold = T_cold.to('K')
        self.T_hot = T_hot.to('K')
        self.V_sys = V_sys.to('m**3')
        self.h_sys = h_sys.to('m')
        self.h_rv = h_rv.to('m')
        self.P_max = P_max.to('Pa')

        self.g = kwargs.get("g", Q_(9.81, 'm / s**2'))
        self.P_stat_margin = kwargs.get("P_stat_margin", Q_(0.3, 'bar'))
        self.V_fill_factor = kwargs.get("V_fill_factor", Q_(1.0, 'pct'))
        self.P_max_margin = kwargs.get("P_max_margin", Q_(0.5, 'bar'))

        self.water_cold = Water(T=self.T_cold, P=P_atm)
        self.water_hot = Water(T=self.T_hot, P=P_atm)
        self.V_fill = self.V_fill_factor.to('frac').m * self.V_sys

    def expansion_volume(self) -> Quantity:
        """
        Returns the volume increase due to thermal expansion of the water in
        the system when heated from its initial cold temperature to its maximum
        operating temperature.
        """
        v_cold = (1 / self.water_cold.rho).to('m**3 / kg')
        v_hot = (1 / self.water_hot.rho).to('m**3 / kg')
        dV_exp = (v_hot / v_cold - 1) * self.V_sys
        return dV_exp

    def static_system_pressure(self) -> Quantity:
        """
        Returns the static pressure of the water at the inlet of the expansion
        tank due to the height of the system above the inlet and including a
        safety margin (default margin is 0.3 bar).
        """
        rho = self.water_cold.rho.to('kg / m**3')
        g = self.g.to('m / s**2')
        P_stat = rho * g * self.h_sys + self.P_stat_margin.to('Pa')
        return P_stat.to('Pa')  # gage pressure

    def pre_air_pressurization(self) -> Quantity:
        """
        Returns the required air-side pressurization of the tank taking account
        of the static water pressure at the inlet of the tank. The air-side
        pre-pressurization cannot be less than 0.5 bar.
        """
        P_stat = self.static_system_pressure().to('bar')
        P_pre_min = Q_(0.5, 'bar')
        P_pre = max(P_pre_min, P_stat)
        return P_pre.to('Pa')  # gage pressure

    def final_air_pressurization(self) -> Quantity:
        """
        Returns the maximum pressure available at the inlet of the expansion
        tank based on the maximum allowable system pressure and the possible
        height difference between the inlet of the expansion tank and the
        pressure-relief valve of the system.
        """
        rho = self.water_cold.rho.to('kg / m**3')
        g = self.g.to('m / s**2')
        P_fin = self.P_max - self.P_max_margin.to('Pa') - rho * g * self.h_rv
        return P_fin.to('Pa')  # gage pressure

    def minimum_required_tank_volume(self) -> Quantity:
        """
        Returns the minumum required volume of the expansion tank. A
        commercially available expansion tank can then be selected whose volume
        is just larger than the calculated minimum required volume.
        """
        P1 = self.pre_air_pressurization()
        P2 = self.final_air_pressurization()
        V = self.V_fill + self.expansion_volume()
        V_tank = (P2 + P_atm) / (P2 - P1) * V
        return V_tank.to('m**3')

    def filling_pressure(self, V_tank: Quantity) -> Quantity:
        """
        Returns the gage pressure to which the system should be filled with
        water once a commercially available expansion tank has been selected.

        Parameters
        ----------
        V_tank:
            Total volume of the selected expansion tank.
        """
        P1 = self.pre_air_pressurization()
        P2 = self.final_air_pressurization()
        dV_exp = self.expansion_volume()
        V_tank = V_tank.to('m**3')
        V_fill = (P2 - P1) / (P2 + P_atm) * V_tank - dV_exp
        P_fill = V_tank / (V_tank - V_fill) * (P1 + P_atm) - P_atm
        P_fill_rec = P1 + Q_(0.3, 'bar').to('Pa')
        P_fill = max(P_fill.to('Pa'), P_fill_rec)
        return P_fill  # gage pressure
