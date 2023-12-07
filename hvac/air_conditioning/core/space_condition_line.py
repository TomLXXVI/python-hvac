from typing import Optional
from hvac import Quantity
from hvac.fluids import HumidAir
from hvac.air_conditioning.core.abc_process import Equation, Variable
from hvac.charts.psychrometric_chart import StatePoint


Q_ = Quantity


class SpaceConditionLine:

    def __init__(
        self,
        space_air: HumidAir,
        Q_sen: Optional[Quantity] = None,
        Q_lat: Optional[Quantity] = None,
        Q: Optional[Quantity] = None,
        SHR: Optional[Quantity] = None
    ):
        """
        Define space condition line.

        Parameters
        ----------
        space_air: HumidAir
            condition of space air
        Q_sen: Quantity, optional
            sensible load of the space
        Q_lat: Quantity, optional
            latent load of the space
        Q: Quantity, optional
            total load of the space
        SHR: Quantity, optional
            sensible heat ratio of the load
        """
        self.space_air = space_air

        if (SHR is None) and (Q_lat is not None):
            self.SHR = Q_sen / (Q_sen + Q_lat)
        elif (SHR is None) and (Q is not None):
            self.SHR = Q_sen / Q_lat
        else:
            self.SHR = SHR

        self.c_pa = Q_(1.02e3, 'J / (kg * K)')
        self.h_wg = Q_(2555e3, 'J / kg')

        # Slope of space condition line:
        self.K = self.c_pa / self.h_wg * (1 / self.SHR - 1)

        self.equation = Equation(
            variables=[
                Variable(
                    name='K',
                    value=self.K.m,
                    unit=self.K.units
                ),
                Variable(
                    name='T_ao',
                    value=self.space_air.Tdb.m,
                    unit=self.space_air.Tdb.units
                ),
                Variable(
                    name='W_ao',
                    value=self.space_air.W.m,
                    unit=self.space_air.W.units
                ),
                Variable(
                    name='T_ai',
                    unit='K'
                ),
                Variable(
                    name='W_ai',
                    unit='kg / kg'
                )
            ],
            lhs="K * (T_ao - T_ai) - (W_ao - W_ai)"
        )

    def W_ai(self, T_ai: Quantity) -> Quantity:
        """Solve space condition line for the humidity ratio `W_ai` if temperature
        `T_ai` of supply air to the space is given, so that the desired space
        air state is maintained under the given sensible and latent space load.
        """
        self.equation['T_ai'] = T_ai
        self.equation['W_ai'] = None
        W_ai = self.equation.solve()
        return W_ai.value

    def T_ai(self, W_ai: Quantity) -> Quantity:
        """Solve space condition line for temperature `T_ai` if the humidity ratio
        `W_ai` of supply air-to-air is given, so that the desired space
        air state is maintained under the given sensible and latent space load.
        """
        self.equation['W_ai'] = W_ai
        self.equation['T_ai'] = None
        T_ai = self.equation.solve()
        return T_ai.value

    def start_point(self, dT_db: Quantity = Q_(50, 'delta_degC')) -> StatePoint:
        """
        Return start point of space condition line on the psychrometric chart.

        Parameters
        ----------
        dT_db: Quantity
            difference between dry-bulb temperature of start point and the
            space air dry-bulb temperature

        Returns
        -------
        StatePoint
        """
        T_ai_start = self.space_air.Tdb - dT_db
        W_ai_start = self.W_ai(T_ai_start)
        return StatePoint(T_ai_start, W_ai_start)

    def end_point(self, dT_db: Quantity = Q_(50, 'delta_degC')) -> StatePoint:
        """
        Return end point of space condition line on the psychrometric chart.

        Parameters
        ----------
        dT_db: Quantity
            difference between dry-bulb temperature of end point and the space
            air dry-bulb temperature

        Returns
        -------
        StatePoint
        """
        T_ai_end = self.space_air.Tdb + dT_db
        W_ai_end = self.W_ai(T_ai_end)
        return StatePoint(T_ai_end, W_ai_end)

    def space_point(self) -> StatePoint:
        """
        Return the state point of the space air on the psychrometric chart.

        Returns
        -------
        StatePoint
        """
        return StatePoint(self.space_air.Tdb, self.space_air.W)
