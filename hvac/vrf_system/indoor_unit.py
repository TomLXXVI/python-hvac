from scipy import interpolate
from hvac import Quantity


Q_ = Quantity


class IndoorUnit:
    default_units = {'temperature': 'degC'}

    def __init__(
        self,
        name: str = '',
        type_: str = '',
        model_size: int = 0,
        Q_rated: Quantity = Q_(0.0, 'kW'),
        iu_temp_corr: list[tuple[float, float]] | None = None,
        T_unit: str | None = None
    ):
        """Creates an `IndoorUnit` instance.

        Parameters
        ----------
        name:
            Descriptive and unique name for the indoor unit (e.g. the name of
            the space where it is located, followed by a rank number if there
            are more than one unit in the space).
        type_:
            The type code of the indoor unit.
        model_size:
            The model size of the indoor unit.
        Q_rated:
            The capacity of the indoor unit at rated conditions.
        iu_temp_corr:
            List of 2-elements tuples (T_indoor_air, Q_ratio). Each element is
            the coordinate of a point on the curve that gives the capacity ratio
            of the indoor unit at a given indoor air temperature (dry-bulb when
            heating, wet-bulb when cooling, and without units).
            The capacity ratio is the ratio of the actual capacity to the rated
            capacity.
        T_unit: optional, default None
            The unit of the indoor air temperatures in the list with coordinates.
            If none, the default unit degC is assumed.
        """
        self.name = name
        self.type = type_
        self.model_size = model_size
        self.Q_rated = Q_rated
        if T_unit is not None:
            self.default_units['temperature'] = T_unit

        Tia_lst, Qr_list = zip(*iu_temp_corr)
        self.T_corr_fun = interpolate.interp1d(
            x=Tia_lst,
            y=Qr_list,
            bounds_error=False,
            fill_value='extrapolate'
        )

    def get_capacity(self, Tia: Quantity) -> Quantity:
        """Get the temperature corrected indoor unit capacity.

        Parameters
        ----------
        Tia:
            Indoor air temperature. Dry-bulb when heating, wet-bulb when cooling.

        Returns
        -------
        Quantity
        The actual indoor unit capacity at the given indoor air temperature.
        """
        Tia = Tia.to(self.default_units['temperature']).m
        cf = self.T_corr_fun(Tia)
        Q_corr = cf * self.Q_rated
        return Q_corr


def get_Tia_avg(
    Ql_zones: list[Quantity],
    Tia_zones: list[Quantity]
) -> Quantity:
    """Get the load-weighted average indoor air temperature in a building.

    Parameters
    ----------
    Ql_zones:
        A list with the (rated) zone loads of the building.
    Tia_zones:
        A list with the (design) indoor air temperatures of the zones (wet-bulb
        when cooling, dry-bulb when heating).

    Returns
    -------
    The load-weighted average indoor air temperature (wet-bulb when cooling,
    dry-bulb when heating).
    """
    Ql_tot = sum(Ql_zones)
    Tia_avg = sum(
        Tia * Ql_zone
        for Tia, Ql_zone in zip(Tia_zones, Ql_zones)
    ) / Ql_tot
    return Tia_avg
