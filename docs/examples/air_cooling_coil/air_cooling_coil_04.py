"""
In this script, we want to determine the mass flow rate of water through a given
air-to-water cooling coil that is needed to maintain a given target temperature
of the air leaving the cooling coil and this at different mass flow rates of air
through the cooling coil.
"""
from scipy import optimize
from hvac import Quantity
from hvac.fluids import Fluid, HumidAir, FluidState
from hvac.heat_exchanger.recuperator.fintube.continuous_fin import PlainFinTubeAirToWaterCounterFlowHeatExchanger
from hvac.charts import LineChart, PsychrometricChart


Q_ = Quantity
Water = Fluid('Water')
AirCoil = PlainFinTubeAirToWaterCounterFlowHeatExchanger


class Output:
    """Takes the dictionary returned from `AirCoil.rate()` and transforms it
    into an `Output` class object of which the attribute names are the keys
    of the dictionary. (See below.)
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


# Define the air-cooling coil:
air_coil = AirCoil(
    width=Q_(1.374, 'm'),
    height=Q_(1, 'm'),
    num_rows=5,
    pitch_trv=Q_(25.4, 'mm'),
    pitch_lon=Q_(22.0, 'mm'),
    d_i=Q_(8.422, 'mm'),
    d_o=Q_(10.2, 'mm'),
    t_fin=Q_(0.3302, 'mm'),
    fin_density=1 / Q_(3.175, 'mm'),
    num_circuits=None
)


# Define a set of input parameters that determine the "nominal" performance of
# the air-cooling coil (which just serves as an arbitrary reference).
class NominalConditions:
    water_in = Water(T=Q_(8, 'degC'), P=Q_(2, 'bar'))
    water_m_dot = Q_(8906.359, 'kg / hr')
    air_in = HumidAir(Tdb=Q_(26.36, 'degC'), W=Q_(10.6, 'g / kg'))
    air_m_dot = Q_(11_470.588, 'kg / hr')


def rate(
    air_in: HumidAir = NominalConditions.air_in,
    water_in: FluidState = NominalConditions.water_in,
    air_m_dot: Quantity = NominalConditions.air_m_dot,
    water_m_dot: Quantity = NominalConditions.water_m_dot,
    verbose: bool = False
) -> Output:
    """Returns the performance of the air-cooling coil for the given operating
    conditions. Without parameters, the "nominal" performance is returned.
    """
    air_coil.set_operating_conditions(
        air_in=air_in,
        water_in=water_in,
        air_m_dot=air_m_dot,
        water_m_dot=water_m_dot
    )
    d = air_coil.rate()
    if verbose:
        print(
            f"cooling coil capacity = {d['Q_dot'].to('kW'):~P.3f}",
            f"heat transfer effectiveness = {d['eps'].to('frac'):~P.2f}",
            f"leaving air = {d['air_out'].Tdb.to('degC'):~P.2f} DB, "
            f"{d['air_out'].W.to('g / kg'):~P.1f} AH, "
            f"{d['air_out'].RH.to('pct'):~P.0f} RH",
            f"leaving water temperature = {d['water_out'].T.to('degC'):~P.2f}",
            f"air-side pressure drop = {d['dP_air'].to('Pa'):~P.0f}",
            f"water-side pressure drop = {d['dP_water'].to('kPa'):~P.3f}",
            f"air-side surface condition = {d['air_surface_condition']}",
            sep='\n', end='\n\n'
        )
    return Output(**d)


def get_water_m_dot(
    target_air_out: HumidAir,
    air_in: HumidAir,
    air_m_dot: Quantity
) -> Quantity:
    """Given the state of air entering the cooling coil `air_in` and the mass
    flow rate of air `air_m_dot` through the cooling coil, finds the mass flow
    rate of water to get the state of air leaving the cooling coil specified by
    `target_air_out`.
    """
    T_air_out = target_air_out.Tdb.to('degC').m

    def _eq_(m_dot_w: float) -> float:
        m_dot_w = Q_(m_dot_w, 'kg / hr')
        output = rate(air_in=air_in, air_m_dot=air_m_dot, water_m_dot=m_dot_w)
        T_air_out_try = output.air_out.Tdb.to('degC').m
        dev = T_air_out - T_air_out_try
        return dev

    water_m_dot_max = NominalConditions.water_m_dot.to('kg / hr').m
    water_m_dot_min = 0.25 * water_m_dot_max
    try:
        sol = optimize.root_scalar(
            _eq_,
            bracket=[water_m_dot_min, water_m_dot_max]
        )
    except ValueError:
        print(f"Sorry, but no solution found.")
    else:
        return Q_(sol.root, 'kg / hr')


def main():
    # First determine the cooling coil performance at its nominal conditions
    # (see class `NominalConditions`):
    print("Nominal cooling coil performance: ")
    nom_output = rate(verbose=True)

    # The air flow rate through the cooling coil varies from 60 to 100 % of the
    # nominal flow rate.
    airflow_fractions = range(60, 110, 10)

    # Find the corresponding water mass flow rates to keep the leaving air
    # temperature constant to its nominal value:
    water_m_dot_results = []
    for aff in airflow_fractions:
        aff = round(aff / 100, 1)
        water_m_dot = get_water_m_dot(
            target_air_out=nom_output.air_out,
            air_in=NominalConditions.air_in,
            air_m_dot=aff * NominalConditions.air_m_dot
        )
        if water_m_dot:
            water_m_dot_results.append((
                aff * NominalConditions.air_m_dot,
                water_m_dot
            ))
            print(
                f"required water mass flow rate at airflow fraction {aff} = "
                f"{water_m_dot.to('kg / hr'):~P.1f}"
            )

    air_m_dot_rng, water_m_dot_rng = zip(*water_m_dot_results)

    # Get the full air state at the cooling coil outlet:
    air_out_rng = []
    for air_m_dot, water_m_dot in zip(air_m_dot_rng, water_m_dot_rng):
        output = rate(
            air_in=NominalConditions.air_in,
            air_m_dot=air_m_dot,
            water_in=NominalConditions.water_in,
            water_m_dot=water_m_dot
        )
        air_out_rng.append(output.air_out)

    # Draw a line chart of water mass flow rate versus air mass flow rate to
    # keep the leaving air temperature constant when the air mass flow rate
    # varies:
    air_m_dot_rng = [
        air_m_dot.to('kg / hr').m
        for air_m_dot in air_m_dot_rng
    ]
    water_m_dot_rng = [
        water_m_dot.to('kg / hr').m
        for water_m_dot in water_m_dot_rng
    ]
    chart = LineChart()
    chart.add_xy_data(
        label='',
        x1_values=air_m_dot_rng,
        y1_values=water_m_dot_rng,
        style_props={'marker': 'o'}
    )
    chart.x1.add_title('air mass flow rate, kg/h')
    chart.y1.add_title('required water mass flow rate, kg/h')
    chart.show()

    # Plot the states of the air leaving the cooling coil on a psychrometric
    # chart:
    psy_chart = PsychrometricChart()
    for i, air_out in enumerate(air_out_rng):
        psy_chart.plot_point(f'point {i}', air_out)
    psy_chart.show()


if __name__ == '__main__':
    main()
