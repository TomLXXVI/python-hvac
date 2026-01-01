"""
Determine the performance of a heat pump at different outdoor air temperatures
by simulation.
"""
from hvac import Quantity
from hvac.fluids import Fluid, HumidAir
from hvac.heat_exchanger.recuperator.fintube import PFT_CF_AE as Evaporator
from hvac.heat_exchanger.recuperator.fintube import PFT_CF_AC as Condenser
from hvac.vapor_compression import VariableSpeedCompressor
from hvac.vapor_compression import SS_VCM as HeatPump
from hvac.vapor_compression import Output


Q_ = Quantity
R410A = Fluid("R410A")


def setup_heat_pump(dT_sh: Quantity) -> HeatPump:

    evaporator = Evaporator(
        W_fro=Q_(625, 'mm'),
        H_fro=Q_(208, 'mm'),
        N_rows=3,
        S_trv=Q_(25.4, 'mm'),
        S_lon=Q_(22.0, 'mm'),
        D_int=Q_(4.211, 'mm'),
        D_ext=Q_(5.1, 'mm'),
        t_fin=Q_(0.3302, 'mm'),
        N_fin=1 / Q_(3.175, 'mm'),
        k_fin=Q_(237, 'W / (m * K)'),
        num_circuits=1
    )

    condenser = Condenser(
        W_fro=Q_(636, 'mm'),
        H_fro=Q_(212, 'mm'),
        N_rows=4,
        S_trv=Q_(25.4, 'mm'),
        S_lon=Q_(22.0, 'mm'),
        D_int=Q_(4.211, 'mm'),
        D_ext=Q_(5.1, 'mm'),
        t_fin=Q_(0.3302, 'mm'),
        N_fin=1 / Q_(3.175, 'mm'),
        k_fin=Q_(237, 'W / (m * K)'),
        num_circuits=1
    )

    compressor = VariableSpeedCompressor(
        coeff_file="VRJ028-K_R410A.csv",
        refrigerant=R410A,
        dT_sh=Q_(8, 'K'),
        dT_sc=Q_(0, 'K'),
        units={"n": "rps"}
    )

    heat_pump = HeatPump(
        evaporator, condenser, compressor,
        dT_sh=dT_sh,
        n_cmp_min=Q_(900, 'rpm'),
        n_cmp_max=Q_(4200, 'rpm')
    )

    return heat_pump


def rate_heat_pump(T_ext: Quantity) -> Output:
    dT_sh = Q_(8, 'K')
    heat_pump = setup_heat_pump(dT_sh=dT_sh)
    Tdb_ext = T_ext.to('K')
    Twb_ext = Tdb_ext - Q_(1, 'K')
    T_evp_ini = Tdb_ext - dT_sh.to('K') - Q_(2, 'K')

    res = heat_pump.rate(
        evp_air_in=HumidAir(Tdb=Tdb_ext, Twb=Twb_ext),
        evp_air_m_dot=Q_(1171.439, 'kg / hr'),
        cnd_air_in=HumidAir(Tdb=Q_(20, 'degC'), RH=Q_(50, 'pct')),
        cnd_air_m_dot=Q_(1157.631, 'kg / hr'),
        n_cmp=Q_(1358, 'rpm'),
        T_evp_ini=T_evp_ini,
        T_cnd_ini=Q_(35.0, 'degC')
    )
    return res


if __name__ == '__main__':
    import shelve
    from concurrent.futures import ProcessPoolExecutor

    # Determine the heat pump performance at the specified outdoor air
    # temperatures:
    rng_T_ext = [Q_(12, 'degC'), Q_(7, 'degC'), Q_(2, 'degC'), Q_(-3, 'degC')]

    outputs = []
    with ProcessPoolExecutor() as executor:
        for output in executor.map(rate_heat_pump, rng_T_ext):
            outputs.append(output)

    # Save the simulation outputs to a shelf on disk:
    path = "./heat_pump_rating"
    with shelve.open(path) as shelf:
        for T_ext, output in zip(rng_T_ext, outputs):
            k = str(int(round(T_ext.to('degC').m)))
            shelf[k] = output
