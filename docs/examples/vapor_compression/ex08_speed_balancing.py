import warnings
from pathlib import Path
from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import Fluid, HumidAir
from hvac.vapor_compression import VariableSpeedCompressor
from hvac.heat_transfer.heat_exchanger.fin_tube import air_evaporator, air_condenser
from hvac.vapor_compression.machine_bis import SingleStageVaporCompressionMachine

warnings.filterwarnings('ignore', category=RuntimeWarning)

logger = ModuleLogger.get_logger('hvac.fluids.fluid')
logger.setLevel(ModuleLogger.CRITICAL)

Q_ = Quantity

R134a = Fluid('R134a')
Evaporator = air_evaporator.rating.PlainFinTubeCounterFlowEvaporator
Condenser = air_condenser.rating.PlainFinTubeCounterFlowCondenser


def create_machine() -> SingleStageVaporCompressionMachine:
    compressor = VariableSpeedCompressor(
        coeff_file=Path("./compressor_data/VTZ054-G_R134a.csv"),
        refrigerant_type=R134a,
        units={'m_dot': 'kg / hr', 'speed': '1 / s'}
    )
    evaporator = Evaporator(
        L1=Q_(0.731, 'm'),
        L3=Q_(0.244, 'm'),
        N_r=3,
        S_t=Q_(25.4, 'mm'),
        S_l=Q_(22.0, 'mm'),
        D_i=Q_(8.422, 'mm'),
        D_o=Q_(10.2, 'mm'),
        t_f=Q_(0.3302, 'mm'),
        N_f=1 / Q_(3.175, 'mm')
    )
    condenser = Condenser(
        L1=Q_(0.999, 'm'),
        L3=Q_(0.333, 'm'),
        N_r=5,
        S_t=Q_(25.4, 'mm'),
        S_l=Q_(22.0, 'mm'),
        D_i=Q_(8.422, 'mm'),
        D_o=Q_(10.2, 'mm'),
        t_f=Q_(0.3302, 'mm'),
        N_f=1 / Q_(3.175, 'mm')
    )
    machine = SingleStageVaporCompressionMachine(
        evaporator, condenser, compressor, R134a,
        n_cmp_min=Q_(2100, '1 / min'),
        n_cmp_max=Q_(4000, '1 / min')
    )
    return machine


def balance(
    machine: SingleStageVaporCompressionMachine,
    evp_m_dot_air: Quantity,
    cnd_m_dot_air: Quantity,
    evp_air_in: HumidAir,
    cnd_air_in: HumidAir,
    T_evp: Quantity,
    T_cnd: Quantity
) -> Quantity:
    machine.set_operating_conditions(
        evp_m_dot_air=evp_m_dot_air,
        cnd_m_dot_air=cnd_m_dot_air,
        evp_air_in=evp_air_in,
        cnd_air_in=cnd_air_in,
        dT_sh=Q_(10, 'K')
    )
    n_cmp = machine.balance_by_speed(T_evp, T_cnd)
    return n_cmp


if __name__ == '__main__':
    machine = create_machine()
    n_cmp = balance(
        machine,
        evp_m_dot_air=Q_(1500, 'kg / hr'),
        cnd_m_dot_air=Q_(3000, 'kg / hr'),
        evp_air_in=HumidAir(Tdb=Q_(24, 'degC'), RH=Q_(50, 'pct')),
        cnd_air_in=HumidAir(Tdb=Q_(35, 'degC'), RH=Q_(30, 'pct')),
        T_evp=Q_(6, 'degC'),
        T_cnd=Q_(51, 'degC')
    )
    print(n_cmp)
