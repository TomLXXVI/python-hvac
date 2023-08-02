"""Find the steady-state compressor speed of a single-stage vapor compression
machine at which the given evaporation and condensation will be reached under
given operating conditions.
"""
import warnings
from pathlib import Path
from hvac import Quantity
from hvac.logging import ModuleLogger
from hvac.fluids import Fluid, HumidAir, CoolPropWarning
from hvac.vapor_compression import VariableSpeedCompressor
from hvac.heat_transfer.heat_exchanger.fin_tube import air_evaporator, air_condenser
from hvac.vapor_compression.machine_bis import SingleStageVaporCompressionMachine

warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=CoolPropWarning)

logger = ModuleLogger.get_logger(__name__)

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
        n_cmp_max=Q_(4000, '1 / min'),
        logger=logger
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
    T_evp, T_cnd = Q_(6, 'degC'), Q_(51, 'degC')
    n_cmp = balance(
        machine,
        evp_m_dot_air=Q_(1500, 'kg / hr'),
        cnd_m_dot_air=Q_(3000, 'kg / hr'),
        evp_air_in=HumidAir(Tdb=Q_(24, 'degC'), RH=Q_(50, 'pct')),
        cnd_air_in=HumidAir(Tdb=Q_(35, 'degC'), RH=Q_(30, 'pct')),
        T_evp=T_evp,
        T_cnd=T_cnd
    )
    print(
        f"The evaporation temperature of {T_evp:~P} and the condensation "
        f"temperature of {T_cnd:~P} are reached at a compressor "
        f"speed around {n_cmp.to('1 / min'):~P.0f}."
    )
