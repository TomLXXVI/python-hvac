"""
CREATE THE INDOOR UNITS OF THE VRF-SYSTEM
"""
from hvac import Quantity
from hvac.vrf_system.indoor_unit import IndoorUnit

Q_ = Quantity

# See document `capacity_tables.pdf` p. 12 (2-329) - Indoor unit temperature
# correction. Read points on the capacity ratio vs. indoor temperature curve.
points = [
    (15.0, 1.02),
    (16.0, 1.016),
    (17.0, 1.012),
    (18.0, 1.008),
    (19.0, 1.004),
    (20.0, 1.000),
    (21.0, 0.959),
    (22.0, 0.917),
    (23.0, 0.876),
    (24.0, 0.834),
    (25.0, 0.793),
    (26.0, 0.751),
    (27.0, 0.710)
]


def create_indoor_unit_15(name: str) -> IndoorUnit:
    iu_15 = IndoorUnit(
        name=name,
        type_='PLFY-P15VFM',
        model_size=15,
        Q_rated=Q_(1.9, 'kW'),
        iu_temp_corr=points
    )
    return iu_15


def create_indoor_unit_25(name: str) -> IndoorUnit:
    iu_25 = IndoorUnit(
        name=name,
        type_='PLFY-P25VFM',
        model_size=25,
        Q_rated=Q_(3.2, 'kW'),
        iu_temp_corr=points
    )
    return iu_25


def create_indoor_unit_32(name: str) -> IndoorUnit:
    iu_32 = IndoorUnit(
        name=name,
        type_='PLFY-P32VFM',
        model_size=32,
        Q_rated=Q_(4.0, 'kW'),
        iu_temp_corr=points
    )
    return iu_32


def create_indoor_unit_40(name: str) -> IndoorUnit:
    iu_40 = IndoorUnit(
        name=name,
        type_='PLFY-P40VFM',
        model_size=40,
        Q_rated=Q_(5.0, 'kW'),
        iu_temp_corr=points
    )
    return iu_40


def create_indoor_unit_50(name: str) -> IndoorUnit:
    iu_50 = IndoorUnit(
        name=name,
        type_='PLFY-P50VFM',
        model_size=50,
        Q_rated=Q_(6.3, 'kW'),
        iu_temp_corr=points
    )
    return iu_50


def get_indoor_units() -> list[tuple[IndoorUnit, Quantity]]:
    """Creates and returns the indoor units of the VRF-systems.

    Notes
    -----
    Next to each indoor unit, the space indoor air setpoint temperature is
    provided to determine the load-weighted average indoor air setpoint
    temperature of the entire building.
    """
    iu_list = [
        (create_indoor_unit_25('room_A'), Q_(23, 'degC')),
        (create_indoor_unit_25('room_B'), Q_(23, 'degC')),
        (create_indoor_unit_25('room_C'), Q_(23, 'degC')),
        (create_indoor_unit_50('room_D'), Q_(23, 'degC')),
        (create_indoor_unit_15('room_E'), Q_(23, 'degC')),
        (create_indoor_unit_40('room_F_1'), Q_(23, 'degC')),
        (create_indoor_unit_40('room_F_2'), Q_(23, 'degC')),
        (create_indoor_unit_25('room_G'), Q_(23, 'degC')),
        (create_indoor_unit_32('room H'), Q_(23, 'degC')),
        (create_indoor_unit_25('room_I'), Q_(23, 'degC')),
        (create_indoor_unit_50('room_J_1'), Q_(23, 'degC')),
        (create_indoor_unit_50('room_J_2'), Q_(23, 'degC')),
        (create_indoor_unit_50('room_J_3'), Q_(23, 'degC'))
    ]
    return iu_list
