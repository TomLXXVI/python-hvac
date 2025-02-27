"""
EXAMPLE 2
---------
EXTERIOR SURFACE - SOL-AIR TEMPERATURE.
"""
from hvac import Quantity
from hvac.sun import Location
from hvac.sun.radiation import ClimateType, ReferenceDates
from hvac.cooling_load_calc_old.core.weather_data import WeatherData
from hvac.cooling_load_calc_old.core.building_element import ExteriorSurface

Q_ = Quantity

# ------------------------------------------------------------------------------
# LOCATION
location = Location(
    fi=Q_(51.183, 'deg'),
    L_loc=Q_(3.8, 'deg'),
    altitude=Q_(8, 'm'),
    climate_type=ClimateType.MID_LATITUDE_SUMMER,
    timezone='Etc/GMT-1'
)

# ------------------------------------------------------------------------------
# EXTERIOR SURFACE WITH CLEAR-SKY MODEL:
weather_data_01 = WeatherData.create_from_climatic_design_data(
    location=location,
    date=ReferenceDates.get_date_for('Jun'),
    T_db_des=Q_(24.2, 'degC'),
    T_db_rng=Q_(11.8, 'K'),
    T_wb_mc=Q_(18.0, 'degC'),
    T_wb_rng=Q_(5.6, 'K')
)
ext_surf_01 = ExteriorSurface(
    weather_data=weather_data_01,
    gamma=Q_(0, 'deg'),
    beta=Q_(60, 'deg')
)
print(
    "clear-sky model: ",
    f"dry-bulb temperature = {weather_data_01.T_db(12):~P.1f}",
    f"irradiation on tilted surface = {weather_data_01.I_T(12, ext_surf_01)[0].to('kJ / m**2'):~P.1f}",
    f"sol-air temperature = {ext_surf_01.T_sa(12 * 3600):~P.1f}",
    sep='\n'
)

# ------------------------------------------------------------------------------
# EXTERIOR SURFACE WITH A RANDOM DAILY SOLAR RADIATION DISTRIBUTION DERIVED FROM
# THE MONTHLY AVERAGE DAILY IRRADIATION
weather_data_02 = WeatherData.create_from_climatic_design_data(
    location=location,
    date=ReferenceDates.get_date_for('Jun'),
    T_db_des=Q_(24.2, 'degC'),
    T_db_rng=Q_(11.8, 'K'),
    T_wb_mc=Q_(18.0, 'degC'),
    T_wb_rng=Q_(5.6, 'K'),
    H_avg=Q_(5.38, 'kWh / m**2')
)
ext_surf_02 = ExteriorSurface(
    weather_data=weather_data_02,
    gamma=Q_(0, 'deg'),
    beta=Q_(60, 'deg')
)
print(
    "correlation: ",
    f"dry-bulb temperature = {weather_data_02.T_db(12):~P.1f}",
    f"irradiation on tilted surface = {weather_data_02.I_T(12, ext_surf_02)[0].to('kJ / m**2'):~P.1f}",
    f"sol-air temperature = {ext_surf_02.T_sa(12 * 3600):~P.1f}",
    sep='\n'
)

# ------------------------------------------------------------------------------
# EXTERIOR SURFACE WITH TMY DATA:
weather_data_03 = WeatherData.create_from_tmy_data(
    location=location,
    date=ReferenceDates.get_date_for('Jun'),
    tmy_file='tmy_50.911_3.192_2005_2020.csv'
)
ext_surf_03 = ExteriorSurface(
    weather_data=weather_data_03,
    gamma=Q_(0, 'deg'),
    beta=Q_(60, 'deg')
)
print(
    "tmy: ",
    f"dry-bulb temperature = {weather_data_03.T_db(12):~P.1f}",
    f"irradiation on tilted surface = {weather_data_03.I_T(12, ext_surf_03)[0].to('kJ / m**2'):~P.1f}",
    f"sol-air temperature = {ext_surf_03.T_sa(12 * 3600):~P.1f}",
    sep='\n'
)
