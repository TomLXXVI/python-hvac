"""
EXAMPLE 1
---------
WEATHER DATA - TEMPERATURE PROFILES AND SOLAR RADIATION PROFILE.
"""
from hvac import Quantity
from hvac.sun import Location
from hvac.sun.radiation import ClimateType, ReferenceDates
from hvac.cooling_load_calc.core.weather_data import WeatherData


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
# WEATHER DATA FROM CLEAR-SKY MODEL:
weather_data_01 = WeatherData.create_from_climatic_design_data(
    location=location,
    date=ReferenceDates.get_date_for('Jun'),
    T_db_des=Q_(24.2, 'degC'),
    T_db_rng=Q_(11.8, 'K'),
    T_wb_mc=Q_(18.0, 'degC'),
    T_wb_rng=Q_(5.6, 'K')
)
weather_data_01.plot_daily_temperature_profile()
weather_data_01.plot_daily_solar_radiation_profile()

# ------------------------------------------------------------------------------
# WEATHER DATA FROM DAILY SOLAR RADIATION DISTRIBUTION DERIVED FROM MONTHLY
# AVERAGE:
weather_data_02 = WeatherData.create_from_climatic_design_data(
    location=location,
    date=ReferenceDates.get_date_for('Jun'),
    T_db_des=Q_(24.2, 'degC'),
    T_db_rng=Q_(11.8, 'K'),
    T_wb_mc=Q_(18.0, 'degC'),
    T_wb_rng=Q_(5.6, 'K'),
    H_avg=Q_(5.38, 'kWh / m**2')
)
weather_data_02.plot_daily_temperature_profile()
weather_data_02.plot_daily_solar_radiation_profile()

# ------------------------------------------------------------------------------
# WEATHER DATA FROM TMY DATA:
weather_data_03 = WeatherData.create_from_tmy_data(
    location=location,
    date=ReferenceDates.get_date_for('Jun'),
    tmy_file='tmy_50.911_3.192_2005_2020.csv'
)
weather_data_03.plot_daily_temperature_profile()
weather_data_03.plot_daily_solar_radiation_profile()
