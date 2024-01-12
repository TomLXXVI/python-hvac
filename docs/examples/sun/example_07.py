"""
EXAMPLE 7
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.8.2, p. 72:
Estimate the standard clear-day radiation on a horizontal surface for Madison on
August 22.
"""
import typing
from datetime import date, time, datetime, timedelta
from hvac import Quantity
from hvac.sun.surface import Location
from hvac.sun.radiation import ClimateType

Q_ = Quantity

location = Location(
    fi=Q_(43, 'deg'),
    altitude=Q_(270, 'm'),
    climate_type=ClimateType.MID_LATITUDE_SUMMER,
    date=date(2023, 8, 22)
)
# altitude and climate type are needed to determine the atmospheric
# transmittance for clear-sky beam radiation

print(
    f"sunrise at {location.sunrise}",
    f"sunset at {location.sunset}",
    sep='\n', end='\n\n'
)

# METHOD 1:

sr_hr = location.sunrise.hour
ss_hr = location.sunset.hour
t_lst = [time(h, 30, 0) for h in range(sr_hr, ss_hr)]

# determine clear-sky beam radiation and clear-sky diffuse radiation for every
# hour between sunrise and sunset:
I_cb_lst, I_cd_lst = [], []
for t in t_lst:
    dt = datetime.combine(location.date, t)
    t1 = (dt - timedelta(hours=0.5)).time()
    t2 = (dt + timedelta(hours=0.5)).time()
    I_cb_lst.append(location.sun.I_cb(t1, t2))
    I_cd_lst.append(location.sun.I_cd(t1, t2))

# determine daily clear-sky radiation:
H_c = sum(I_cb + I_cd for I_cb, I_cd in zip(I_cb_lst, I_cd_lst))
H_c = typing.cast(Quantity, H_c)
print(f"daily radiation = {H_c.to('MJ/m**2'):~P.3f} (method 1)")

# METHOD 2:

H_cb = location.sun.I_cb(location.sunrise, location.sunset)
H_cd = location.sun.I_cd(location.sunrise, location.sunset)
H_c = H_cb + H_cd
print(f"daily radiation = {H_c.to('MJ/m**2'):~P.3f} (method 2)")

# METHOD 3:

I_cb_lst, I_cd_lst = [], []
for t in t_lst:
    location.solar_time = t
    I_cb_lst.append(location.sun.G_cb * Q_(1, 'hr'))
    I_cd_lst.append(location.sun.G_cd * Q_(1, 'hr'))

H_c = sum(I_cb + I_cd for I_cb, I_cd in zip(I_cb_lst, I_cd_lst))
H_c = typing.cast(Quantity, H_c)
print(f"daily radiation = {H_c.to('MJ/m**2'):~P.3f} (method 3)")
