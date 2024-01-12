"""
EXAMPLE 12
----------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.13.3, p. 86
The average daily June total radiation on a horizontal plane in Madison is 23
MJ/m². Estimate the average diffuse, the average beam, and the average total
radiation for the hours 10 to 11 and 1 to 2.
"""
from datetime import time
from hvac import Quantity
from hvac.sun.radiation import (
    ReferenceDates,
    estimate_r_t,
    estimate_r_d,
    estimate_avg_H_d_fraction
)
from hvac.sun.surface import Location

Q_ = Quantity

location = Location(
    fi=Q_(43, 'deg'),
    date=ReferenceDates.get_date_for('Jun'),
    solar_time=time(10, 30, 0)
)

H_daily_avg = Q_(23, 'MJ / m**2')

K_T_avg = H_daily_avg / location.sun.H_o

r_H_d_avg = estimate_avg_H_d_fraction(
    K_T_avg.to('frac').m,
    location.omega_ss.to('rad').m
)
H_daily_d_avg = r_H_d_avg * H_daily_avg

r_t = estimate_r_t(location.omega.to('rad').m, location.omega_ss.to('rad').m)
I_t = r_t * H_daily_avg

r_d = estimate_r_d(location.omega.to('rad').m, location.omega_ss.to('rad').m)
I_d = r_d * H_daily_d_avg

I_b = I_t - I_d

print(
    f"average beam radiation: {I_b.to('MJ / m**2'):~P.2f}",
    f"average diffuse radiation: {I_d.to('MJ / m**2'):~P.2f}",
    f"average total radiation: {I_t.to('MJ / m**2'):~P.2f}",
    sep='\n'
)
