"""
EXAMPLE 9
---------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.12.1, p. 82:
Estimate the fraction of the average June radiation on a horizontal surface
that is diffuse in Madison, WI.
"""
from hvac import Quantity
from hvac.sun.surface import Location
from hvac.sun.radiation import estimate_avg_H_d_fraction
from hvac.sun.radiation import ReferenceDates

Q_ = Quantity

# From Appendix D, the June average daily radiation for Madison is 23 MJ/m²:
H_avg = Q_(23, 'MJ / m**2')

location = Location(
    fi=Q_(43, 'deg'),
    date=ReferenceDates.get_date_for('June')
)

K_T_avg = H_avg / location.sun.H_o

H_d_avg_frac = estimate_avg_H_d_fraction(
    K_T_avg=K_T_avg.to('frac').m,
    omega_ss=location.omega_ss.to('rad').m
)
print(f"Fraction of diffuse radiation: {H_d_avg_frac:.2f}")
