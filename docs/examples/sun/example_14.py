"""
EXAMPLE 14
----------
From: Duffie, J. A., Beckman, W. A., & Blair, N. (2020).
SOLAR ENGINEERING OF THERMAL PROCESSES, PHOTOVOLTAICS AND WIND.
John Wiley & Sons.

Example 2.15.1, p. 92:
Using the isotropic diffuse model, estimate the beam, diffuse, and
ground-reflected components of solar radiation and the total radiation on a
surface sloped 60° toward the south at a latitude of 40 °N for the hour 9 to
10 a.m. on February 20. Here I = 1.04 MJ/m² on the horizontal surface and
rho_g = 0.60.

Example 2.16.1, p. 94:
Do example 2.15.1 using the HDKR model.

Example 2.16.2, p. 97:
Do example 2.15.1 using the Perez model.
"""
from datetime import date, time
from hvac import Quantity
from hvac.sun.surface import Location, Surface
from hvac.sun.transposition import IsotropicSkyModel, AnisotropicSkyModel

Q_ = Quantity


def print_results(G_Tb, G_Td, G_Tg, G_T):
    print(
        "hourly beam irradiation on tilted surface = "
        f"{G_Tb.to('MJ/m**2'):~P.3f}",
        "hourly diffuse irradiation on tilted surface = "
        f"{G_Td.to('MJ/m**2'):~P.3f}",
        "hourly ground-reflected irradiation on tilted surface = "
        f"{G_Tg.to('MJ/m**2'):~P.3f}",
        "hourly total irradiation on tilted surface = "
        f"{G_T.to('MJ/m**2'):~P.3f}",
        sep='\n'
    )


def main():
    location = Location(
        fi=Q_(40, 'deg'),
        date=date(2023, 2, 20)
    )

    surface = Surface(
        location=location,
        gamma=Q_(0, 'deg'),
        beta=Q_(60, 'deg')
    )
    
    time_interval = (time(9), time(10))
    
    # First, we need to determine the beam and diffuse component of the hourly
    # irradiation on the horizontal surface:
    I_b, I_d = location.sun.split_radiation(
        I=Q_(1.04, 'MJ / m**2'),
        time_interval=time_interval
    )
    
    iso_sky_model = IsotropicSkyModel(
        surface=surface,
        rho_g=Q_(0.6, 'frac'),
        G_b=I_b,
        G_d=I_d,
        time_interval=time_interval
    )
    print('Isotropic sky model:')
    print_results(
        iso_sky_model.G_Tb,
        iso_sky_model.G_Td,
        iso_sky_model.G_Tg,
        iso_sky_model.G_T
    )
    
    HDKR_model = AnisotropicSkyModel.HDKR(
        surface=surface,
        rho_g=Q_(0.6, 'frac'),
        G_b=I_b,
        G_d=I_d,
        time_interval=time_interval
    )
    print('HDKR model')
    print_results(
        HDKR_model.G_Tb,
        HDKR_model.G_Td,
        HDKR_model.G_Tg,
        HDKR_model.G_T
    )

    Perez_model = AnisotropicSkyModel.Perez(
        surface=surface,
        rho_g=Q_(0.6, 'frac'),
        G_b=I_b,
        G_d=I_d,
        time_interval=time_interval
    )
    print('Perez model')
    print_results(
        Perez_model.G_Tb,
        Perez_model.G_Td,
        Perez_model.G_Tg,
        Perez_model.G_T
    )


if __name__ == '__main__':
    main()
