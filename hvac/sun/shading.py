import numpy as np
from hvac import Quantity
from hvac.sun.surface import Surface


def shading_point_coords(
    surface: Surface,
    z: Quantity
) -> tuple[Quantity, Quantity]:
    """Calculates the horizontal and vertical distance x and y of a shading
    point S located on the given surface with respect to its obstacle point P
    that is situated at the given distance z in front of the surface.

    Parameters
    ----------
    surface:
        `Surface` object encapsulating the location of the surface, the
        orientation of the surface, and the position of the sun.
    z:
        The distance between the obstacle point P in front of the surface and
        this surface.

    Returns
    -------
    The horizontal and vertical distance of the shading point S on the surface,
    measured from the obstacle point P.
    """
    x = z * np.tan(surface.location.sun.gamma - surface.gamma)
    y = z * np.tan(surface.alpha_p)
    return x, y


if __name__ == '__main__':

    from datetime import date, time
    from hvac.sun.surface import Location

    Q_ = Quantity

    location = Location(fi=Q_(45, 'deg'))
    location.date = date(2022, 7, 21)
    location.solar_time = time(13, 0, 0)
    window = Surface(location, gamma=Q_(0.0, 'deg'), beta=Q_(90.0, 'deg'))

    print(
        f"solar azimuth angle = {location.sun.gamma.to('deg'):~P.1f}\n"
        f"solar altitude angle = {location.sun.alpha.to('deg'):~P.1f}\n"
        f"window azimuth angle = {window.gamma.to('deg'):~P.1f}\n"
        f"window profile angle = {window.alpha_p.to('deg'):~P.1f}\n"
    )

    z = Q_(0.5, 'm')
    x, y = shading_point_coords(window, z)

    print(
        f"horizontal distance x of shading point = {x.to('m'):~P.2f}\n"
        f"vertical distance y of shading point = {y.to('m'):~P.2f}"
    )
