"""
COPPER PIPING
-------------
Create copper piping shelves.
"""
from pathlib import Path
from hvac import Quantity
from hvac.refrigerant_piping import CopperTubing

Q_ = Quantity


# COPPER TYPE L
# nominal diameter (inch), outside diameter (inch), inside diameter (inch),
# pipe wall thickness (inch).
class CopperL:
    sizes = [
        ('1/4', 0.375, 0.315, 0.030),
        ('3/8', 0.500, 0.430, 0.035),
        ('1/2', 0.625, 0.545, 0.040),
        ('5/8', 0.750, 0.666, 0.042),
        ('3/4', 0.875, 0.785, 0.045),
        ('1', 1.125, 1.025, 0.050),
        ('1 1/4', 1.375, 1.265, 0.055),
        ('1 1/2', 1.625, 1.505, 0.060),
        ('2', 2.125, 1.985, 0.070),
        ('2 1/2', 2.625, 2.465, 0.080),
        ('3', 3.125, 2.945, 0.090),
        ('3 1/2', 3.625, 3.425, 0.100),
        ('4', 4.125, 3.905, 0.110),
        ('5', 5.125, 4.875, 0.125),
        ('6', 6.125, 5.845, 0.140),
        ('8', 8.125, 7.725, 0.200),
        ('10', 10.125, 9.625, 0.250),
        ('12', 12.125, 11.565, 0.280)
    ]


# COPPER TYPE ACR A
# nominal diameter (inch), outside diameter (inch), inside diameter (inch),
# pipe wall thickness (inch).
class CopperACRAnnealed:
    sizes = [
        ('1/8', 0.125, 0.065, 0.030),
        ('3/16', 0.187, 0.128, 0.030),
        ('1/4', 0.250, 0.190, 0.030),
        ('5/16', 0.312, 0.248, 0.032),
        ('3/8', 0.375, 0.311, 0.032),
        ('1/2', 0.500, 0.436, 0.032),
        ('5/8', 0.625, 0.555, 0.035),
        ('3/4', 0.750, 0.680, 0.035),
        ('7/8', 0.875, 0.785, 0.045),
        ('1 1/8', 1.125, 1.025, 0.050),
        ('1 3/8', 1.375, 1.265, 0.055),
        ('1 5/8', 1.625, 1.505, 0.060),
    ]


# COPPER TYPE ACR D
# nominal diameter (inch), outside diameter (inch), inside diameter (inch),
# pipe wall thickness (inch).
class CopperACRDrawnTemper:
    sizes = [
        ('1/4', 0.250, 0.200, 0.025),
        ('3/8', 0.375, 0.315, 0.030),
        ('1/2', 0.500, 0.430, 0.035),
        ('5/8', 0.625, 0.545, 0.040),
        ('3/4', 0.750, 0.666, 0.042),
        ('7/8', 0.875, 0.785, 0.045),
        ('1 1/8', 1.125, 1.025, 0.050),
        ('1 3/8', 1.375, 1.265, 0.055),
        ('1 5/8', 1.625, 1.505, 0.060),
        ('2 1/8', 2.125, 1.985, 0.070),
        ('2 5/8', 2.625, 2.465, 0.080),
        ('3 1/8', 3.125, 2.945, 0.090),
        ('3 5/8', 3.625, 3.425, 0.100),
        ('4 1/8', 4.125, 3.905, 0.110)
    ]


def create_shelf(sizes_: list[tuple[str, float, float, float]]) -> None:
    """Create a shelf with the copper tube dimensions in list `sizes_`."""
    for record in sizes_:
        DN = record[0]
        D_ext = Q_(record[1], 'inch')
        D_int = Q_(record[2], 'inch')
        t = Q_(record[3], 'inch')
        CopperTubing.add_record(DN, D_ext, D_int, t)


def main():
    # Create a shelf with copper tubes of type ACR A:
    CopperTubing.db_path = Path("./copper-tube-ACR-A")
    create_shelf(CopperACRAnnealed.sizes)
    
    # Create a shelf with copper tubes of type ACR D:
    CopperTubing.db_path = Path("./copper-tube-ACR-D")
    create_shelf(CopperACRDrawnTemper.sizes)
    
    # Create a shelf with copper tubes of type L:
    CopperTubing.db_path = Path("./copper-tube-L")
    create_shelf(CopperL.sizes)


if __name__ == '__main__':
    main()
