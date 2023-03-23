# Heating Load Calculation of a House

The script `house.py` demonstrates how package `heating_load_calc` can be
used to program the heating load calculation of a house. Opening `floor_plan.pdf`
will show the floor plan of this example house.

## Construction assembly
The first step is to determine the construction assemblies that constitute the
building elements of the house. 

A construction assembly is a layered assembly of flat "thermal components". 
A "thermal component" represents any physical object that has thermal resistance 
and through which heat can flow from the building's interior environment to the 
outdoor environment or vice versa. There can be three different types of 
"thermal components" in a construction assembly:
- a **building component** is an opaque layer of building material 
- an **airspace** or air layer inside a construction assembly
- a **surface film** layer at both the exterior and the interior side of a 
construction assembly


### Building Component
A building component is made of a **building material**. A building material is
thermally characterized by three material properties:
- thermal conductivity *k*
- mass density *rho*
- specific heat capacity *c*

However, for some building materials it may not be possible to specify the 
thermal conductivity unambiguously, e.g. hollow concrete blocks. In such a
case the manufacturer specifies the unit thermal resistance *R* and the
thickness *t* of the building material.

Besides the building material, a building component is further characterized by
its **geometry**. Thickness *t* has already been mentioned. Other possible 
geometrical properties are: width, height, and area.

If width and height are specified, the area can be calculated. Otherwise, the
area of the building material needs to be specified. However, as we are mostly 
interested in the heat flux, i.e. the heat transfer per unit of area, through a 
construction assembly, the area of the building material will be 1 unit of area 
(e.g. 1 m²).

An exception to this general rule applies to cases where different building
materials are combined in a single layer of the construction assembly, e.g.
a wooden frame filled up with insulation material. When considering a relevant
area of this layer, we could determine the fraction taken up by the wood and the 
fraction taken up by the insulation material, e.g. 11 % wood and 89 % insulation.
In that case, if we consider an area of the layer of 1 m², 0.11 m² is wood and 
0.89 m² is insulation.

A building component is represented by the class `BuildingComponent`. Its
building material is represented by the class `Material`, and its geometry by
the class `Geometry`. These classes and their docstrings can be found in the 
module `construction_assembly.py` of `hvac.cooling_load_calc.core` (note that 
the heating load calculation package `hvac.heating_load_calc` uses the same 
`ConstructionAssembly` class and its underlying elements as the cooling load 
calculation package `hvac.cooling_load_calc`).
When a `BuildingComponent` object is instantiated, its thermal resistance `R` 
(SI-units: m².K/W) and thermal capacity `C` (SI-units: J/(m².K)) will be 
determined from the input data.

Example 1:
```python
from hvac import Quantity
from hvac.heating_load_calc import Geometry, Material, BuildingComponent

Q_ = Quantity

gypsum_layer = BuildingComponent.create(
    ID='gypsum_layer',
    geometry=Geometry(t=Q_(1.5, 'cm')),
    material=Material(
        k=Q_(0.56, 'W / (m * K)'),
        rho=Q_(1300, 'kg / m ** 3'),
        c=Q_(840, 'J / (kg * K)')
    )
)
```
Note that only the thickness *t* of the geometry has been specified. The area
*A* has a default value of 1 m².

Example 2:
```python
from hvac import Quantity
from hvac.heating_load_calc import Geometry, Material, BuildingComponent

Q_ = Quantity

insulation_ = BuildingComponent.create(
    ID='insulation_',
    geometry=Geometry(t=Q_(12, 'cm'), A=Q_(0.88, 'm ** 2')),
    material=Material(
        k=Q_(0.04, 'W / (m * K)'),
        rho=Q_(135.0, 'kg / m ** 3'),
        c=Q_(840.0, 'J / (kg * K)')
    )
)

wood_ = BuildingComponent.create(
    ID='wood_',
    geometry=Geometry(t=Q_(12, 'cm'), A=Q_(0.12, 'm ** 2')),
    material=Material(
        k=Q_(0.16, 'W / (m * K)'), 
        rho=Q_(550.0, 'kg / m ** 3'),
        c=Q_(1630.0, 'J / (kg * K)')
    )
)

insulation = insulation_ // wood_
insulation.ID = 'insulation'
```
In this second example the usage of parameter `A` in the geometry is 
demonstrated. In the example, a building component composed of a wooden
frame filled up with insulation is being considered. The area fraction of wood 
is 12 % and therefore the area fraction of insulation is 88 %. Consequently, we 
assign the parameter `A` of wood to be equal to 0.12 m² and the area of 
insulation to be equal to 0.88 m².

Class `BuildingComponent` inherits from class `ThermalComponent`. Objects
derived from this class and its child classes `BuildingComponent`, `AirSpace`,
and `SurfaceLayer` can be added in series (using the +-operator) or in parallel
(using the //-operator). In the example above, the wooden frame and the 
insulation form two thermal resistances connected in parallel. As such, the final
insulation layer `insulation` is created from a parallel combination of the 
wooden frame and the insulation material.

#### Building Materials Shelf
As it obviously may become a burden to re-enter the same building material over 
and over again, some kind of database is provided in which building materials 
can be stored. This database is called a "shelf" here, because it is implemented 
using the `shelve` module from Python's Standard Library. You can use the class 
`MaterialsShelf` in the module `materials` of `hvac.cooling_load_calc.shelves` 
to store building materials at a disk location of your choice. An example of the
usage of the class `MaterialsShelf` and its methods can be found in the 
sub-package `hvac.cooling_load_calc.wtcb_catalog`. (WTCB, nowadays Buildwise, is
a Belgian scientific and technical research centre for the construction sector.
They've published a catalog of construction assemblies frequently used in the 
Belgian construction sector, that I have implemented in the sub-package 
`wtcb_catalog`. If you run the script `materials.py` in this package for the 
first time, a file `wtcb_materials.db` (actually there will be 3 files) is 
created in a directory `wtcb_database` that will be located in the user's home 
directory).


### Airspace
The thermal resistance of an airspace is calculated according to standard 
ISO 6946 (2017). An airspace or air layer is represented by the class `AirSpace`,
which can also be found in the module `construction_assembly.py` of 
`hvac.cooling_load_calc.core`.

Example:
```python
from hvac import Quantity
from hvac.heating_load_calc import (
    AirSpace, 
    Geometry, 
    HeatFlowDirection,
    BuildingComponent,
    Material
)

Q_ = Quantity

airspace_ = AirSpace.create(
    ID='airspace_',
    geometry=Geometry(t=Q_(18, 'cm'), w=Q_(50, 'cm'), A=Q_(0.89, 'm ** 2')),
    dT=Q_(5, 'K'),
    heat_flow_direction=HeatFlowDirection.UPWARDS,
    Tmn=Q_(10, 'degC')
)

wood_ = BuildingComponent.create(
    ID='wood_',
    geometry=Geometry(t=Q_(18, 'cm'), A=Q_(0.11, 'm ** 2')),
    material=Material(
        k=Q_(0.16, 'W / (m * K)'), 
        rho=Q_(550.0, 'kg / m ** 3'),
        c=Q_(1630.0, 'J / (kg * K)')
    )
)
airspace = airspace_ // wood_
airspace.ID = 'airspace'
```
The example above applies again to a construction assembly layer composed of a 
wooden frame, however without insulation between the beams, but only air. It
should be especially noticed here that the geometry of the airspace also 
contains a parameter `w`, which refers to the width of the airspace (i.e. the 
smallest dimension of the air pocket between two wooden beams). If one of the 
planar dimensions of an airspace is smaller than 10 times its thickness, the 
thermal resistance of the airspace (air pocket or air void) needs to be 
calculated in a different way according to standard ISO 6946, annex D.4. So, in
case of an air void or pocket, the smallest planar dimension (the width) of the
airspace must also be specified in the geometry, so that the thermal resistance
can be calculated correctly. On the other hand, if this condition does not 
apply, parameter `w` can be omitted.
