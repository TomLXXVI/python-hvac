<!-- TOC -->
* [Heating Load Calculation of a House](#heating-load-calculation-of-a-house)
  * [Construction assembly](#construction-assembly)
    * [Building Component](#building-component)
      * [Building Materials Shelf](#building-materials-shelf)
    * [Airspace](#airspace)
    * [Surface Film](#surface-film)
    * [Creation of a Construction Assembly](#creation-of-a-construction-assembly)
      * [Construction Assemblies Shelf](#construction-assemblies-shelf)
  * [Window Thermal Properties](#window-thermal-properties)
    * [Window Properties Shelf](#window-properties-shelf)
  * [Building Elements](#building-elements)
<!-- TOC -->

# Heating Load Calculation of a House

The script `house.py` demonstrates how package `heating_load_calc` can be
used to program the heating load calculation of a house. Opening `floor_plan.pdf`
will show the floor plan of this example house.

## Construction assembly
The first step in a heating load calculation is to determine the construction 
assemblies that constitute the building elements of the house.

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

A building component is represented by the class **`BuildingComponent`**. Its
building material is represented by the class **`Material`**, and its geometry by
the class **`Geometry`**. These classes and their docstrings can be found in the 
module `construction_assembly.py` of `hvac.cooling_load_calc.core` (note that 
the heating load calculation package `hvac.heating_load_calc` uses the same 
`ConstructionAssembly` class and its underlying elements as the cooling load 
calculation package `hvac.cooling_load_calc`).
When a `BuildingComponent` object is instantiated, its thermal resistance `R` 
(SI-units: m².K/W) and thermal capacity `C` (SI-units: J/(m².K)) will be 
determined from the input data.

**Example 1:**
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

**Example 2:**
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
ISO 6946 (2017). An airspace or air layer is represented by the class **`AirSpace`**,
which can also be found in the module `construction_assembly.py` of sub-package 
`hvac.cooling_load_calc.core`.

**Example:**
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
wooden frame, however without insulation between the beams; only air is present. 
It should especially be noticed here that the geometry of the airspace also 
contains a parameter `w`, which refers to the width of the airspace (i.e. the 
smallest dimension of the air cavity between two wooden beams). If one of the 
planar dimensions of an airspace is smaller than 10 times its thickness, the 
thermal resistance of the airspace (air void) needs to be calculated in a 
different way according to standard ISO 6946, annex D.4. So, in case of an air 
void, the smallest planar dimension (the width) of the airspace must also be 
specified in the geometry for the thermal resistance to be calculated correctly. 
On the other hand, if this condition does not apply, parameter `w` can be 
omitted.

### Surface Film
The thermal resistance at the exterior and interior surface of a construction
assembly is also calculated according to standard ISO 6946 (2017). The surface 
film is represented by the class **`SurfaceLayer`**, which can also be found in the
module `construction_assembly.py` of sub-package `hvac.cooling_load_calc.core`.

**Example:**
```python
from hvac import Quantity
from hvac.heating_load_calc import SurfaceLayer, Geometry, HeatFlowDirection

Q_ = Quantity

ext_surf_film = SurfaceLayer.create(
    ID='ext_surf_film',
    geometry=Geometry(),
    heat_flow_direction=HeatFlowDirection.HORIZONTAL,
    Tmn=Q_(10, 'degC'),
    internal_surface=False,
    wind_speed=Q_(4, 'm / s')
)
```
In this example an exterior surface film has been created. To indicate that the
`SurfaceLayer` object represents the surface film at the exterior side of a 
construction assembly, the parameter `internal_surface` must be set to `False`
(the default value of this parameter is `True`). The default value for parameter
`wind_speed` is 4 m/s in accordance with ISO 6946 (2017), § 6.8, so we could 
have omitted this parameter in the example.

### Creation of a Construction Assembly
Objects of classes `BuildingComponent`, `AirSpace`, and `SurfaceLayer` form the
basis of a `ConstructionAssembly` object. The class **`ConstructionAssembly`** is 
located in the same module, `construction_assembly.py` of 
`hvac.cooling_load_calc.core`.

The example below demonstrates how a construction assembly is created. The 
example is taken from the module `exterior_walls.py` in sub-package 
`hvac.cooling_load_calc.wtcb_catalog`, that was already mentioned above and where 
more examples can be found by examining the modules in that package.

**Example 1:**
```python
from hvac import Quantity
from hvac.heating_load_calc import (
    ConstructionAssembly,
    SurfaceLayer,
    Geometry,
    HeatFlowDirection,
    BuildingComponent,
    MaterialsShelf,
    AirSpace,
    MechanicalFastening
)

Q_ = Quantity

def create_ca_ext_wall_wtcbF1(
    t_ins: Quantity,
    T_ext: Quantity = Q_(0.0, 'degC'),
    T_int: Quantity = Q_(20.0, 'degC'),
    T_asp: Quantity = Q_(10, 'degC'),
    dT_asp: Quantity = Q_(5, 'K'),
    v_wind: Quantity = Q_(4, 'm / s')
) -> ConstructionAssembly:
    """Creates and returns the construction assembly for an exterior wall 
    according to card 1 in the WTCB-catalog.
    
    Parameters
    ----------
    t_ins:
        Insulation thickness.
    T_ext:
        Outdoor temperature.
    T_int:
        Indoor temperature.
    T_asp:
        Mean thermodynamic temperature of the airspace.
    dT_asp:
        Temperature difference across the airspace.
    v_wind:
        Wind speed.
    """
    ext_surf_film = SurfaceLayer.create(
        ID='ext_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_ext,
        internal_surface=False,
        wind_speed=v_wind
    )
    outer_leaf = BuildingComponent.create(
        ID='outer_leaf',
        geometry=Geometry(t=Q_(10.0, 'cm')),
        material=MaterialsShelf.load('bakstenen gebakken aarde, 1500 kg/m3'),
    )
    air_space = AirSpace.create(
        ID='air_space',
        geometry=Geometry(t=Q_(5.0, 'cm')),
        dT=dT_asp,
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp,
        surface_emissivities=(Q_(0.9, 'frac'), Q_(0.9, 'frac')),
        inclination_angle=Q_(0.0, 'deg')
    )
    insulation = BuildingComponent.create(
        ID='insulation',
        geometry=Geometry(t=t_ins),
        material=MaterialsShelf.load('minerale wol, plaat + glasvlies, 22 kg/m3')
    )
    inner_leaf = BuildingComponent.create(
        ID='inner_leaf',
        geometry=Geometry(t=Q_(15.0, 'cm')),
        material=MaterialsShelf.load('blokken cellenbeton, gelijmd, 600 kg/m3')
    )
    gypsum_layer = BuildingComponent.create(
        ID='gypsum_layer',
        geometry=Geometry(t=Q_(1.5, 'cm')),
        material=MaterialsShelf.load('gipspleister')
    )
    int_surf_film = SurfaceLayer.create(
        ID='int_surf_film',
        geometry=Geometry(),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_int
    )
    ext_wall = ConstructionAssembly.create(
        ID=f'ext_wall_wtcbF1 (t_ins={t_ins.to("cm"):~P.0f})',
        layers=[
            ext_surf_film,
            outer_leaf,
            air_space,
            insulation,
            inner_leaf,
            gypsum_layer,
            int_surf_film
        ]
    )
    ext_wall = ext_wall.apply_ventilated_layer_correction(
        ventilated_layer_ID=air_space.ID,
        area_of_openings=Q_(1000, 'mm ** 2'),
        heat_flow_direction=HeatFlowDirection.HORIZONTAL,
        Tmn=T_asp
    )
    ext_wall = ext_wall.apply_insulation_correction(
        insulation_layer_ID=insulation.ID,
        insulation_level=2,
        mechanical_fastening=MechanicalFastening.create(
            diameter=Q_(2, 'mm'),
            number_per_unit_area=Q_(4, '1 / m ** 2'),
            insulation_thickness=t_ins,
            length=t_ins
        )
    )
    return ext_wall
```
The example above also demonstrates that it is possible to correct the thermal 
resistance of a construction under certain circumstances, namely when it contains
a slightly or well ventilated airspace, and/or when a correction is needed for the
thermal resistance of the insulation layer, to take the quality of installation 
into account or due to mechanical fasteners that penetrate the insulation layer. 
In the latter case, the mechanical fastening is specified through an object of 
class `MechanicalFastening`, defined in module `construction_assembly.py` of
`hvac.cooling_load_calc.core`. Both corrections are calculated according to 
standard ISO 6946 (2017).

It may happen that the thermal resistance *R* or conductance *U* of a 
construction assembly is already given, but without knowing the layered 
composition of it. In such case, the construction assembly can be created 
simply as shown in the next example, which applies to the construction assembly 
of an interior door.

**Example 2:**
```python
from hvac import Quantity
from hvac.heating_load_calc import (
    ConstructionAssembly,
    Geometry
)

Q_ = Quantity

interior_door = ConstructionAssembly.create(
    ID='interior_door',
    U=Q_(4.0, 'W / (m ** 2 * K)'),
    geometry=Geometry(t=Q_(4, 'cm'))
)
```

#### Construction Assemblies Shelf
For the same reason as for building materials, `ConstructionAssembly` objects can,
once they are created, be "shelved" for re-using them afterwards, using the class 
`ConstructionAssembliesShelf` in module 
`hvac.cooling_load_calc.shelves.construction_assemblies.py`. Usage examples of
class `ConstructionAssembliesShelf` can be found in the modules of sub-package 
`hvac.cooling_load_calc.wtcb_catalog`.

## Window Thermal Properties
The thermal properties of windows are gathered in objects instantiated from class
**`WindowThermalProperties`**, which resides in the module `fenestration.py` of
`hvac.cooling_load_calc.core`. These properties are the U-value of the entire 
window, including its frame and edge effects, and the solar heat gain coefficient 
(SHGC). In heating load calculations only the U-value of the window is of 
importance. The SHGC is used in cooling load calculations only to calculate the 
solar heat gain trough windows. In a heating load calculation on the other hand
it is assumed that no solar radiation is present. In fact, there are three values 
of SHGC to be distinguished when performing a cooling load calculation:
- the SHGC at the center-of-glazing related to direct solar radiation at 
different incidence angles (0°, 40°, 50°, 60°, 70°, and 80°)
- the SHGC at the center-of-glazing related to diffuse solar radiation
- the global SHGC of the entire window at normal incidence

Values for different types of glazing and window frames can be found in *ASHRAE
Fundamentals 2017*, Chapter 15, tables 4 and 10.

**Example**
```python
from hvac import Quantity
from hvac.heating_load_calc import WindowThermalProperties

Q_ = Quantity

wnd_props = WindowThermalProperties(
    ID='window_5a_operable_wood/vinyl',
    U=Q_(2.86, 'W / (m ** 2 * K)'),
    SHGC_cog_dir={
        0.00: 0.76,
        40.0: 0.74,
        50.0: 0.71,
        60.0: 0.64,
        70.0: 0.50,
        80.0: 0.26
    },
    SHGC_cog_dif=0.66,
    SHGC_wnd=0.62
)
```
In this example the thermal properties of "window 5-a" (identified the same way
as in ASHRAE Fundamentals 2017, table 4 and 10) with an operable frame made from 
wood or vinyl are gathered in a `WindowThermalProperties` object.

### Window Properties Shelf
For the same reason as for building materials and building materials,
`WindowThermalProperties` objects can, once they are created, be "shelved" for 
re-using them afterwards, using the class `WindowPropertiesShelf` in module 
`hvac.cooling_load_calc.shelves.window_properties.py`. A usage example can be 
found in the module `window_properties.py` of sub-package 
`hvac.cooling_load_calc.wtcb_catalog`.


## Building Elements
Once the construction assemblies are created or loaded from their respective 
shelf, we can create the building elements that constitute the building spaces.
In a heating load calculation there are three types of building elements to
be distinguished:
- exterior building elements separate the indoor environment of a heated space 
from the outdoor environment,
- adjacent building elements separate a heated space from another adjacent 
heated or unheated space, an adjacent building entity in the same building, or an
adjacent building.
- building elements which are in contact with the ground

Each type of building element is represented by its own class:
- `ExteriorBuildingElement`
- `AdjacentBuildingElement`
- `GroundBuildingElement`

These classes are defined in the module `building_element.py` of sub-package
`hvac.heating_load_calc.core`. 

Although it is possible, we don't need to instantiate these classes directly.
The building elements will be created when we add them to a space as explained
in the next section.
