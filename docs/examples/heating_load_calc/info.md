TODO: the information in this file is not fully up to date anymore! 

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
  * [Building Hierarchy](#building-hierarchy)
  * [Programming a Heat Load Calculation](#programming-a-heat-load-calculation-)
    * [1. Imports](#1-imports)
    * [2. Set up the Shelves](#2-set-up-the-shelves)
    * [3. Create the Construction Assemblies](#3-create-the-construction-assemblies)
    * [4. Create the Building](#4-create-the-building)
    * [5. Run the Heating Load Calculation](#5-run-the-heating-load-calculation)
<!-- TOC -->

# Heating Load Calculation of a House

This document explains the build-up of package `hvac.heating_load_calc` and
how it should be used to program the heating load calculation of a building. 
This package implements the "standard method" (§ 6) of standard EN 12831-1 
(2017) *Energy performance of buildings - Method for calculation of the 
design heat load - Part 1: Space heating load.*
In the text that follows, we first go over the main concepts and classes that 
are needed to program the heat loss calculation of a building.

The script `house.py` in the folder `example_01` demonstrates 
how package `hvac.heating_load_calc` can be used to program the heating load 
calculation of a simple house. Opening `floor_plan.pdf` in the folder `example_01` 
will show the floor plan of this example house.

## Construction assembly
The first step in a heating load calculation is to determine the construction 
assemblies that constitute the building elements of the house.

A construction assembly is a layered assembly of flat "thermal components". 
A "thermal component" represents any physical object that has thermal resistance 
and through which heat can flow from the building's interior environment to the 
outdoor environment or vice versa. There can be three different types of 
"thermal components" in a construction assembly:

- a **building component** is an opaque layer made of some building material 
- an **airspace** or air layer can be present inside a construction assembly
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
construction assembly, the default area of the building material will be 1 unit 
of area (e.g. 1 m²).

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

**It is important to keep in mind that the layers of a construction assembly
are always ordered from outdoors to indoors!**

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
For the same reason as for building materials and construction assemblies,
`WindowThermalProperties` objects can, once they are created, be "shelved" for 
re-using them afterwards, using the class `WindowPropertiesShelf` in module 
`hvac.cooling_load_calc.shelves.window_properties.py`. A usage example can be 
found in the module `window_properties.py` of sub-package 
`hvac.cooling_load_calc.wtcb_catalog`.


## Building Elements
Once the construction assemblies are created or loaded from their respective 
shelf, we can create the building elements that will surround the building spaces.
In a heating load calculation there are three types of building elements to
be distinguished, each represented by its own class, located in the module 
`building_element.py` of sub-package `hvac.heating_load_calc.core`:
- Class **`ExteriorBuildingElement`** represents an exterior building element that
separates the indoor environment of a heated space from the outdoor environment.
- Class **`AdjacentBuildingElement`** represents an adjacent building element that
separates a heated space from another adjacent heated or unheated space, an 
adjacent building entity in the same building, or an adjacent building.
- Class **`GroundBuildingElement`** represents a building element in contact with 
the ground.

Although it is possible, we don't need to instantiate these classes directly.
The building elements will be created when we add them to a space as explained
in the following section.

## Building Hierarchy
The standard method §6 of EN 12831-1 divides a building into one or more building
entities, e.g. apartments in an apartment building. Each building entity can have
one or more ventilation zones. A ventilation zone is defined in the standard as: 
*group of rooms that are air-connected by design, either directly or indirectly 
(through other rooms there between); e.g. through internally mounted air transfer
devices / shortened door leafs, etc. By design, there is no air transfer between 
ventilation zones. Usually, each building entity is considered a separate zone.*
Each ventilation zone can have one or more heated or unheated spaces.

The table below shows which class represents each "physical object" and in 
which module each of these classes is implemented.

| physical object  | class             | module                                                |
|------------------|-------------------|-------------------------------------------------------|
| building         | `Building`        | `hvac.heating_load_calc.building.building.py`         |
| building entity  | `BuildingEntity`  | `hvac.heating_load_calc.building.building_entity.py`  |
| ventilation zone | `VentilationZone` | `hvac.heating_load_calc.building.ventilation_zone.py` |
| heated space     | `HeatedSpace`     | `hvac.heating_load_calc.building.space.py`            |
| unheated space   | `UnheatedSpace`   | `hvac.heating_load_calc.building.space.py`            |

It is not the purpose of this document to explain the calculation method in 
detail and to explain all the associated input parameters required to perform the 
calculations. More explanation about the meaning of the parameters can be found 
in the docstrings that document the classes and methods. It's also recommended
to have the standard EN 12831-1 (2017) at your own disposal. Many input parameters
have already been assigned reasonable default values that can be found in this 
standard.

## Programming a Heat Load Calculation 
In the following, an example will be used to explain how a program for the 
heating load calculation of a building can be structured and which classes and 
functions from the package are required or can be used for this. 
In the script `house.py`, in the directory `example_01`, the heating load
calculation of a simple two-storey house has been programmed. Also in this 
directory is the floor plan of this house (`floor_plan.pdf`). The bottom of the 
floor plan is oriented to the South, the right side of the floor plan points to 
the East. At the east, the north, and south side the house is open. At the west 
side a neighboring house is built against it. The toilet on the first floor and 
the hallway are considered to be unheated spaces. The house is mechanically 
ventilated, as well for the supply, as the exhaust of air. The supply, exhaust 
and transfer volume flow rates are indicated on the floor plan.

### 1. Imports
Of course the first step when writing a Python script is to import everything
we'll need:

```python
from hvac import Quantity
from hvac.heating_load_calc import (
    ClimateDesignData,
    Building,
    BuildingEntity,
    VentilationZone,
    HeatedSpace,
    UnheatedSpace,
    ConstructionAssembliesShelf,
    WindowPropertiesShelf,
    ConstructionAssembly,
    HeatFlowDirection,
    Geometry,
    BuildingComponent,
    SurfaceLayer,
    MaterialsShelf
)

import hvac.cooling_load_calc_old.wtcb_catalog.exterior_walls as ew
import hvac.cooling_load_calc_old.wtcb_catalog.interior_walls as iw
import hvac.cooling_load_calc_old.wtcb_catalog.floors as fl
import hvac.cooling_load_calc_old.wtcb_catalog.roofs as rf

Q_ = Quantity
```
I will use the sub-package `wtcb_catalog` in `hvac.cooling_load_calc` that contains
a number of modules in which I have programmed the construction assemblies from
the WTCB catalog, already mentioned above. Each construction assembly is 
implemented inside a function with some arguments that allow adapting the U-value
of a construction assembly depending on some specific parameters, e.g. to set
the insulation thickness of an exterior wall.

### 2. Set up the Shelves
I will use the shelves installed in my home directory to retrieve building 
materials, construction assemblies and window properties. To use them, the paths
to the appropriate db-files must be set first:
````python
ConstructionAssembliesShelf.path = "C:/Users/Tom/wtcb_database/wtcb_construction_assemblies.db"
WindowPropertiesShelf.path = "C:/Users/Tom/wtcb_database/window_properties.db"
MaterialsShelf.path = "C:/Users/Tom/wtcb_database/wtcb_materials.db"
````
If you ran the scripts in the `hvac.cooling_load_calc.shelves` package, the
'wtcb_database' directory will be installed in your own home directory, and you
would have to change the file paths in your program script accordingly.

### 3. Create the Construction Assemblies
Now, the construction assemblies of the building elements that surround
the spaces of the building can be created. This is done in a separate class.
(For larger projects, instead of class, you could also use a separate module for
this.) The construction assemblies are instance attributes created in 
the `__init__()`-method of the class. The script also demonstrates how to create 
project-specific construction assemblies or construction assemblies that are
not already available on the shelf, inside a static method of the class.
````python
class ConstructionAssemblies:
    # prepare the construction assemblies used for the building
    
    def __init__(self):
        # exterior walls
        self.ext_wall_WTCB1 = ew.create_ca_ext_wall_wtcbF1(
            t_ins=Q_(12, 'cm')
        )
        # interior walls
        self.int_wall_WTCB1 = iw.create_ca_int_wall_wtcbF1(
            t_ins=Q_(6, 'cm')
        )
        self.int_wall_PROJ1 = self._create_ca_int_wall_projF1()
        # floors
        self.floor_WTCB4 = fl.create_ca_floor_wtcbF4(
            t_ins=Q_(12, 'cm'),
            heat_flow_direction=HeatFlowDirection.DOWNWARDS
        )
        # ceilings
        self.ceiling_WTCB13 = rf.create_ca_ceiling_wtcbF13(
            t_ins=Q_(12, 'cm'),
            heat_flow_direction=HeatFlowDirection.UPWARDS,
            T_ext=Q_(10, 'degC')
        )
        # windows
        self.window_ASHRAE5a = WindowPropertiesShelf.load('window_5a_operable_wood/vinyl')
        # interior doors
        self.int_door = ConstructionAssembly.create(
            ID='door_int',
            U=Q_(4.0, 'W / (m ** 2 * K)'),
            geometry=Geometry(t=Q_(4, 'cm'))
        )
        # exterior doors
        self.ext_door = ConstructionAssembly.create(
            ID='door_ext',
            U=Q_(3.0, 'W / (m ** 2 * K)'),
            geometry=Geometry(t=Q_(8, 'cm'))
        )

    @staticmethod
    def _create_ca_int_wall_projF1():
        # create interior wall construction assembly "à la carte"
        ext_surf_film = SurfaceLayer.create(
            ID='ext_surf_film',
            geometry=Geometry(),
            heat_flow_direction=HeatFlowDirection.HORIZONTAL,
            Tmn=Q_(10, 'degC')
        )
        gypsum_layer_ext = BuildingComponent.create(
            ID='gypsum_layer_ext',
            geometry=Geometry(t=Q_(1.5, 'cm')),
            material=MaterialsShelf.load('gipspleister')
        )
        brick_layer = BuildingComponent.create(
            ID='brick_layer',
            geometry=Geometry(t=Q_(9, 'cm')),
            material=MaterialsShelf.load('blokken gebakken aarde, 1200 kg/m3')
        )
        gypsum_layer_int = BuildingComponent.create(
            ID='gypsum_layer_int',
            geometry=Geometry(t=Q_(1.5, 'cm')),
            material=MaterialsShelf.load('gipspleister')
        )
        int_surf_film = SurfaceLayer.create(
            ID='int_surf_film',
            geometry=Geometry(),
            heat_flow_direction=HeatFlowDirection.HORIZONTAL,
            Tmn=Q_(20, 'degC')
        )
        int_wall_projF1 = ConstructionAssembly.create(
            ID='int_wall_projF1',
            layers=[
                ext_surf_film,
                gypsum_layer_ext,
                brick_layer,
                gypsum_layer_int,
                int_surf_film
            ]
        )
        return int_wall_projF1
````

### 4. Create the Building
We will encapsulate the creation of the building inside a class called `House`. 
The following grand steps need to be taken:

1. In the `__init__`-method we load the construction assemblies by 
instantiating class `ConstructionAssemblies` and referring it to an instance 
attribute `self.constr_assem`. 

2. Next we declare the building and all of its sub-parts (namely the building 
entity, the ventilation zone and all the spaces of the building) as instance 
attributes of the class. 

3. Then, the building and all of its sub-parts are created inside the method 
`_create_building`, which will be called when instantiating class `House`.

4. Finally, all heated spaces inside the building are further configured in 
separate methods, that will also be called when instantiating class `House`. In
this step the building elements are added to each of the heated spaces. For 
unheated spaces it won't be necessary to add building elements, as these 
spaces don't lose heat; they only accept it from heated spaces. However, they
do play a role in the determination of the ventilation heat loss.

```python
    def __init__(self):
        # load construction assemblies
        self.constr_assem = ConstructionAssemblies()

        # declare the attributes of the building
        self.building = Building()
        self.building_entity = BuildingEntity()
        self.ventilation_zone = VentilationZone()
        self.kitchen_and_dining_room = HeatedSpace()
        self.living_room = HeatedSpace()
        self.bedroom_1 = HeatedSpace()
        self.bedroom_2 = HeatedSpace()
        self.bathroom = HeatedSpace()
        self.hall_way = UnheatedSpace()
        self.toilet = UnheatedSpace()

        # create building and configure all heated spaces
        self._create_building()
        self._config_kitchen_and_dining_room()
        self._config_living_room()
        self._config_bedroom_1()
        self._config_bedroom_2()
        self._config_bathroom()
```

To create the building layout in step 3, we have to go from top to bottom:

3.1. Create a `Building` instance. This instance takes the climatic design data
that is needed to perform the heating load calculations.

3.2. Add a single building entity to the building.

3.3. Add a single ventilation zone to this building entity.

3.4. Add the heated and unheated spaces to this ventilation zone.

```python
    def _create_building(self):
        # create building
        self.building = Building.create(
            ID='house',
            climate_data=ClimateDesignData(
                T_ext_d=Q_(-7.0, 'degC'),
                T_ext_an=Q_(10.0, 'degC'),
                T_ext_min=Q_(0.0, 'degC')
            )
        )
        # add building entity to building
        self.building_entity = self.building.add_building_entity(ID='house')
        # add ventilation zone to building entity
        self.ventilation_zone = self.building_entity.add_ventilation_zone(ID='house')
        # add heated spaces to ventilation zone
        self.kitchen_and_dining_room = self.ventilation_zone.add_heated_space(
            ID='kitchen_and_dining_room',
            height=Q_(3.0, 'm'),
            area=Q_(27.0, 'm ** 2'),
            T_int_d=Q_(20, 'degC'),
            V_exh=Q_(75.0, 'm ** 3 /hr'),
            V_sup=Q_(100.0, 'm ** 3 / hr')
        )
        self.living_room = self.ventilation_zone.add_heated_space(
            ID='living_room',
            height=Q_(3.0, 'm'),
            area=Q_(14.0, 'm ** 2'),
            T_int_d=Q_(20, 'degC'),
            V_sup=Q_(75.0, 'm ** 3 / hr')
        )
        self.bedroom_1 = self.ventilation_zone.add_heated_space(
            ID='bedroom_1',
            height=Q_(3, 'm'),
            area=Q_(18, 'm ** 2'),
            T_int_d=Q_(18.0, 'degC'),
            V_sup=Q_(72, 'm ** 3 / hr')
        )
        self.bedroom_2 = self.ventilation_zone.add_heated_space(
            ID='bedroom_2',
            height=Q_(3, 'm'),
            area=Q_(11, 'm ** 2'),
            T_int_d=Q_(18.0, 'degC'),
            V_sup=Q_(72, 'm ** 3 / hr')
        )
        self.bathroom = self.ventilation_zone.add_heated_space(
            ID='bathroom',
            height=Q_(3, 'm'),
            area=Q_(7.5, 'm ** 2'),
            T_int_d=Q_(24.0, 'degC'),
            V_exh=Q_(50.0, 'm ** 3 / hr'),
            V_trf=Q_(50.0, 'm ** 3 / hr'),
            T_trf=Q_(10.0, 'degC')
        )
        self.hall_way = self.ventilation_zone.add_unheated_space(
            ID='hallway',
            height=Q_(6, 'm'),
            area=Q_(11.5, 'm ** 2'),
            T_int_d=Q_(10, 'degC'),
            V_exh=Q_(169, 'm ** 3 / hr'),
            V_trf=Q_(75, 'm ** 3 / hr')
        )
        self.toilet = self.ventilation_zone.add_unheated_space(
            ID='toilet',
            height=Q_(3, 'm'),
            area=Q_(1.75, 'm ** 2'),
            T_int_d=Q_(10, 'degC'),
            V_exh=Q_(25, 'm ** 3 / hr'),
            V_trf=Q_(25, 'm ** 3 / hr')
        )
```

After all spaces are created, as the final grand step 4, we need to further 
configure the heated spaces by adding the building elements that surround them. 
Here, only the programming code for the configuration of the kitchen and dining 
room is shown.
```python
    def _config_kitchen_and_dining_room(self) -> None:
        # add exterior wall at the north side of the kitchen and dining room
        self.kitchen_and_dining_room.add_exterior_building_element(
            ID='ext_wall_north',
            area=(Q_(5.6, 'm'), Q_(3.0, 'm')),
            construction_assembly=self.constr_assem.ext_wall_WTCB1
        )
        # add wall adjacent to neighbouring house
        self.kitchen_and_dining_room.add_adjacent_building_element(
            ID='adj_wall_west',
            area=(Q_(4.65, 'm'), Q_(3.0, 'm')),
            construction_assembly=self.constr_assem.int_wall_WTCB1,
            kind_of_adjacent_space='unheated',
            T_adj=Q_(10, 'degC')
        )
        # add exterior wall at the east side of the kitchen and dining room
        ext_wall_east = self.kitchen_and_dining_room.add_exterior_building_element(
            ID='ext_wall_east',
            area=(Q_(6.225, 'm'), Q_(3, 'm')),
            construction_assembly=self.constr_assem.ext_wall_WTCB1
        )
        # add window and backdoor to east wall
        ext_wall_east.add_building_element(
            ID='window_east',
            area=(Q_(3.7, 'm'), Q_(2.7, 'm')),
            construction_assembly=self.constr_assem.window_ASHRAE5a
        )
        ext_wall_east.add_building_element(
            ID='door_east',
            area=(Q_(0.8, 'm'), Q_(2.7, 'm')),
            construction_assembly=self.constr_assem.ext_door
        )
        # add interior wall adjacent to toiletroom
        int_wall_toilet = self.kitchen_and_dining_room.add_adjacent_building_element(
            ID='int_wall_toilet',
            area=(Q_(3.025, 'm'), Q_(3, 'm')),
            construction_assembly=self.constr_assem.int_wall_PROJ1,
            T_adj=Q_(10, 'degC')
        )
        # add door to interior wall
        int_wall_toilet.add_building_element(
            ID='door_toilet',
            area=(Q_(0.8, 'm'), Q_(2.1, 'm')),
            construction_assembly=self.constr_assem.int_door
        )
        # add floor to kitchen and dining room
        self.kitchen_and_dining_room.add_ground_building_element(
            ID='floor',
            area=Q_(27.0, 'm ** 2'),
            construction_assembly=self.constr_assem.floor_WTCB4,
            A_slab=Q_(59.36, 'm ** 2'),
            P_slab=Q_(32.4, 'm'),
            z=Q_(0, 'm')
        )
        # add ceiling to kitchen and dining room
        self.kitchen_and_dining_room.add_adjacent_building_element(
            ID='ceiling',
            area=Q_(27.0, 'm ** 2'),
            construction_assembly=self.constr_assem.ceiling_WTCB13,
            kind_of_adjacent_space='heated',
            T_adj=Q_(18.60, 'degC')  # area-weighted average adjacent space temperature
        )
```

### 5. Run the Heating Load Calculation
To run the heating load calculation, it suffices to instantiate class `House`.
As soon as the instantiation is finished, the results of the heating load 
calculations are available.

```python
def main():
    house = House()
    with pd.option_context(
        'display.max_rows', None,
        'display.max_columns', None,
        'display.width', None
    ):
        print(house.ventilation_zone.get_summary())
        print()
        print(house.building_entity.get_summary())


if __name__ == '__main__':
    main()
```

To get the results of the heated spaces, you can call `get_summary()` on the 
`ventilation_zone` attribute of the `House` instance. This will return the 
results of all heated spaces that are part of the ventilation zone in a 
Pandas `DataFrame` object.

```text
heated space                Q transmission [kW]  Q ventilation [kW]  Q heating-up [kW]  Q total [kW]
0  kitchen_and_dining_room                2.366               0.422                0.0         2.787
1              living_room                2.217               0.308                0.0         2.525
2                bedroom_1                1.349               0.290                0.0         1.639
3                bedroom_2                0.581               0.231                0.0         0.812
4                 bathroom                0.626               0.288                0.0         0.914
```

Getting the global results of the ventilation zone, happens in a similar way:

```text
  ventilation zone  Q transmission [kW]  Q ventilation [kW]  Q heating-up [kW]  Q total [kW]
0            house                7.139                 1.2                0.0         8.339
```
