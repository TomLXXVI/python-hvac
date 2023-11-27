# python-hvac

A Python library for HVAC engineering:
- cooling and heating load calculations of buildings
- estimation of energy needs and energy consumption
- design and analysis of air and hydronic distribution systems
- air conditioning process calculations
- simulation and analysis of heat exchangers (limited for now to plain fin-and-tube heat exchangers)
- simulation and analysis of single-stage vapor compression machines and VRF systems
- sizing of refrigerant piping
- correlations for convection heat transfer coefficients and friction factors

Example scripts and notebooks can be found under `docs/examples`. More examples
will be added over time.

> [!NOTE]
> This project replaces the older project `HVAC`, which is no longer under 
> development and will be no longer maintained.

## Overview
At this moment the following subpackages are part of `hvac`:

**air_conditioning**<br>
A package about air-conditioning processes and systems.

**cooling_load_calc**<br>
A package for doing cooling load calculations of a building, based upon ASHRAE's
Radiant-Times-Series (RTS) method. 
It uses a lumped linear thermal network model to take the thermal inertia of 
flat, opaque building components into account. The same network model is also 
used to implement a simplified model for the internal thermal mass of a space. 

This package uses the subpackage **climate**, also included in the main package
`hvac`, for the development of daily temperature profiles from a limited number 
of climate design parameters. 
Package **climate** also contains a subpackage **sun** for calculating solar 
irradiance and irradiation on surfaces, either based on ASHRAE's 
*Clear-Sky Model*, or based on TMY-data. For calculating the position of the 
sun, the third-party library `astral` is used.

**heating_load_calc**<br>
A package for doing heating load calculations of a building based on the method 
of European standard EN 12831-1 (July 2017).

**energy_estimation**<br>
Estimation of the energy consumption of a building using the bin table method.

**fluid_flow**<br>
A package about fluid flow in pipes and ducts and for designing and analyzing 
pipe and duct network systems. It includes a large number of fittings taken 
from Crane Technical Paper No. 410M and SMACNA's HVAC Systems Duct Design Manual.

**refrigerant_piping**<br>
Package for sizing refrigerant lines (suction line, discharge line, and liquid
line) between an outdoor unit and indoor unit of an air conditioning system.

**heat_transfer**<br>
This package implements a number of correlations for calculating convection
heat transfer coefficients and friction factors for different geometries. Most 
correlations were taken from the book "Introduction to Engineering Heat Transfer" 
by G. F. Nellis and S. A. Klein (Cambridge University Press).

**heat_exchanger**<br>
Originally (and still) part of the subpackage **heat_transfer**, but a new, 
rewritten package has also been added directly as a subpackage of the main package
**hvac**. It implements the effectiveness-NTU method for both dry and wet 
external heat transfer surfaces, and it contains rating and sizing procedures 
for some fin-tube heat exchangers, which are applicable to single-phase fluid 
flow, and which are implementations of the procedures described in the book 
"Fundamentals of Heat Exchanger Design" by R. K. Shah and D. P. Sekulic 
(John Wiley & Sons, 2003).

**vapor_compression**<br>
Furthermore, subpackage `heat_transfer.heat_exchanger.fin_tube` also contains 
a model for a plain fin-tube, counter-flow air evaporator and similar air 
condenser. These models are used together with a model for a fixed or variable 
speed compressor in the module `machine` of subpackage `vapor_compression`. 
Purpose of this module is to simulate the steady-state performance of a 
single-stage vapor compression machine (air conditioning machine, heat pump). 
Several application examples have been included in the folder `vapor_compression` 
of `docs/examples`.

**vrf_system**<br>
This subpackage can be used to model a VRF-system in relation to the building
in which it is installed, with the purpose to estimate its energy consumption
during a typical heating or cooling season. An example of such an energy 
consumption estimation and a comparison between different candidate VRF-systems 
in relation to the same building have been included in the folder `vrf_system` 
of `docs/examples`.

Besides the aforementioned application-oriented packages, `hvac` also includes a 
number of more basic subpackages which are used throughout the modules of 
`hvac`:

First of all, `hvac` is heavily based on third-party library `pint` for
working with physical quantities in Python. Module `pint_setup.py` inside `hvac`
does the necessary setup for using Pint's `Quantity` class throughout the 
package. When you write a script using package `hvac`, simply write `from hvac
import Quantity` to work with physical quantities in your script.

Subpackage `fluids` contains a number of modules with classes that act like
object-oriented wrappers around third-party libraries `CoolProp` and `aipws`. 
The principal classes of `fluids`, that are used extensively throughout 
the code, are `Fluid` and `HumidAir`, which encapsulate CoolProp's API and
allow to accept `Quantity` objects.

Finally, subpackage `charts` contains a package `matplotlibwrapper`, being a 
tiny wrapper around third-party library `matplotlib`, meant to ease the drawing
of some frequently used chart types. It also contains a module to plot 
refrigeration cycles on the log-p-h diagram of a refrigerant and a module to 
plot air conditioning processes on a psychrometric chart.
