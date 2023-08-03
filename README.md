# Project python-hvac

A Python package for HVAC engineering:
- Cooling and heating load calculations
- Estimation of energy needs and energy consumption
- Design and analysis of air and hydronic distribution circuits
- Air conditioning process calculations
- Simulation and analysis of vapor compression machines and VRF-systems
- Sizing of refrigerant piping

Example scripts and notebooks can be found under `docs/examples`. More examples
will be added over time.

> [!INFO]
> This project replaces the older project `HVAC`, which is no longer under 
> development and will be no longer maintained.

## Overview
At this moment following sub-packages are part of `hvac`:

**air_conditioning**<br>
A package about air-conditioning processes and systems, including VAV-systems.

**cooling_load_calc**<br>
A package for doing cooling load calculations of a building, based upon ASHRAE's
RTS-method. It uses a lumped linear thermal network model to take the thermal 
inertia of flat, opaque building components into account, and also for the simple
modeling of the internal thermal mass of a space. It uses a sub-package
**climate**, also included in main package `hvac`, for developing daily 
temperature profiles from a limited number of climate design parameters. This 
package also contains a sub-package **sun** for calculating solar irradiance and
irradiation on surfaces based on the ASHRAE Clear-Sky Model or based on TMY-data.
For calculating the position of the sun, the third-party library `astral` is 
used.

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
This packages implements a number of correlations for calculating convection
heat transfer coefficients and friction factors for different geometries. Most 
correlations were taken from the book "Introduction to Engineering Heat Transfer" 
by G.F. Nellis and S.A. Klein (Cambridge University Press).

Also part of this package is the sub-package **heat_exchanger**. It implements
the effectiveness-NTU method for both dry and wet external heat transfer 
surfaces. It also contains rating and sizing procedures for some 
fin-tube heat exchangers (see `heat_transfer.heat_exchanger.fin_tube`), which 
are applicable to single-phase fluid flow, and which are implementations of 
the procedures described in the book "Fundamentals of Heat Exchanger Design" by 
R. K. Shah and D. P. Sekulic (John Wiley & Sons, 2003).

**vapor_compression**<br>
Furthermore, sub-package `heat_transfer.heat_exchanger.fin_tube` also contains 
a model for a plain fin-tube, counter-flow air evaporator and similar air 
condenser. These models are used, together with a model for a fixed or variable 
speed compressor, in the module `machine_bis` of sub-package `vapor_compression`. 
Purpose of this module is to simulate the steady-state performance of a 
single-stage air conditioning machine (vapor compression machine, heat pump). 
Multiple examples of such simulations have been included in the folder 
`vapor_compression` of `docs/examples`.

**vrf_system**<br>
This sub-package can be used to model a VRF-system in relation to the building
in which it is installed, with the purpose to estimate its energy consumption
during a typical heating or cooling season. An example of such an energy 
consumption estimation and a comparison between different candidate VRF-systems 
in relation to the same building have been included in the folder `vrf_system` 
of `docs/examples`.

Besides the aforementioned application-oriented packages, `hvac` also includes a 
number of more basic sub-packages which are used throughout the code of `hvac`:

First of all, `hvac` is entirely based on the third-party library *Pint* for
representing physical quantities. Module `pint_setup.py` does the necessary
setup for using Pint's `Quantity` class throughout the package.

Sub-package `fluids` contains a number of modules with classes that act like more
object-oriented wrappers around the third-party libraries `CoolProp` and `aipws`. 
The most principal classes of `fluids`, that are used throughout the code, are 
`Fluid` and `HumidAir`, which encapsulate CoolProp's API, mainly to allow the 
use of `Quantity` objects.

Finally, sub-package `charts` contains a thin wrapper around the third-party 
library `matplotlib`. It also has two modules for plotting refrigeration cycles 
on the log-p-h diagram of a refrigerant and for plotting air conditioning 
processes on a psychrometric chart.
