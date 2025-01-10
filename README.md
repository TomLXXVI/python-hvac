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
A package about air-conditioning processes and systems that can be used to 
design and analyze single-zone or multi-zone CAV and VAV air conditioning
systems.

**air_diffusion**<br>
A package about supplying air to a room. Implements free air jets and some
procedures for designing air supply to a room which come from *Awbi, H. B. (2003).
VENTILATION OF BUILDINGS. Taylor & Francis. Chapter 6.*

**cooling_load_calc**<br>
A package for doing cooling load calculations of a building, based upon ASHRAE's
Radiant-Times-Series (RTS) method. 
It uses a lumped capacitance model to take the thermal inertia of flat, opaque 
building components and the interior thermal mass of a zone into account. 

**sun**<br>
This package is about solar radiation on surfaces. It implements a number of
sky-models (isotropic, anisotropic Perez, anisotropic HDKR, KT) to estimate the 
solar radiation incident on exterior surfaces.
It can generate solar radiation data based on the clear-sky model or it can use
TMY-data of a certain geographic location .
Subpackage *sun* implements the equations in chapters 1 and 2 of *Duffie, 
J. A., Beckman, W. A., & Blair, N. (2020). SOLAR ENGINEERING OF THERMAL 
PROCESSES, PHOTOVOLTAICS AND WIND. John Wiley & Sons.*

**heating_load_calc**<br>
A package for doing the heating load calculation of a building based on the 
method of *European standard EN 12831-1 (July 2017)*.

**energy_estimation**<br>
Estimation of the energy consumption of a building using the bin table method.

**fluid_flow**<br>
A package about fluid flow in pipes and ducts and for designing and analyzing 
pipe and duct network systems. It includes a large number of fittings taken 
from *Crane's TECHNICAL PAPER NO. 410M* and *SMACNA's HVAC SYSTEMS DUCT DESIGN 
MANUAL*.

**refrigerant_piping**<br>
Package for sizing refrigerant lines (suction line, discharge line, and liquid
line) between an outdoor unit and indoor unit of an air conditioning system.

**heat_transfer**<br>
This package implements a number of correlations for calculating convection
heat transfer coefficients and friction factors for different geometries. Most 
correlations were taken from *Nellis G. F. , & Klein S. A.  (2021). 
INTRODUCTION TO ENGINEERING HEAT TRANSFER. Cambridge University Press*.

**heat_exchanger**<br>
Contains two subpackages: **recuperator** and **regenerator**.

Subpackage *recuperator* is about heat exchangers in which the hot and cold 
fluid are separated by a heat transfer wall. 
It implements the effectiveness-NTU method for both dry and wet air 
cooling/heating coils. 
It contains (at this moment) some models of continuous fin-tube heat exchangers
(air condenser, air evaporator, air-to-water cooling coil) for analysis and 
simulation. Implementations are based on solving methods described in *Shah, 
R. K. , & Sekulic, D. P.  (2003). FUNDAMENTALS OF HEAT EXCHANGER 
DESIGN. John Wiley & Sons*.

Subpackage *regenerator* implements the effectiveness-NTU method for a 
counterflow rotary regenerator (e.g. a sensible heat recovery wheel) as outlined
in *Shah, R. K. , & Sekulic, D. P.  (2003). FUNDAMENTALS OF HEAT 
EXCHANGER DESIGN. John Wiley & Sons*, Chapter 5: *Thermal Design Theory for
Regenerators*.

**vapor_compression**<br>
Contains models to represent fixed- and variable speed compressors, and a model
of a single-stage vapor compression machine, that integrates the models for the
air evaporator and condenser in subpackage *heat_exchanger*.
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

**radiant_emitter**<br>
This subpackage implements two types of radiant heat emitters. Class 
`PanelRadiator` is for panel radiators and class `RadiantFloorPanel` is for 
floor heating panels. Both classes can be used for design and analysis.

**control**<br>
This subpackage implements different types of controllers: on/off, PID, and PWM.
Together with class `VariableTemperatureZone` in subpackage `cooling_load_calc`,
these controller classes can be used to simulate a zone with a temperature 
controlled heating/cooling system. Script `example_02.py` in 
`docs/examples/heating_load_calc` shows an example with an on/off controller.

**kitchen_ventilation**<br>
This subpackage implements the detailed calculation method according to European
standard EN 16282-1:2017 *Equipment for commercial kitchens - Components for 
ventilation in commercial kitchens - Part 1: General requirements including 
calculation method* for the design calculation of the extraction air volume flow
rates in commercial kitchens. It also includes two quick methods that can be 
used for a preliminary calculation.

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
A new module `fluid_experimental` has been added to subpackage `fluids` which has a better structured and understandable user interface, especially for defining mixtures. This module has not been integrated yet in other modules or subpackages, while it's still in an experimental stage.

Finally, subpackage `charts` contains a package `matplotlibwrapper`, being a 
tiny wrapper around third-party library `matplotlib`, meant to ease the drawing
of some frequently used chart types. It also contains a module to plot 
refrigeration cycles on the log-p-h diagram of a refrigerant and a module to 
plot air conditioning processes on a psychrometric chart.
