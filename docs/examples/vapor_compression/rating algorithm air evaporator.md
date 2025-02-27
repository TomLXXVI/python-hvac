# Description of the rating algorithm for an air evaporator

Source code:

```hvac.heat_exchanger.fintube.continuous_fin.air_evaporator```

Example:

```docs/examples/vapor_compression/02_air_evaporator.ipynb```

---

The air evaporator is modeled as a **plain fin-tube counterflow heat exchanger**
and is characterized by the following parameters:

- the width of the heat exchanger core or length of the tubes
- the height of the heat exchanger core
- the number of rows, which, together with the longitudinal pitch, determine 
  the total flow length of the evaporator
- the transversal pitch, i.e. the distance between 2 tubes of the same row
- the longitudinal pitch, i.e. the distance between 2 tubes of adjacent rows
- the inner diameter of the tubes
- the outer diameter of the tubes
- the thickness of the plain fins
- the fin density, i.e. the number of fins per unit length of tube
- the thermal conductivity of the fin material

Two regions can be distinguished on the refrigerant side of the evaporator :

- the boiling region, where a liquid/vapor mixture of refrigerant enters the 
  evaporator and refrigerant leaves as a saturated vapor

- the superheating region, where refrigerant enters as saturated vapor and 
  leaves as superheated vapor. 


The performance of the evaporator depends on:

- the state of air at the inlet of the evaporator (temperature and humidity)
- the mass flow rate of air through the evaporator
- the state of refrigerant at the inlet of the evaporator (a saturated two-phase
  mixture)
- the degree of superheat of the refrigerant, which is a setting of the 
  expansion device

The expansion device will determine the mass flow rate of refrigerant through 
the evaporator to maintain the set degree of superheat. Once the mass flow rate
of refrigerant is determined, the performance of the evaporator can be fully 
determined. The determination of the mass flow rate of refrigerant is an 
iterative process that starts with an initial guess of this mass flow rate.

## Superheating Region

First, for the mass flow rate of the refrigerant in the current iteration, the 
flow length of refrigerant is determined needed to superheat the refrigerant 
from its saturated vapor state to the degree of superheat set on the expansion 
device. 

The refrigerant enters the superheating region of the evaporator as a saturated 
vapor (vapor quality 1.0). The pressure of the saturated vapor equals the known 
pressure of the refrigerant at the evaporator inlet (i.e. the evaporation 
pressure). 
At the evaporator inlet the refrigerant is a known saturated two-phase mixture of 
liquid and vapor, so the evaporation temperature/pressure is also known. Between
the evaporator inlet section and the section where all saturated liquid has been
transformed into saturated vapor, the refrigerant is boiling. It is assumed that
the pressure in the boiling region of the evaporator remains constant. 
At the evaporator outlet the refrigerant is superheated to the degree set on the
expansion device. Again, it is assumed that in the superheating region of the 
evaporator the refrigerant pressure remains constant. Knowing the state of the 
saturated vapor at the inlet of the superheating region of the evaporator and 
the state of superheated vapor leaving the evaporator (the degree of superheat 
is given), the heat transfer rate in the superheating region can be determined 
for the given mass flow rate of refrigerant.

On the air-side of the superheating region of the evaporator, the state of the 
entering air is given, and it is further assumed that heat transfer is only 
sensible. This implies that the humidity ratio of the air across the 
superheating region remains constant and equal to the humidity ratio of the air 
entering the evaporator. As the heat transfer rate and the air mass flow rate 
are known, the enthalpy of the air leaving the superheating region can
be determined. Humidity ratio and enthalpy fix the state of air leaving the 
superheating region of the evaporator.

Based on the inlet and outlet states of the refrigerant and the air and the mass
flow rates of the refrigerant and the air, the mean refrigerant and mean air 
state across the superheating region can be determined. 
Now, for a given flow length of refrigerant, the overall heat transfer 
conductance across the superheating region of the evaporator can be determined. 
Using the effectiveness-NTU method applied to the counterflow heat exchanger, 
the heat transfer rate across the superheating region can be determined for the 
given flow length of refrigerant. So, it comes down to finding the flow length 
for which the heat transfer rate becomes equal to the heat transfer rate we 
already know. This can be implemented with a root finding algorithm.

## Boiling Region

Once the flow length of the superheating region is known, the available flow 
length of the boiling region of the evaporator is also determined. 

At the refrigerant-side of the evaporator the state of refrigerant at the inlet 
is given. At the exit of the boiling region the refrigerant is always saturated 
vapor. 

Now, a new value for the mass flow rate of refrigerant can be determined such 
that along the available flow length for boiling the refrigerant is transformed 
from the known inlet state into the saturated vapor state. 

To get at this new value for the refrigerant mass flow rate an iterative 
(sub)routine is used. Starting with the initially guessed value of the mass flow
rate and knowing the states of the refrigerant at the inlet and outlet of the 
boiling region, the heat transfer rate across the boiling region of the 
evaporator can be calculated. 

On the air-side of the boiling region the inlet air state equals the state of 
air leaving the superheating region. Knowing the mass flow rate of air and the 
heat transfer rate, the enthalpy of air leaving the evaporator can be determined.
Initially, we can assume that this leaving air is fully saturated (RH 100 %). 

Knowing the refrigerant and air states at the inlet and outlet of the boiling 
region and the mass flow rates of refrigerant and air, the mean state of 
refrigerant and air across the boiling region can be determined and the heat 
transfer characteristics for the boiling region can be calculated. It is assumed
that on the air-side of the boiling region, water vapor from the humid air is 
condensing on the external heat transfer surface. Using the effectiveness-NTU 
method for a wet surface applied to a counterflow heat exchanger, a new value of
the heat transfer rate and a new state of the air leaving the evaporator can be 
determined. As the inlet state and the saturated vapor state of refrigerant 
remain fixed, the new value of the heat transfer rate also leads to a new value 
for the mass flow rate of refrigerant. 

The iterative subroutine is repeated until the new value becomes (almost) equal 
to the previous value. Once the difference has become small enough (i.e. the 
difference is equal or smaller than the accepted tolerance), the last calculated
value can be considered as the mass flow rate of refrigerant that is needed to 
transform the refrigerant entering the evaporator into saturated vapor within 
the available flow length for boiling. 

This new mass flow rate is now used to recalculate the flow length of the 
superheating region of the evaporator. The iterative main routine is repeated in
the same way as the iterative subroutine for the boiling region. Once the 
difference between the last and the previous value of the refrigerant mass flow 
rate has become small enough (i.e. the difference is equal or smaller than the 
accepted tolerance), the operation of the air evaporator is fully determined. 
The result of the calculations contains, besides the final mass flow rate of 
refrigerant, the state of air leaving the evaporator, the total heat transfer 
rate across the boiling and superheating region of the evaporator, the 
effectiveness of the air evaporator, the flow length of the evaporator needed to
superheat the refrigerant, and also the pressure drop of air across the external
heat transfer surface of the evaporator.
