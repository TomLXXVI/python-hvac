# Sizing or Rating a Single-Stage Vapor Compression Machine

## Preliminary Calculations for Sizing

 Example: 

```docs/examples/vapor_compression/01_prelim_calc.ipynb```

---

Problem outline: An air stream with a given state at the entrance and a given 
mass flow rate must be cooled to a certain final temperature in the evaporator 
of a single-stage vapor compression machine.

1. Select the refrigerant.

2. From the problem statement, the heat rate to be absorbed by the refrigerant 
   in the evaporator can immediately be determined.

3. Select an **evaporation temperature** and determine the corresponding 
   saturation pressure (the lower this value with respect to the entering air 
   temperature, the smaller the evaporator can be; however a value below the 
   freezing point of water, will cause frost formation on the air-side heat 
   transfer surface of the evaporator, which will block the air flow and 
   requires a periodic defrosting cycle).

4. Select a **condensation temperature** and determine the corresponding 
   saturation pressure (the higher this value with respect to the entering 
   outdoor air temperature, the smaller the condenser can be).

5. To find the mass flow rate of refrigerant needed to absorb the heat from the 
   air stream, we need to know the state (especially the enthalpy) of the 
   refrigerant entering and leaving the evaporator. Therefore, we need to select
   a certain **degree of superheat** of the refrigerant leaving the evaporator 
   and a certain **degree of subcooling** of the refrigerant leaving the 
   condenser. A sufficient degree of superheat will prevent that any liquid is 
   still present in the refrigerant which could damage the compressor. Some 
   degree of subcooling is needed to prevent that any vapor is still present in 
   the refrigerant which could disturb the operation of the expansion device. 
   Moreover, subcooling also has a beneficial effect on the refrigeration effect. 
   The expansion process is considered to be isenthalpic, so that the enthalpy 
   of the refrigerant entering the evaporator will be equal to the enthalpy of 
   refrigerant leaving the condenser. Dividing the required heat absorption rate
   (i.e. the required refrigeration capacity) by the enthalpy difference of the 
   leaving and entering refrigerant in the evaporator, gives us the required 
   mass flow rate of refrigerant.

6. To determine the mass flow rate of air through the condenser, we first need 
   to know the heat rate that must be rejected to the outdoor air stream, being 
   the sum of the heat rate absorbed in the evaporator and the mechanical power 
   delivered by the compressor to the refrigerant. To determine the compressor 
   power, a compressor selection software-program can be used to select a 
   compressor that might be suitable for the given application.

7. To determine the mass flow rate of air through the condenser, we first select
   a target temperature rise of the air stream through the condenser and then we
   can determine the temperature of air leaving the condenser. Next, we can 
   determine the corresponding mass flow rate of air for the given heat 
   rejection rate based on the mean specific heat of humid air.

8. Finally, based on a selected air face velocity, the frontal areas of the 
   evaporator and condenser can be determined.

## Rating of a Plain Fin-Tube Counterflow Air Evaporator

Example: 

```docs/examples/vapor_compression/02_air_evaporator.ipynb```

---

From the problem statement or our preliminary calculations we have:

- the state of air entering the evaporator
- the mass flow rate of air through the evaporator
- an estimation of the state of refrigerant entering the evaporator
- the degree of refrigerant superheating set on the expansion device
- an initial guess for the needed mass flow rate of refrigerant

We select a candidate air evaporator with known dimensions (the width and height
of the frontal area and the number of rows) and with a known core geometry. 
Package `hvac.heat_exchanger.fintube.continuous_fin` contains one model of an 
air evaporator and one model of an air condenser, both being plain fin-tube 
counter-flow heat exchangers, where all refrigerant is distributed into the heat
exchanger across the tubes of the first row and leaves the heat exchanger 
through the corresponding tubes of the last row (in reality, heat transfer will 
be a combination of cross- and counter-flow).

We set the known operating conditions to the air evaporator model. Using an 
iterative solving technique, the mass flow rate of refrigerant through the 
evaporator is calculated, such that the degree of superheat of the refrigerant, 
which is set by the expansion device, is attained. Once the mass flow rate of 
refrigerant through the evaporator is determined, all other performance 
parameters of the evaporator can also be determined.

Now we can play with the dimensions (number of rows, frontal area) and the core 
geometry parameters to find an appropriate evaporator, i.e. to get the desired 
state of the air leaving the evaporator. However, in a vapor compression machine,
it will be the simultaneous interaction between evaporator, compressor, 
condenser, and expansion device that will ultimately determine the actual state 
of air leaving the evaporator.

## Rating of a Plain Fin-Tube Counterflow Air Condenser

Example:

```docs/examples/vapor_compression/03_air_condenser.ipynb```

---

From the preliminary calculations and the rating of the evaporator we have:

- the state of air entering the condenser
- the mass flow rate of air through the condenser
- the state of refrigerant entering the condenser
- the mass flow rate of refrigerant through the condenser, determined when 
  rating the evaporator

The same procedure as with the air evaporator can be followed. By changing the 
number of rows, the core geometry and/or changing the mass flow rate of air 
through the condenser, we can influence the performance of the condenser so 
that we get at the required heat rejection rate.

## Rating of a Single-Stage Vapor Compression Machine

Examples:

```
docs/examples/vapor_compression/05_vcm_rating.ipynb
docs/examples/vapor_compression/07_multi_vcm_analysis.py
```
---

A single-stage vapor compression machine is composed of a:

- an air evaporator
- an air condenser
- a compressor
- an expansion device

To configure the single-stage vapor compression machine model, we first need a 
model of the evaporator, the condenser, and the compressor. The expansion device 
is not modeled. It is assumed that the expansion device is able to regulate the 
mass flow rate of refrigerant into the evaporator, so that the set degree of 
superheat is maintained under any possible set of operating conditions. The 
model of an air evaporator and condenser were already introduced above. We can 
now use these models to configure the model of the single-stage vapor 
compression machine.  How the model of the compressor can be created, is 
demonstrated in notebook `04_compressor_ipynb`. 

Once the single-stage vapor compression machine model is configured, we can set 
the operating conditions on this machine and determine the steady-state 
performance of the machine under these operating conditions. To do this, an 
iterative solving technique is used similar to the analysis procedure explained 
in Chapter 14 of the book *Stoecker, W. F., & Jones, J. W. (1982). Refrigeration
and air conditioning. McGraw-Hill International Editions, Mechanical Technology 
Series.* 
(However, in the implementation of the `SingleStageVaporCompressionMachine` 
model class, the compressor correlation for the refrigeration capacity, obtained
from the selection software-program of the compressor manufacturer, is not used, 
as it is implicitly based on a certain, fixed degree of superheat of refrigerant
leaving the evaporator and a certain, fixed degree of subcooling of the 
refrigerant leaving the condenser. Only the compressor correlation for the mass 
flow rate of refrigerant and the compressor power are used.)  

Another purpose of rating (analyzing) a single-stage vapor compression machine 
could be to investigate for which configuration the desired state of air leaving 
the evaporator can be accomplished. To find this configuration for a given core 
geometry of the evaporator and the condenser, and a given compressor, we can 
play with:

- the number of rows of the evaporator
- the number of rows of the condenser
- the mass flow rate of air through the condenser
- the compressor speed in case the compressor has a variable speed drive

## Balancing a Variable Speed, Single-Stage Vapor Compression Machine

Example:

```docs/examples/vapor_compression/06_vcm_rating_bis.py``` 

---

Suppose we have a single-stage vapor compression machine with a variable speed 
drive working under a given set of operating conditions (state and mass flow 
rate of air entering the evaporator, state and mass flow rate of air entering 
the condenser, fixed superheat setting on the expansion device). With the model 
of the vapor compression machine, it is also possible to determine the 
compressor speed at which a certain evaporation and condensation temperature 
will be attained. This is implemented in the `SingleStageVaporCompressionMachine`
class by the method `balance_by_speed`.
