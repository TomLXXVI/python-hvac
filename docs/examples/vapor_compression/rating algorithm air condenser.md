# Description of rating algorithm for an air condenser

Source code:

```hvac.heat_exchanger.fintube.continuous_fin.air_condenser.rating```

Example:

```docs/examples/vapor_compression/03_air_condenser.ipynb```

---

The air condenser is modeled as a **plain fin-tube counter-flow heat exchanger** 
and is characterized by the following parameters:

- the frontal width of the heat exchanger core
- the frontal height of the heat exchanger core
- the number of rows, which, together with the longitudinal pitch, determine the 
  total flow length of the condenser
- the transversal pitch, i.e. the distance between 2 tubes of the same row
- the longitudinal pitch, i.e. the distance between 2 tubes of adjacent rows
- the inner diameter of the tubes
- the outer diameter of the tubes
- the thickness of the plain fins
- the fin density, i.e. the number of fins per unit length of tube
- the thermal conductivity of the fin material

Three regions can be distinguished on the refrigerant side of the condenser:

- the desuperheating region, where superheated refrigerant vapor from the 
  compressor is cooled to saturated vapor
- the condensing region, where saturated vapor is condensed to saturated liquid
- the subcooling region, where liquid refrigerant is cooled further below the
  saturation temperature that corresponds with the condenser pressure

The performance of the condenser depends on:

- the state of air at the inlet of the condenser
- the mass flow rate of air through the condenser
- the state of refrigerant at the inlet of the condenser
- the mass flow rate of refrigerant through the condenser

To fully determine the performance of the condenser, we start at the subcooling 
region of the condenser. The solving process is based on a root finding 
algorithm. The aim is to find the flow lengths of the subcooling, the condensing,
and the desuperheating region of the condenser, such that the sum of these flow 
lengths equals the total flow length of the condenser, i.e. the difference 
between the sum of these three flow lengths and the actual, total flow length of
the condenser must be (nearly) zero. The algorithm starts with an initial guess 
for the subcooling flow length. When a solution is found for the subcooling 
region, the inlet and outlet states of air at the condensing and desuperheating 
regions can be determined. It is then, for the given mass flow rate of 
refrigerant and air, also possible to find the condensing flow length, needed to
condense the refrigerant from saturated vapor into saturated liquid, and the 
desuperheating flow length, needed for desuperheating the hot, compressed 
refrigerant vapor entering the condenser into saturated vapor. If no solution 
can be found based on the assumption that subcooling of refrigerant happens in 
the condenser, the solving routine will set the subcooling flow length to zero 
and will continue to try to find the flow lengths of only the condensing and 
the desuperheating region of the condenser such that the sum of these flow 
lengths equals the total flow length of the condenser. If again no solution can 
be found, an error is raised to signal that the condenser cannot be rated under 
the given set of operating conditions.

## Subcooling Region

The state of refrigerant at the inlet of the subcooling region is known. It is 
saturated liquid (vapor quality 0 %) of which the pressure is assumed to be 
equal to the pressure of the refrigerant at the inlet of the condenser. In other
words, it is assumed that the refrigerant pressure in the condenser remains 
constant (i.e. the condenser pressure). Air enters the condenser at the 
subcooling region. The state of air at the condenser inlet is known. The mass 
flow rate of air and the mass flow rate of refrigerant are known. With these 
data and the current value for the subcooling flow length it is possible to 
solve for the heat transfer rate across the subcooling region and to determine 
the state of the air leaving the subcooling region and entering the condensing
region. An iterative solution technique is used to accomplish this:

1. An initial guess is made of the state of the refrigerant leaving the 
   condenser.

2. With this the heat transfer rate in the subcooling region can be determined.

3. Knowing the heat transfer rate to the air, the temperature of the air leaving
   the subcooling region and entering the condensing region can be determined.
   As the air heating process is only a sensible process, the humidity ratio of
   the leaving air remains equal to the humidity ratio of the air entering the
   condenser.

4. With the known inlet states of air and refrigerant and the outlet states of 
   refrigerant and air being determined, the mean state of the air and the 
   refrigerant across the subcooling region can be determined.

5. Based on the current value for the subcooling flow length, the mass flow 
   rates and the mean states of air and refrigerant, the overall heat transfer 
   conductance of the heat exchanger core for the subcooling region can be 
   determined.

6. Using the effectiveness-NTU method, a new value for the heat transfer rate 
   across the subcooling region can be determined.

7. With this new value of the heat transfer rate, also a new state of the 
   refrigerant leaving the condenser can be determined.

8. Now, the deviation is calculated between the refrigerant temperature in the
   new state and in the previous state.

9. If the deviation is greater than the allowed, small deviation (i.e. the 
   tolerance), steps 2 to 8 are repeated. When the tolerance is met, the state
   of the air leaving the subcooled region is determined and returned.

## Condensing Region

The state of refrigerant at the inlet and at the outlet of the condensing region
are already known from the start: saturated vapor and saturated liquid. 
Therefore, as the mass flow rate of refrigerant is given, the heat transfer rate
across the condensing region can immediately be solved for. Using the mean 
specific heat of humid air, and as the mass flow rate of air is given, it is 
also possible to immediately calculate the temperature rise of air across the 
condensing region of the condenser.

After the subcooling region has been solved, we know the state of air entering 
the condensing region. As the temperature rise of air across the condensing 
region is already determined, and as the humidity ratio of air across the 
condenser remains constant and equal to the known humidity ratio of air at the 
condenser inlet, we can also determine the state of air leaving the condensing 
region and entering the desuperheating region. This allows to determine the mean
state of air and refrigerant across the condensing region.

Now, the required flow length of the condensing region can be determined, so 
that for the given mass flow rate of refrigerant, the refrigerant has been 
transformed from saturated vapor into saturated liquid. This is done with 
a root finding algorithm. The root finding algorithm searches for a flow length
so that the corresponding heat transfer rate equals the already known heat 
transfer rate in the condensing region of the condenser.

1. The root finding algorithm determines a value for the condensing flow length
   between a given minimum and maximum value.

2. With this flow length, the heat transfer characteristics in the condensing 
   region are determined.

3. Using the effectiveness-NTU method, the heat transfer rate across the 
   subcooling region can be determined.

4. This heat transfer rate should be equal to the already known heat transfer
   rate. In other words, the deviation between these two values should be close
   to zero (this is what's root finding is about).

5. When the root finding algorithm finds a flow length for which the deviation 
   is zero (or very close to it), it returns this flow length, which is then the
   flow length of the condensing region that corresponds with the entering air
   state coming from solving the subcooling region first.


## Desuperheating Region

The same procedure as for the condensing region can be followed. The state of
refrigerant at the inlet and at the outlet of the desuperheating region are 
already known from the start. At the outlet the refrigerant must be saturated 
vapor (vapor quality 100 %). So, as in the case of the condensing region, the 
heat transfer rate and the temperature rise of the air across the desuperheating
region can be immediately solved for. Consequently, the state of air leaving the
condenser is now also determined, as the state of air entering the 
desuperheating region has been determined when solving for the condensing region
first. Again, the mean state of air and refrigerant across the desuperheating 
region are determined and the required flow length of the desuperheating region 
is determined by applying a similar root finding algorithm.
