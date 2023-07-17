# Rating of an air condenser

Source code:

```hvac.heat_transfer.heat_exchanger.fin_tube.air_condenser.rating```

Example:

```docs/examples/vapor_compression/03_air_condenser.ipynb```

---

The air condenser is modeled as a **plain fin-tube, counterflow heat exchanger** and is characterized by the following parameters:

- the width of the heat exchanger core or length of the tubes
- the height of the heat exchanger core
- the number of rows, which, together with the longitudinal pitch, determine the total flow length of the condenser
- the transversal pitch, i.e. the distance between 2 tubes of the same row
- the longitudinal pitch, i.e. the distance between 2 tubes of adjacent rows
- the inner diameter of the tubes
- the outer diameter of the tubes
- the thickness of the plain fins
- the fin density, i.e. the number of fins per unit length of tube
- the thermal conductivity of the fin material

Three regions can be distinguished on the refrigerant side of the condenser:

- the desuperheating region, where superheated refrigerant vapor from the compressor is cooled to saturated vapor

- the condensing region, where saturated vapor is condensing to saturated liquid

- the subcooling region, where liquid refrigerant is further cooled

The operation of the condenser depends on:

- the state of air at the inlet of the condenser
- the mass flow rate of air through the condenser
- the state of refrigerant at the inlet of the condenser
- the mass flow rate of refrigerant through the condenser

To fully determine the operation of the condenser, we start at the subcooling region of the condenser. The solving process is based on a root finding algorithm. The aim is to find the flow lengths of the subcooling, condensing, and desuperheating regions of the condenser, such that the sum of these flow lengths equals the total flow length of the condenser, i.e. the difference between the sum of these three flow lengths and the actual, total flow length of the condenser must be zero. The algorithm searches between a minimum and maximum value for the subcooling flow length. The maximum value can certainly not be larger than the actual flow length of the condenser; the minimum value is the smallest subcooling flow length for which the heat transfer characteristics of the subcooling region can be solved. When a solution is found for the subcooling region, the inlet and outlet states of air at the condensing and desuperheating regions can be determined. It is then, for the given mass flow rate of refrigerant and air, also possible to find the condensing flow length needed to condense the refrigerant from saturated vapor to saturated liquid and the desuperheating flow length needed to desuperheat the compressed refrigerant vapor to saturated vapor. The root finding algorithm will first check the difference at the minimum value and at the maximum value. If the solutions have a different sign, it means that somewhere between the minimum and maximum value for the subcooling flow length, a solution exists for which the difference between the sum of the three flow lengths and the actual, total flow length of the condenser will be equal to zero.

## Subcooling Region

The state of refrigerant at the inlet of the subcooling region is known. It is saturated liquid (vapor quality 0 %) of which the pressure is assumed to be equal to the pressure of the refrigerant at the inlet of the condenser. In other words, it is assumed that the refrigerant pressure in the condenser remains constant (i.e. the condenser pressure). Air enters the condenser at the subcooling region. The state of air at the condenser inlet is known. The mass flow rate of air and the mass flow rate of refrigerant are known. With these data and the current value for the subcooling flow length it is possible to solve for the heat transfer characteristics across the subcooling region. An iterative solution technique is used to accomplish this:

1. Reasonable initial values are guessed for the specific heats of air and refrigerant. With these and the known mass flow rates of air and refrigerant, the capacitance rates can be calculated.

2. Also, a reasonable initial value is guessed for the heat transfer effectiveness across the subcooling region.

3. Now the heat transfer rate across the subcooling region can be calculated and the outlet temperatures of air and refrigerant can be determined.

4. With the known inlet states of air and refrigerant, the calculated outlet temperatures and capacitance rates, the mean state of air and refrigerant across the subcooling region can be determined.

5. Based on the current value for the subcooling flow length, the mass flow rates and the mean states of air and refrigerant, the overall heat transfer conductance of the subcooling region can be determined.

6. Using the effectiveness-NTU method, the heat transfer characteristics of the subcooling region can be determined: a new value for the heat transfer effectiveness, the heat transfer rate and the outlet temperatures of air and refrigerant.
7. With the new heat transfer rate and outlet temperatures, the capacitance rates of air and refrigerant can be recalculated.
8. Steps 4 to 7 are repeated until the difference between the last and previous calculated value of the heat transfer effectiveness becomes small enough (the difference is equal or smaller than the accepted tolerance). 

## Condensing Region

The state of refrigerant at the inlet and at the outlet of the condensing region are already known from the start: saturated vapor and saturated liquid. Therefore, as the mass flow rate of refrigerant is given, the heat transfer rate across the condensing region can immediately be solved for. Using the mean specific heat of humid air, and as the mass flow rate of air is given, it is also possible to immediately calculate the temperature rise of air across the condensing region of the condenser.

After the subcooling region has been solved, we know the state of air entering the condensing region. As the temperature rise of air across the condensing region is already determined, and as the humidity ratio of air across the condenser remains constant and equal to the known humidity ratio of air at the condenser inlet, we can also determine the state of air leaving the condensing region and entering the desuperheating region. This allows to determine the mean state of air and refrigerant across the condensing region.

Now, the flow length of the condensing region can be determined such that for the given mass flow rate of refrigerant the refrigerant has been transformed from saturated vapor to saturated liquid. Again this is done with an iterative solution technique:

1. An initial value for the condensing flow length is guessed.
2. The internal (refrigerant-side) convection heat transfer coefficient and the mean wall temperature are determined.
3. With the mean refrigerant temperature, mean wall temperature, and the already known heat transfer rate, the thermal resistance on the refrigerant-side can be calculated (in other words: we search for the internal thermal resistance for which the heat transfer rate will equal the already known heat transfer rate).
4. Knowing the internal thermal resistance and convection heat transfer coefficient, the internal heat transfer surface area can be calculated that is required to attain the already known heat transfer rate.
5. From the internal heat transfer surface area, a new value of the condensing flow length can be deduced.
6. Steps 2 to 5 are repeated until the difference between the last and previous calculated value for the condensing flow length becomes sufficiently small (i.e. the difference is equal or smaller than the accepted tolerance).

## Desuperheating Region

The same procedure as for the condensing region can be followed. The state of refrigerant at the inlet and at the outlet of the desuperheating region are already known from the start. At the outlet the refrigerant must be saturated vapor (vapor quality 100 %). So, as in the case of the condensing region, the heat transfer rate and temperature rise of air across the desuperheating region can be immediately solved for. Consequently, the state of air leaving the condenser is now also determined, as the state of air entering the desuperheating region has already been determined when solving for the condensing region. Again, the mean state of air and refrigerant across the desuperheating region are determined and the flow length of the desuperheating region is calculated by iteration.

