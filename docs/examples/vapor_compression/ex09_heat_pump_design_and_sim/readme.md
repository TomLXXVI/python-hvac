# Air-To-Air Heat Pump Design and Simulation

## Jupyter Notebook `heat_pump_design.ipynb`

Illustrates the design of an air-to-air heat pump model with python-hvac.
Uses the classes defined in Python-file `design_parameters.py` to determine
design parameters for the condenser and the evaporator (air mass flow rates,
evaporating temperature, condensing temperature, face area).

## Python Script `heat_pump_simulation.py`

Illustrates the simulation of the air-to-air heat pump model. Runs the 
simulations with a number of different outdoor air temperatures. Simulations
are executed in parallel using a `ProcessPoolExecutor`. Simulation results  
(`Output` objects) are stored on disk using `shelve` (files `heat_pump_rating.dat`,
`heat_pump_rating.dir` and `heat_pump_rating.bak`).

## Jupyter Notebook `heat_pump_simulation.ipynb`

Loads the shelved simulation results from disk and analyzes them: heating
capacity of the heat pump at different outdoor air temperatures, compressor
power, COP.
