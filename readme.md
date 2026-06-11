# Read me for accompanying code to  “Hypothalamic recurrent inhibition regulates functional states of stress effector neurons” 

Wilten Nicola, June 11th, 2026. 

The accompanying code consists of 8 scripts, consisting of 4 simulation loops and 4 plotting scripts.  The four simulation loops are:

**TRACKING_MODEL_1**:  Runs a simulation and saves the results of the network in Figure 1O from [1].  

**TRACKING_MODEL_2**: Runs a simulation and saves the results of the network in Figure 1P from [1].  

**DISINHIBITION_MODEL_1**:  Runs a simulation and saves the results of the network in Figure 1Q from [1].  

**DISINHIBTION_MODEL_2**:  Runs a simulation and saves the results of the network in Figure 1R from [1].  

The simulation scripts have to be run before the plotting scripts are run, and approximately 20GB of storage is required to store the results of the 4 simulations (in total).  Each simulation should take ~5 minutes at the most, depending on computing resources.    With the simulations concluded, the appropriate plotting script for the models above can be run to visualize the results, similarly to the sub-panels from Figure 1 O-R.   Note that while the RNG seed is fixed, the results may differ slightly from the original run and figures from [1] due to differences in MATLAB iterations over the years.  

[1] Ichiyama, A., Mestern, S., Fuzesi, T., Allman, B.L., Nicola, W., Bains, J., Muller, L. and Inoue, W., 2025. Hypothalamic recurrent inhibition regulates functional states of stress effector neurons. *bioRxiv*, pp.2025-09.
