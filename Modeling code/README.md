To generate figures for the mTFP1-cRel paper, one could use the simulated data that we generated for the paper (in the `simulated data` dir), 
or generate new data from scratch using the scripts in the `simulation scripts` dir.
Notice that the simulation script `mainHaricRel.jl` use parameters in the `params` dir.

### Simulation Scripts:
- `mainHaricRel.jl`: function for running the simulations and saving results
- `ConstantParams.jl`: constant parameters, including stimulus doses, stimulus delay, simulation time, scaling factors, etc.
- `ReactionRates3_old.jl`: reaction rate parameters for all module
- `ODE_NFkB2.jl`: ODE equations for NFkB module (before integrating noncanonical)
- `ODE_Apoptosis2_old.jl`: ODE equations for Apoptosis module (removed the presence of A52, C52, and ABCR)
- `ODE_Differentiation.jl`: ODE equations for Differentiation module
- `ODE_Proliferation2_old.jl`: ODE equations for Cell Cycle module (modifed from Mitchell et al. 2018 PNAS)
- `SimulateFunctionsHaricRel.jl`: pre-simulation and simulation functions
- `HelperFunctions_old.jl`: helper functions for Michaelis-Menten and Hill functions, as well as parameter distributions

### Parameter sets:
Only includes parameters that are modified from the constant parameter file (`ConstantParams.jl` that is shared among all conditions):
- `Param1_cRel_IKBEKO.jl`: IkBe-KO w/ IFFL
- `Param1_cRel_WT.jl`: WT w/ IFFL
- `Param2_cRel_IKBEKO.jl`: IkBe-KO in direct condition (no IFFL)
- `Param2_cRel_WT.jl`: WT in direct condition (no IFFL)

### Plotting scripts:
- `NFkB trajectories (cRel paper).ipynb`: plot the trajectory of nuclear p50:p50, nuclear cRel, and whole-cell cMyc (Figure 7I) between WT and IkBe-KO
- `cRel WT vs EKO gen histogram.ipynb`: plot the generation histogram and ratio (Figure 7J) between WT and IkBe-KO
