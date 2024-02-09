To generate figures for the mTFP1-cRel paper, one could use the simulated data that we generated for the paper (Link to simulated data?), 
or generate new data from scratch using the scripts in the `simulation scripts` dir.
Notice that the simulation script `mainHaricRel.jl` use parameters in the `params` dir. To be able to run the scripts, you need to have installed julia (> v1.7) and the required packages each script specifies.
- To plot the figures from our simulated data (Link to simulated data?), one can modify the `cRel WT vs EKO gen histogram.ipynb` jupyter notebook to include the correct working directory (`homedir` and `outdir` on line 1 and 2, block 7), then run all the blocks in the notebook to generate Figure 7J. The same applies to the `NFkB trajectories (cRel paper).ipynb` notebook for Figure 7I, where one needs to modify line 1 and 2 in block 5. Note that `homedir` is the directory where the simulated data is stored, and outdir is the directory where you would like the figures to be outputted.

- To run the simulation from scratch, one can modify the `BCR_CD40_cRelAbundance.sh` shell script to include the correct working directory, then run the shell script:
```
export JULIA_NUM_THREADS=50     # modify this line to the number of threads available on your computer or server
home_dir="/home/helen/BCELL_PROJECT/"     # modify this line to your working directory, that should include a simulation_scripts/, data/, and params/ subdirectories

# To run the simulation, type the following in a bash environment (i.e. Command Line on Mac):
./BCR_CD40_cRelAbundance.sh
```

### Simulation Scripts:
- `mainHaricRel.jl`: function for running the simulations and saving results
- `ConstantParams.jl`: constant parameters, including stimulus doses, stimulus delay, simulation time, scaling factors, etc.
- `ReactionRates3_old.jl`: reaction rate parameters for all module
- `ODE_NFkB2.jl`: ODE equations for NFkB module (before integrating noncanonical)
- `ODE_Apoptosis2_old.jl`: ODE equations for Apoptosis module (removed the presence of A52, C52, and ABCR)
- `ODE_Differentiation.jl`: ODE equations for Differentiation module (Roy et al. 2019 Immunity, included but not used)
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
