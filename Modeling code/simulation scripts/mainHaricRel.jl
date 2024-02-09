# Hari's cRel abundance paper: Test the correlation between basal and peak cRel and Tdiv0 and Ndiv
# Using old version of NFkB and proliferation modules
# Include libraries / packages required
#--------------------------------------------
using DifferentialEquations;

using ArgParse; # Argument parsing

# Benchmarking & profiling packages
using BenchmarkTools;
using StatProfilerHTML;

# Visualization packages
using DataFrames;
using Plots;
using Gadfly;
using Cairo;
using JLD;

# Include source files
#--------------------------------------------
# include("ConstantParams.jl");
include("ReactionRates3_old2.jl");
include("HelperFunctions_old.jl");

# include("ODE_Receptor3.jl");
include("ODE_NFkB2.jl");
include("ODE_Apoptosis2_old.jl");
include("ODE_Differentiation.jl");
include("ODE_Proliferation2_old.jl"); # CHANGED this to ODE_Proliferation_old.jl from ODE_Proliferation2_old.jl

include("SimulateFunctionsHaricRel2.jl");

# Argument parsing w/ command line input
#--------------------------------------------
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--version", "-v"
            help = "nonthread, thread (static schedule), or spawn (dynamic schedule) version of lineage simulation"
            arg_type = String
            default = "nonthread"
        "--initial", "-i"
            help = "reaction rates and steady states for each starting cell after initialization (.jld format)"
            arg_type = String
            default = "initial.jld"
        "--cells", "-c"
            help = "simulated outputs from each cell (.jld format)"
            arg_type = String
            default = "cells.jld"
        "--output", "-o"
            help = "output file name for simulated cell lineages"
            arg_type = String
            default = "output.txt"
        "--param", "-p"
            help = "whether to load a parameter set from a .jl file"
            arg_type = String
            default = ""
            required = false
        "--reload", "-r"
            help = "whether to reload from previous steady states, provided by the --initial command"
            arg_type = Bool
            action = :store_true
            default = false
            required = false
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
const version = get(parsed_args, "version", "nonthread");
const steady_fn = get(parsed_args, "initial", "initial.jld");
const cells_fn = get(parsed_args, "cells", "cells.jld");
const output_fn = get(parsed_args, "output", "output.txt");
const param_fn = get(parsed_args, "param", "");
const reload = get(parsed_args, "reload", false);

if (param_fn != "")
    include(param_fn);
    print("Parameter set loaded!");
end

# const version = "spawn";
# const steady_fn = "/Users/helenhuang/Desktop/test_steady.jld";
# const cells_fn = "/Users/helenhuang/Desktop/test_cells.jld";
# const output_fn = "/Users/helenhuang/Desktop/test_output.txt";

# Define standard reaction rates
#--------------------------------------------
# rates = Matrix{Float64}(undef, TOTAL_SPECIES, TOTAL_PARAMS);
rates = Matrix{Float64}(undef, TOTAL_SPECIES, TOTAL_PARAMS); # CHANGED
setAllRates!(rates);
const Srates = rates;

# Define ODEs for network
#--------------------------------------------
# Attempt to speed up the solver by splitting into 2 parts
# time-independent ODE (pre-simulation: phase = 1, only receptor & NFkB needs steady state simulation)
function computeNetworkNettFluxes!(nettFlux, concentration, (Srates, reactionFlux), time)
    # computeReceptorNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    # computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1); # CHANGED to prestimulation to vary cMyc SS level
    nothing
end

function computeNetworkNettFluxes_AP1!(nettFlux, concentration, (Srates, reactionFlux), time)
    # computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 1);
    nothing
end

# time-dependent ODE (simulation, w/ delay: phase = 2)
function computeNetworkNettFluxes!(nettFlux, concentration, delay, (birthday, IKKCurve, Srates, reactionFlux, historicFlux), time)
    # computeReceptorNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, time);
    computeNFkBNettFluxes!(nettFlux, concentration, delay, reactionFlux, Srates, 2, time, birthday, IKKCurve, historicFlux)
    # computeNFkBNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, delay, historicFlux, time);
    # computeDiffNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    nothing
end

function computeNetworkNettFluxes_AP2!(nettFlux, concentration, (birthday, Srates, reactionFlux, inputCurves), time)
    computeApoptosisNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2, time, birthday, inputCurves);
    computeProlifNettFluxes!(nettFlux, concentration, reactionFlux, Srates, 2);
    nothing
end


# Callback function to detect cell death, mitotic and differentiation events
#----------------------------------------------------------------------
function condition(out, u, t, integrator)
    out[1] = u[CPARP] - 2500;
    out[2] = u[CDH1] - 0.2;
    # out[3] = u[IRF4] - 0.65*u[BCL6] - 1.2;
    # out[4] = t - CD40L_DELAY;
end

function affect!(integrator, index)
    if (index == 1)
        if (integrator.u[DEATH_TIME] == 0) # death time
            integrator.u[DEATH_TIME] = integrator.t;
        end
    end
    if (index == 2)
        if (integrator.u[CYCB] > 2.0)
            integrator.u[MASS] /= 2.0;
            integrator.u[GEN] += 1.0; # generation (ignoring cell death)
            # if (integrator.u[DEATH_TIME] == 0) # cell is not dead
            #     integrator.u[TOTAL_SPECIES] += 1.0; # generation (considering cell death)
            # end # CHANGED removed
            if (integrator.u[DIV0_TIME] == 0) # 1st division
                integrator.u[DIV0_TIME] = integrator.t;
            end
            if (integrator.t <= 24) # CHANGED added
                integrator.u[GEN60] += 1.0;
                integrator.u[GEN48] += 1.0;
                integrator.u[GEN36] += 1.0;
                integrator.u[GEN24] += 1.0;
            elseif (integrator.t <= 36) # CHANGED added
                integrator.u[GEN60] += 1.0;
                integrator.u[GEN48] += 1.0;
                integrator.u[GEN36] += 1.0;
            elseif (integrator.t <= 48)
                integrator.u[GEN60] += 1.0;
                integrator.u[GEN48] += 1.0;
            elseif (integrator.t <= 60)
                integrator.u[GEN60] += 1.0;
            end
        end
    # elseif (index == 3)
    #     integrator.u[TOTAL_SPECIES] = 3;
    #     terminate!(integrator);
    # elseif (index == 4)
    #     integrator.u[CD40L] = CD40L_DOSE;
    #     integrator.u[ANTIGEN] = 0;
        # integrator.u[ABCR] = 0;
    end
    nothing
end

cellFate = VectorContinuousCallback(condition, affect!, nothing, 2, save_positions=(false, false));

# Set up a structure to hold cells
#-------------------------------------------------------------------------
allCells = Vector{Cell}(undef, FOUNDER_CELL_NUM);
if reload
    allCells = load(steady_fn, "allCells");
    print("Old steady states loaded!")
else
    initializeFounderCells!(Srates, allCells);
    JLD.save(steady_fn, "allCells", allCells);
end

# Define default delay functions (t < tau)
#--------------------------------------------
# const delay(historicFlux, p, t) = (historicFlux .= 0.0);
const delay(p, t; idxs=nothing) = typeof(idxs) <: Number ? 0.0 : zeros(TOTAL_SPECIES); # CHANGED

# Define input IKK curve (according to stimulus type)
#--------------------------------------------
const IKKCurve = specifyIKKCurve(IKK_shape = IKK_TYPE);


# print("current_idx", '\t', "div0_time", '\t', "div1_time", '\n');

# Simulate cell lineages
#--------------------------------------------
if version == "nonthread"
    @time Simulate_nonthreaded!(allCells, Srates, IKKCurve, delay);
elseif version == "thread"
    @time Simulate_threaded!(allCells, Srates, IKKCurve, delay);
elseif version == "spawn"
    @time Simulate_spawned!(allCells, Srates, IKKCurve, delay);
else
    print("Please input the correct version -v: nonthread, thread, or spawn.");
end

# Output information about all cells (for visualization)
#--------------------------------------------
JLD.save(cells_fn, "allCells", allCells);

out = open(output_fn, "w");
write(out, "current_idx", '\t', "div0_time", '\t', "death_time", '\t', "WCcRelSS", '\t', "NcRelSS", '\t', "cMycSS", '\t', 
    # "IkBaC50", '\t', "IkBbC50", '\t', "IkBeC50", '\t', "IkBdC50", '\t', "IkBaC50n", '\t', 
    # "IkBbC50n", '\t', "IkBeC50n", '\t', "IkBdC50n", '\t', "cRel", '\t', "cReln", '\t', "C50", '\t', "C50n", '\t', 
    "gen72", '\t', "gen60", '\t', "gen48", '\t', "gen36", '\t', "gen24", '\n');
for i in 1:length(allCells)
    write(out, string(allCells[i].current_idx), '\t', string(allCells[i].div0_time), '\t', string(allCells[i].death_time), '\t', string(allCells[i].cRelSS), '\t', 
    string(allCells[i].NcRelSS), '\t', string(allCells[i].cMycSS), '\t', #string(allCells[i].IkBaC50), '\t', string(allCells[i].IkBbC50), 
    # '\t', string(allCells[i].IkBeC50), '\t', string(allCells[i].IkBdC50), '\t', string(allCells[i].IkBaC50n), '\t', 
    # string(allCells[i].IkBbC50n), '\t', string(allCells[i].IkBeC50n), '\t', string(allCells[i].IkBdC50n), '\t', 
    # string(allCells[i].cRel), '\t', string(allCells[i].cReln), '\t', string(allCells[i].C50), '\t', string(allCells[i].C50n), '\t', 
    string(allCells[i].generation72), '\t', string(allCells[i].generation60), '\t', string(allCells[i].generation48), '\t', string(allCells[i].generation36), '\t', string(allCells[i].generation24), '\n');
end
close(out);
