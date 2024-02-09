# reaction rate parameters (added support for fate id as a species) for sensitivity analysis
using Base.Threads;
const SPLOCK = SpinLock(); # spin lock for multi-threading security

# Define structure to hold Cell parameters
#-------------------------------------------------------------------------
mutable struct Cell
    birthday::Float64             # hrs
    current_idx::Int64
    generation72::Float64
    generation60::Float64
    generation48::Float64
    generation36::Float64
    generation24::Float64
    death_time::Float64                     # (0: surviving, 1: death, 2: division, 3: differentiation)
    div0_time::Float64               # integrator.t at fate (hrs)
    WTcRel::Array{Float64}
    Np50p50::Array{Float64}
    cRelSS::Float64
    NcRel::Array{Float64}
    cMyc::Array{Float64}
    NcRelSS::Float64
    cMycSS::Float64
    IkBaC50::Float64
    IkBbC50::Float64
    IkBeC50::Float64
    IkBdC50::Float64
    IkBaC50n::Float64
    IkBbC50n::Float64
    IkBeC50n::Float64
    IkBdC50n::Float64
    cRel::Float64
    cReln::Float64
    C50::Float64
    C50n::Float64
    rates::Matrix{Float64}
    reactionFlux::Matrix{Float64}
    historicFlux::Array{Float64}
    finalConcentration::Array{Float64}
    daughter_1_idx::Int64           # for keeping track of daughter cells when doing depth-first computation
    daughter_2_idx::Int64
end

# Define function to set Cell parameters
#-------------------------------------------------------------------------
function setCell(; birthday = 0.0, current_idx = 0, generation72 = 0.0, generation60 = 0.0, generation48 = 0.0, generation36 = 0.0, generation24 = 0.0, death_time = 0.0, div0_time = 0.0, 
    WTcRel = zeros(Float64, 720), Np50p50 = zeros(Float64, 720), cRelSS = 0.0, NcRel = zeros(Float64, 720), cMyc = zeros(Float64, 720), 
    NcRelSS = 0.0,cMycSS = 0.0,IkBaC50 = 0.0,IkBbC50 = 0.0,IkBeC50 = 0.0,IkBdC50 = 0.0,
    IkBaC50n = 0.0,IkBbC50n = 0.0,IkBeC50n = 0.0,IkBdC50n = 0.0,cRel = 0.0,cReln = 0.0,C50 = 0.0,C50n = 0.0, 
    rates = Matrix{Float64}(undef, TOTAL_SPECIES, TOTAL_PARAMS), reactionFlux = zeros(Float64, TOTAL_SPECIES, TOTAL_PARAMS), historicFlux = zeros(Float64, TOTAL_SPECIES), finalConcentration = zeros(Float64, TOTAL_SPECIES), daughter_1_idx = 0, daughter_2_idx = 0)
    return Cell(birthday, current_idx, generation72, generation60, generation48, generation36, generation24, death_time, div0_time, WTcRel, 
    Np50p50, cRelSS, NcRel, cMyc, NcRelSS, cMycSS, IkBaC50, IkBbC50, IkBeC50, IkBdC50, IkBaC50n, IkBbC50n, IkBeC50n, IkBdC50n, 
    cRel, cReln, C50, C50n, rates, reactionFlux, historicFlux, finalConcentration, daughter_1_idx, daughter_2_idx);
end

# Define function to initialize founder cells
#-------------------------------------------------------------------------
function initializeFounderCells!(Srates, allCells)
    @inbounds Threads.@threads for i in 1:FOUNDER_CELL_NUM
        # allCells[i] = setCell(current_idx = i, rates = distribute_params(Srates, i)); # CHANGE to this when doing sensitivity analysis
        allCells[i] = setCell(current_idx = i, rates = distribute_params(Srates)); # CHANGE to this when doing distributed simulation
        # Burn-in period for distributed cells to reach steady-state before simulation
        Pre_simulate!(allCells, i);
        # allCells[i].rates[ANTIGEN, INITIALCONC] = ANTIGEN_DOSE;
    end
    nothing
end



# Define function to pre-simulate a founder cell to its steady-state
#-------------------------------------------------------------------------
function Pre_simulate!(allCells, i)
    Bcell = SteadyStateProblem(computeNetworkNettFluxes!, (@view allCells[i].rates[1:end, INITIALCONC]), p=(allCells[i].rates, allCells[i].reactionFlux));
    solution = solve(Bcell, DynamicSS(method, tspan=BURN_IN_PERIOD), abstol=abserr, reltol=relerr, maxiters=1e7);
    allCells[i].rates[1:end, INITIALCONC] = solution.u;
    # a different burn-in period for the proliferation module (before it reaches full steady state)
    Bcell2 = SteadyStateProblem(computeNetworkNettFluxes_AP1!, (@view allCells[i].rates[1:end, INITIALCONC]), p=(allCells[i].rates, allCells[i].reactionFlux));
    solution = solve(Bcell2, DynamicSS(method, tspan=BURN_IN_PERIOD_PROLIF), abstol=abserr, reltol=relerr, maxiters=1e7);
    allCells[i].rates[1:E2F, INITIALCONC] = solution.u[1:E2F]; # CHANGED to include only initial proliferation species
    # allCells[i].rates[1:MYC, INITIALCONC] = solution.u[1:MYC]; # CHANGED to include only initial proliferation species
    solution = nothing;
    GC.gc();
    nothing
end

# Define function to simulate a single cell
#-------------------------------------------------------------------------
function Simulate!(allCells, Srates, i, IKKCurve, delay)
    @inbounds begin
        has_children = false;
        # if CD40L_DELAY == 0
        #     allCells[i].rates[CD40L, INITIALCONC] = CD40L_DOSE;
        # end

        Bcell = DDEProblem(computeNetworkNettFluxes!, allCells[i].rates[1:end, INITIALCONC], delay, (0.0, GLOBAL_END_TIME-allCells[i].birthday), (allCells[i].birthday, IKKCurve, allCells[i].rates, allCells[i].reactionFlux, allCells[i].historicFlux); constant_lags=[0.25, 0.75, 1.0, 4.0, 12.0]);
        NFkBsolution = solve(Bcell, MethodOfSteps(method), callback=cellFate, abstol=abserr, reltol=relerr, save_everystep=false, saveat = 0.1, maxiters=1e10);
        Bcell2 = ODEProblem(computeNetworkNettFluxes_AP2!, allCells[i].rates[1:end, INITIALCONC], (0.0, GLOBAL_END_TIME-allCells[i].birthday), (allCells[i].birthday, allCells[i].rates, allCells[i].reactionFlux, NFkBsolution));
        solution = solve(Bcell2, method, callback=cellFate, abstol=abserr, reltol=relerr, save_everystep=false, saveat = 0.1, maxiters=1e10);
        allCells[i].WTcRel = NFkBsolution[NC50, :] + NFkBsolution[NCREL, :] + NFkBsolution[C50, :] + NFkBsolution[CREL, :] + 
        NFkBsolution[IKBAC50, :] + NFkBsolution[IKBBC50, :] + NFkBsolution[IKBEC50, :] + NFkBsolution[IKBDC50, :] +
        NFkBsolution[NIKBAC50, :] + NFkBsolution[NIKBBC50, :] + NFkBsolution[NIKBEC50, :] + NFkBsolution[NIKBDC50, :];
        allCells[i].Np50p50 = NFkBsolution[NP50P50, :]
        allCells[i].NcRel = NFkBsolution[NC50, :] + NFkBsolution[NCREL, :] + NFkBsolution[NIKBAC50, :] + NFkBsolution[NIKBBC50, :] + NFkBsolution[NIKBEC50, :] + NFkBsolution[NIKBDC50, :];

        allCells[i].cMyc = solution[MYC, :];
        allCells[i].cRelSS = NFkBsolution[NC50, 1] + NFkBsolution[NCREL, 1] + NFkBsolution[C50, 1] + NFkBsolution[CREL, 1] + 
        NFkBsolution[IKBAC50, 1] + NFkBsolution[IKBBC50, 1] + NFkBsolution[IKBEC50, 1] + NFkBsolution[IKBDC50, 1] +
        NFkBsolution[NIKBAC50, 1] + NFkBsolution[NIKBBC50, 1] + NFkBsolution[NIKBEC50, 1] + NFkBsolution[NIKBDC50, 1]; # whole-cell cRel at Steady State
        allCells[i].NcRelSS = (allCells[i].NcRel)[1]; # nuclear cRel at Steady State
        allCells[i].cMycSS = (allCells[i].cMyc)[1]; # cMyc at Steady State
        allCells[i].IkBaC50 = NFkBsolution[IKBAC50, 1]; # IkBa-bound cRel at Steady State
        allCells[i].IkBbC50 = NFkBsolution[IKBBC50, 1]; # IkBb-bound cRel at Steady State
        allCells[i].IkBeC50 = NFkBsolution[IKBEC50, 1]; # IkBe-bound cRel at Steady State
        allCells[i].IkBdC50 = NFkBsolution[IKBDC50, 1]; # IkBd-bound cRel at Steady State
        allCells[i].IkBaC50n = NFkBsolution[NIKBAC50, 1]; # nuclear IkBa-bound cRel at Steady State
        allCells[i].IkBbC50n = NFkBsolution[NIKBBC50, 1]; # nuclearIkBb-bound cRel at Steady State
        allCells[i].IkBeC50n = NFkBsolution[NIKBEC50, 1]; # nuclearIkBe-bound cRel at Steady State
        allCells[i].IkBdC50n = NFkBsolution[NIKBDC50, 1]; # nuclearIkBd-bound cRel at Steady State
        allCells[i].cRel = NFkBsolution[CREL, 1];
        allCells[i].cReln = NFkBsolution[NCREL, 1];
        allCells[i].C50 = NFkBsolution[C50, 1];
        allCells[i].C50n = NFkBsolution[NC50, 1];
        # allCells[i].cMyc24 = (length(allCells[i].cMyc)>=240) ? (allCells[i].cMyc)[240] : 0.0; # cMyc at 24hrs

        allCells[i].finalConcentration[1:NFKB_SPECIES+DIFF_SPECIES] = NFkBsolution(solution.t[end])[1:NFKB_SPECIES+DIFF_SPECIES];
        allCells[i].finalConcentration[NFKB_SPECIES+DIFF_SPECIES+1:TOTAL_SPECIES] = solution.u[end][NFKB_SPECIES+DIFF_SPECIES+1:TOTAL_SPECIES];
        # allCells[i].generation96 = allCells[i].finalConcentration[GEN] # ignoring cell death
        # allCells[i].generation2 = allCells[i].finalConcentration[TOTAL_SPECIES] # considering cell death # CHANGED removed
        allCells[i].generation72 = allCells[i].finalConcentration[GEN] # ignoring cell death
        allCells[i].generation60 = allCells[i].finalConcentration[GEN60] # ignoring cell death
        allCells[i].generation48 = allCells[i].finalConcentration[GEN48] # ignoring cell death
        allCells[i].generation36 = allCells[i].finalConcentration[GEN36] # ignoring cell death
        allCells[i].generation24 = allCells[i].finalConcentration[GEN24] # ignoring cell death

        allCells[i].death_time = (solution.u[end][DEATH_TIME] != 0) ? solution.u[end][DEATH_TIME] : 96.0;   # 2nd division
        allCells[i].div0_time = (solution.u[end][DIV0_TIME] != 0) ? solution.u[end][DIV0_TIME] : 96.0;   # 1st division 
        # value1 = (i-1) % ST_num + 1;
        # value2 = (i-1) ÷ ST_num + 1;
        # print(string(allCells[i].current_idx), '\t', string(value1), '\t', string(value2), '\t', string(allCells[i].div0_time), '\t', string(allCells[i].death_time), '\n');

        NFkBsolution = nothing;
        solution = nothing;
        GC.gc();
    end
    # print(has_children, " ", allCells[i].daughter_1_idx, " ", allCells[i].daughter_2_idx)
    return nothing
end

# Define function to simulate a cell lineage (depth-first)
#-------------------------------------------------------------------------
# function Simulate_lineage!(allCells, Srates, i, delay)
#     stack = Int[i]
#     while !isempty(stack)
#         j = pop!(stack)
#         has_children, daughter_1_idx, daughter_2_idx = Simulate!(allCells, Srates, j, delay)
#         if has_children
#             push!(stack, daughter_1_idx)
#             push!(stack, daughter_2_idx)
#         end
#         lock(SPLOCK);

#         cellsSaver = jldopen(cells_fn, "r+");
#         write(cellsSaver, string(j), allCells[j]);
#         close(cellsSaver);

#         print(allCells[j].birthday, '\t', allCells[j].current_idx, '\t', allCells[j].parent_idx, '\t', allCells[j].generation, '\t', allCells[j].fate, '\t', allCells[j].fate_t, '\t', allCells[j].abs_fate_t, '\t', allCells[j].daughter_1_idx, '\t', allCells[j].daughter_2_idx, '\n');
#         unlock(SPLOCK);
#         # GC.gc(true)
#         # ccall(:malloc_trim, Cvoid, (Cint,), 0)
#     end
# end

# Define function to simulate all cell lineages serially
#-------------------------------------------------------------------------
function Simulate_nonthreaded!(allCells, Srates, IKKCurve, delay)
    @inbounds for i in 1:FOUNDER_CELL_NUM
        Simulate!(allCells, Srates, i, IKKCurve, delay)
    end
end

# Define function to simulate all cell lineages w/ multi-threading (static schedule)
#-------------------------------------------------------------------------
function Simulate_threaded!(allCells, Srates, IKKCurve, delay)
    Threads.@threads for i in 1:FOUNDER_CELL_NUM
        Simulate!(allCells, Srates, i, IKKCurve, delay)
    end
end

# Define function to simulate all cell lineages w/ multi-threading (dynamic schedule)
#-------------------------------------------------------------------------
function Simulate_spawned!(allCells, Srates, IKKCurve, delay)
    tasks = [Threads.@spawn(Simulate!(allCells, Srates, i, IKKCurve, delay)) for i in 1:FOUNDER_CELL_NUM]
    out = [fetch(t) for t in tasks]
end
