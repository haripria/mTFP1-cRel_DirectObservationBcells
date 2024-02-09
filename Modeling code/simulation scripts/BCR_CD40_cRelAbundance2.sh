source /home/helen/.bashrc
export JULIA_NUM_THREADS=50

home_dir="/home/helen/BCELL_PROJECT/"

for i in {6..10}
do
    # With competitive inhibition
    modifier="cRelAbundance_KOF${i}"
    julia $home_dir'scripts/mainHaricRel.jl' -v "spawn" -o $home_dir'results/'$modifier'.txt' -i $home_dir'data/steady_'$modifier'.jld' -c $home_dir'data/cells_'$modifier'.jld' -p $home_dir'params/Param1_cRel_IKBEKO.jl' >> $home_dir'job-logs/'$modifier'.out'

    modifier="cRelAbundance_WTF${i}"
    julia $home_dir'scripts/mainHaricRel.jl' -v "spawn" -o $home_dir'results/'$modifier'.txt' -i $home_dir'data/steady_'$modifier'.jld' -c $home_dir'data/cells_'$modifier'.jld' -p $home_dir'params/Param1_cRel_WT.jl' >> $home_dir'job-logs/'$modifier'.out'

    # No competitive inhibition
    modifier="cRelAbundance_WTF${i}n"
    julia $home_dir'scripts/mainHaricRel.jl' -v "spawn" -o $home_dir'results/'$modifier'.txt' -i $home_dir'data/steady_'$modifier'.jld' -c $home_dir'data/cells_'$modifier'.jld' -p $home_dir'params/Param2_cRel_WT.jl' >> $home_dir'job-logs/'$modifier'.out'

    modifier="cRelAbundance_KOF${i}n"
    julia $home_dir'scripts/mainHaricRel.jl' -v "spawn" -o $home_dir'results/'$modifier'.txt' -i $home_dir'data/steady_'$modifier'.jld' -c $home_dir'data/cells_'$modifier'.jld' -p $home_dir'params/Param2_cRel_IKBEKO.jl' >> $home_dir'job-logs/'$modifier'.out'

done
=