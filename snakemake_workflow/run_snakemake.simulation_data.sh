#!/bin/sh
# Jake Yeung
# run_snakemake.sh
#  
# 2021-09-10


indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/simulation_data"
cd $indir

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; snakemake --cluster-config cluster.json --cluster "sbatch --time={cluster.time} --mem-per-cpu={cluster.mem} --cpus-per-task={cluster.N}" --jobs 50 --latency-wait 60

