#!/bin/sh
# Jake Yeung
# run_snakemake.sh
#  
# 2021-09-10


indir="/home/hub_oudenaarden/jyeung/projects/snakemake_workflow_gastru_scchix"
cd $indir

# sbatch --time=  --mem-per-cpu=<+jmem+> --output=<+BNAME+>_%j.log --ntasks=1 --nodes=1 --cpus-per-task=<+NCORES+> --job-name=<+jobname+> --wrap "$cmd"
# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; snakemake --cluster-config cluster.json --cluster "sbatch --time={cluster.time} --mem-per-cpu={cluster.mem} --cpus-per-task={cluster.N} --output={cluster.log}" --jobs 50
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; snakemake --cluster-config cluster.json --cluster "sbatch --time={cluster.time} --mem-per-cpu={cluster.mem} --cpus-per-task={cluster.N}" --jobs 50 

