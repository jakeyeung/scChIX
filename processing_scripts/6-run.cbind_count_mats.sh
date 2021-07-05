#!/bin/sh
# Jake Yeung
# 6-run.cbind_count_mats.sh
#  
# 2021-06-27

jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/cbind_count_mats.R"
# jmark="K9m3"
# jmark="K36"
jmark="K36-K9m3"
# jmark="K27"
# jmark="K36-K27"
# jmark="K4m1"
# jmark="K9m3"
# jmark="K36-K4m1"

jsuffix="10000"

indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/counts_tables_${jsuffix}/${jmark}"

outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/counts_tables_${jsuffix}/${jmark}/cbind_out"

[[ ! -d $outdir ]] && mkdir $outdir

outf=${outdir}/${jmark}.countTable.binsize_${jsuffix}.rds

BNAME=${outdir}/${jmark}.sbatch_log
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

infs=${indir}/*.csv

# echo $infs

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $infs -outfile $outf"
echo $cmd
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jmark} --wrap "$cmd"
