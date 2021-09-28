#!/bin/sh
# Jake Yeung
# 2-run_LDA.sh
#  
# 2019-12-06

# sleep 600

# # WRAP UP
# while [[ `squeue -u jyeung | grep unmix | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='16G'
jtime='24:00:00'


workdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline_unmixed_singles_LDA_together"
[[ ! -d $workdir ]] && mkdir $workdir
[[ ! -d $workdir ]] && echo "$workdir not found, exiting" && exit 1

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

jname="var_filtered_manual2nocenterfilt2_K36_K9m3_K36-K9m3"
inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/scchix_unmixing_downstream"
indir="${inmain}/${jname}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outmain="${workdir}/${jname}"
[[ ! -d $outmain ]] && mkdir $outmain
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

for inf in `ls -d $indir/*merged_mat*.rds`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}"
    [[ ! -d $outdir ]] && mkdir -p $outdir

    echo $bname

    BNAME=$outmain/$bname.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd="cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname --RemoveEmptyCells"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${bname} --wrap "$cmd"
done
wait
