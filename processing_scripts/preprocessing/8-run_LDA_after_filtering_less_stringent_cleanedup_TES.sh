#!/bin/sh
# Jake Yeung
# 8-run_LDA_after_filtering.sh
#  
# 2021-06-28
# Filtered using interactive scripts: eg File: ~/projects/scChIX/analysis_scripts/1-explore_raw_count_data_K27me3.R

jmem='16G'
jtime='96:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# filtdir="filtered_count_tables_filtered_counts_varcutoffmin_0.03"
# jmarkstr="K36_K9m3_K36-K9m3"
# jmarkstr="K36_K27_K36-K27"
# jmarkstr="K36_K4m1_K36-K4m1"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables/for_LDA"
# indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/${filtdir}/${jmarkstr}"

# # jsuffix="dbl_k36_k9m3_cleaned"
# jsuffix="dbl_cleaned"
jsuffix="TES_k36_cleaned"
indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_analysis/rds_objs_celltyping_NN_15_check_plates/${jsuffix}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

prefix="50000_${jsuffix}"
outmain0="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_outputs"
[[ ! -d $outmain0 ]] && mkdir $outmain0

outmain="${outmain0}/ldaAnalysis_${prefix}"
[[ ! -d $outmain ]] && mkdir $outmain
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

for inf in `ls -d $indir/countmat*.rds`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}"
    [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd="cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=LDA_${bname} --wrap "$cmd"
done
