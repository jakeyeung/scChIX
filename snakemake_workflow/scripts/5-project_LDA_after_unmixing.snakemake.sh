#!/bin/sh
# Jake Yeung
# 1-project_new_samples_on_LDA.H3K27me3.StemCells.sh
#  
# 2019-10-24

inlda=$1
inmat=$2
outf=$3

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

[[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
[[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $inmat $outf --RemoveEmptyCells

