#!/bin/sh
# Jake Yeung
# 0-run_LDA_initial.sh
#  
# 2021-06-28
# Filtered using interactive scripts: eg File: ~/projects/scChIX/analysis_scripts/1-explore_raw_count_data_K27me3.R


topics=30
rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/run_LDA_model3.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inf=$1
outf=$2

/hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf --topics $topics
