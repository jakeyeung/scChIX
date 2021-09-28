#!/bin/sh
# Jake Yeung
# 2b-run.setup_objs_for_unmixing.sh
# Visually look for interesting topics, create umap_annotated... RData objects which 
# are input to setting up for unmixing
# 2020-02-05


jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/setup_objs_for_unmixing_general_from_LDA_or_countmat.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1


mark1=$1
mark2=$2
inf1meta=$3
inf2meta=$4
infdblmeta=$5
inf1lda=$6
inf2lda=$7
infdblcountmat=$8
outprefix=$9

[[ ! -e $inf1meta ]] && echo "$inf1meta not found, exiting" && exit 1
[[ ! -e $inf2meta ]] && echo "$inf2meta not found, exiting" && exit 1
[[ ! -e $infdblmeta ]] && echo "$infdblmeta not found, exiting" && exit 1

[[ ! -e $infdblcountmat ]] && echo "$infdblcountmat not found, exiting" && exit 1

[[ ! -e $inf1lda ]] && echo "$inf1lda not found, exiting" && exit 1
[[ ! -e $inf2lda ]] && echo "$inf2lda not found, exiting" && exit 1


. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -mark1 $mark1 -mark2 $mark2 -inf1meta $inf1meta -inf2meta $inf2meta -infdblmeta $infdblmeta -inf1lda $inf1lda -inf2lda $inf2lda -infdblcountmat $infdblcountmat -outprefix $outprefix

