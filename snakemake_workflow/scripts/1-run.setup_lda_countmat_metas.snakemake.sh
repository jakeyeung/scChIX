#!/bin/sh
# Jake Yeung
# 0-setup_lda_countmat_metas.sh
#  
# 2021-07-16


jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/cluster_make_lda_countmat_metadata_objects.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1


mark1=$1
mark2=$2
markdbl=$3

inlda1=$4
inlda2=$5
inldadbl=$6

outdir=$7

nneighbors=$8

if [ -z ${nneighbors} ]; then nneighbors="30"; else echo "nneighbors is manually set"; fi

echo "Nneighbors set to: $nneighbors"

[[ ! -e $inlda1 ]] && echo "$inlda1 not found, exiting" && exit 1
[[ ! -e $inlda2 ]] && echo "$inlda2 not found, exiting" && exit 1
[[ ! -e $inldadbl ]] && echo "$inldadbl not found, exiting" && exit 1


marks="$mark1 $mark2 $markdbl"
inldas="$inlda1 $inlda2 $inldadbl"

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -marks $marks -infiles $inldas -outdir ${outdir} -n_neighbors ${nneighbors}

# BNAME=${outdir}/setup_sbatch_output
# DBASE=$(dirname "${BNAME}")
# [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
# cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -marks $marks -infiles $inldas -outdir ${outdir}"
# sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jstr}_${midname}_${jquant} --wrap "$cmd"
# 
