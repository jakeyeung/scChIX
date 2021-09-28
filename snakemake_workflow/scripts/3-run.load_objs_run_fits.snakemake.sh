#!/bin/sh
# Jake Yeung
# 2b-run.load_objs_run_fits.from_server.sh
#  
# 2019-08-08

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/load_objs_run_fits.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

indat=$1
outdat=$2
ncores=$3
jmethod="Brent"

[[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
# [[ -e $outdat ]] && echo "$outdat found, exiting" && exit 1
[[ $ncores != [0-9]* ]] && echo "Must be integer: $ncores" && exit 1

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $indat $outdat --ncores $ncores --method $jmethod
