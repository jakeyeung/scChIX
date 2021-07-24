#!/bin/sh
# Jake Yeung
# 7-move_outputs.sh
#  
# 2021-06-28

# indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/RZ_counts"
jsuffix="counts_tables_topfeatures_K36_genebodies_K9m3_bins"
indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/tagged_bams/${jsuffix}"
outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/${jsuffix}"
[[ ! -d $outdir ]] && mkdir $outdir
cd $outdir

dnames="K36-K27 K36-K9m3 K36-K4m1 K9m3 K4m1 K36 K27"

for d in $dnames; do
    mkdir $d
    for f in `ls -d $indir/*-${d}*.txt`; do 
      # ln -s $f $d  # does not work you will double move doubleincubated
      mv $f $d
    done
done
