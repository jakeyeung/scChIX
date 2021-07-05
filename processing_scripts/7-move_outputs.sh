#!/bin/sh
# Jake Yeung
# 7-move_outputs.sh
#  
# 2021-06-28

indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/RZ_counts"
cd $indir

dnames="K36-K27m3 K36-K9m3 K36-K4m1 K9m3 K4m1 K36m3"

for d in $dnames; do
    mkdir $d
    mv *${d}*.csv $d
done
