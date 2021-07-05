#!/bin/sh
# Jake Yeung
# 7-move_outputs.sh
#  
# 2021-06-28

# indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/RZ_counts"
indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/tagged_bams/counts_tables_TSS10kb"
outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/counts_tables_TSS10kb"
[[ ! -d $outdir ]] && mkdir $outdir
cd $outdir

dnames="K36-K27 K36-K9m3 K36-K4m1 K9m3 K4m1 K36m3 K27"

for d in $dnames; do
    mkdir $d
    mv $indir/*${d}*.txt $d
done
