#!/bin/sh
# Jake Yeung
# 0-copy_counts_from_MF.sh
#  
# 2021-06-28

indir="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_bwa/embryo_B6CastMaskedGenome_noAlleles/embryo_bwa_B6Cmask"

outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/from_MF"

cd $outdir

# all countTables
cp -r $indir/countTables* $outdir

