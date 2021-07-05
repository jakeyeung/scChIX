#!/bin/sh
# Jake Yeung
# 6-make_RZ_counts.sh
#  
# 2019-09-04

# sleep 3600

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed"

jmark="K27"
inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/tagged_bams"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/RZcounts.NewFilters"
outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/RZ_counts/${jmark}"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='8G'
jtime='2:00:00'

for inbam in `ls -d $inmain/*${jmark}*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}

    # remove if contains K36
    if `echo $bname | grep -q "K36"`; then
        echo "K36 found in $inbambase, skipping"
        continue
    else
        echo "K36 found in $inbambase, keeping..."
    fi

    outf=$outdir/${bname}.LH_counts.csv

    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$outdir/${bname}.LHcounts.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam -sampleTags SM -featureTags lh -o $outf --dedup --filterXA -minMQ 40 --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
    echo $cmd
    # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam -sampleTags SM -featureTags lh -o $outf --dedup --filterXA -minMQ 40 --proper_pairs_only"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done
wait
