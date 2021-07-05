#!/bin/sh
# Jake Yeung
# 5-make_count_tables_from_tagged_bams.sh
#  
# 2021-06-27


# # WRAP UP
# while [[ `squeue -u jyeung | wc -l` > 1 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='16G'
jtime='24:00:00'

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams"
inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/tagged_bams"

mapq=40
# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nospikeins.nochromo.bed"
# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.SpikeIns.nochromo.bed"

binsizes="10000 50000"

# binsize=10000
for binsize in $binsizes; do
    stepsize=$binsize

    # echo $jmark
    outdir="${inmain}/counts_tables_${binsize}"
    [[ ! -d $outdir ]] && mkdir $outdir

    for inbam in `ls -d $inmain/*.bam`; do
        inbambase=$(basename $inbam)
        inbambase=${inbambase%.*}
        echo $inbambase
        bname=${inbambase}

        BNAME=$outdir/${bname}.sbatchout
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "dbase $DBASE not found, exiting" && exit 1

        outf1=$outdir/${bname}.countTable.binsize_${binsize}.csv
        [[ -e $outf1 ]] && echo "outf1 $outf1 found, continuing" && continue

        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py $inbam -sliding $stepsize --filterXA -minMQ $mapq -o $outf1 -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup --r1only -blacklist $bl --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
    done

done
