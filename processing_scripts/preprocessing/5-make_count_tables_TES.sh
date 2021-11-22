#!/bin/sh
# Jake Yeung
# 6-make_count_tables_TES.sh
#  
# 2021-06-27

jmem='16G'
jtime='24:00:00'

inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/tagged_bams"

mapq=40
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nospikeins.nochromo.bed"

# bedfile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/TES_genomewide.neg_strand_bug_fixed.2021-06-27.txt"
bedfile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/TES_genomewide.neg_strand_bug_fixed.2021-06-27.50kb_max_length.txt"
[[ ! -e $bedfile ]] && echo "$bedfile not found, exiting" && exit 1

# outdir="${inmain}/counts_tables_TSStoTES_50kb_max"
outdir="${inmain}/counts_tables_genebody50kbmax"
[[ ! -d $outdir ]] && mkdir $outdir

for inbam in `ls -d $inmain/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname}
    outftab=${outdir}/${bname}.count_table_genebody50kbmax.txt
    [[ -e $outftab ]] && echo "$outftab found, continuing" && continue
    BNAME=${outdir}/${bname}_sbatch_TES_genomewide_log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamToCountTable.py --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed ${bedfile} -blacklist $bl --r1only --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=TSS_${bname} --wrap "$cmd"
done
