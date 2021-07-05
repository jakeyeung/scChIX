#!/bin/sh
# Jake Yeung
# 4-sort_index_tag_bam.sh
#  
# 2019-12-19

# # # WRAP UP
# # while [[ `qstat | grep PZ |  wc -l` > 0 ]]; do
# #         echo "sleep for 60 seconds"
# #         sleep 60
# # done
# 
# echo "Sleep for 7200 seconds"
# sleep 7200

# WRAP UP
while [[ `squeue -u jyeung | grep Bwa |  wc -l` > 0 ]]; do
    echo "sleep for 66 seconds"
    sleep 66
done

jmem='64G'
jtime='48:00:00'
ncores=8

inbase="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation"
dnames="trimmed_demux_fastqs_singles trimmed_demux_fastqs_doubles"
outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/tagged_bams"
[[ ! -d $outdir ]] && mkdir $outdir

wd="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/workdir"
cd $wd

# clusterdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/clusterfiles"
for d in $dnames; do
    inmain=${inbase}/${d}

    [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
    tmpdir=$inmain/tmpdirmulti
    [[ ! -d $tmpdir ]] && mkdir $tmpdir
    # echo $tmpdir
    # exit 0

    for indir in `ls -d $inmain/E*`; do
        bname=$(basename $indir)
        inbam=$indir/bwaMapped.bam
        [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
        BNAME=$outdir/$bname.sort_index_tag.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        sortedbam=$indir/bwaMapped.sorted.bam
        outbamtagged=$outdir/${bname}.sorted.tagged.bam
        [[ -e $outbamtagged ]] && echo "$outbamtagged found, continuing" && continue
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam; cd $wd; bamtagmultiome.py -method chic --cluster -clusterdir $tmpdir -o $outbamtagged -mem 64 -time 48 --multiprocess $sortedbam"
        echo $cmd
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
    done
done

