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

jmem='128G'
jtime='48:00:00'
ncores=1
jmem2=128
jtime2=48

# inbase="/hpc/hub_oudenaarden/jyeung/data/dblchic/differentiation/round1"
# inbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data"
inbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/raw_data/OUD7117_to_OUD7119_OUD7150"
inmain="${inbase}/raw_demultiplexed"

[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
tmpdir=$inmain/tmpdirmulti
[[ ! -d $tmpdir ]] && mkdir $tmpdir

outdir="${inbase}/tagged_bams"
[[ ! -d $outdir ]] && mkdir $outdir

wd="${inbase}/workdir"
[[ ! -d $wd ]] && mkdir $wd
cd $wd

for indir in `ls -d $inmain/PZ*`; do
    bname=$(basename $indir)
    inbam=$indir/bwaMapped.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    BNAME=$outdir/$bname.sort_index_tag.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    sortedbam=$indir/bwaMapped.sorted.bam
    outbamtagged=$outdir/${bname}.sorted.tagged.bam
    [[ -e $outbamtagged ]] && echo "$outbamtagged found, continuing" && continue
    cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam; cd $wd; bamtagmultiome.py -method chic --cluster -clusterdir $tmpdir -o $outbamtagged -mem $jmem2 -time $jtime2 --multiprocess $sortedbam"
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

