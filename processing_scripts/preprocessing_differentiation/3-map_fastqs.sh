#!/bin/sh
# Jake Yeung
# 3-map_fastq.sh
#  
# 2022-01-06

# jmem='64G'
# jtime="30:00:00"
# ncores=8

jmem='72G'
jtime="24:00:00"
ncores=16

# WRAP UP
while [[ `squeue -u jyeung | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

bwabin="bwa"
indxbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/genomes/primary_assembly_97_129B6Masked_ERCC92.fa"
# inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/raw_demultiplexed"
inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/raw_data/OUD7117_to_OUD7119_OUD7150/raw_demultiplexed"

[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

for indir in `ls -d $inmain/PZ*`; do
    bname=$(basename $indir)
    outdir=$indir

    f1=$indir/trimmed.R1.fastq.gz
    f2=$indir/trimmed.R2.fastq.gz
    [[ ! -e $f1 ]] && echo "$f1 not found, continuing" && continue
    [[ ! -e $f2 ]] && echo "$f2 not found, continuing" && continue

    BNAME=$outdir/$bname.mapping.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=$outdir/bwaMapped.bam
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; $bwabin mem -t $ncores $indxbase $f1 $f2 | samtools view -Sb - > $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=$ncores --nodes=1 --job-name=BwaMap_${bname} --wrap "$cmd"
done
