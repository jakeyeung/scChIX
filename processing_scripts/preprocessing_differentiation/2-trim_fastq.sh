#!/bin/sh
# Jake Yeung
# 2-trim_fastq.sh
#  
# 2022-01-21

# WRAP UP
while [[ `squeue -u jyeung | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/raw_data/OUD7117_to_OUD7119_OUD7150/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

jmem='16G'
jtime='18:00:00'

for indir in `ls -d $inmain/PZ*`; do
    echo $indir
    bname=$(basename $indir)
    BNAME=$indir/${bname}.qsub.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    f1=$indir/demultiplexedR1.fastq.gz
    f2=$indir/demultiplexedR2.fastq.gz
    f1out=$indir/trimmed.R1.fastq.gz
    f2out=$indir/trimmed.R2.fastq.gz
    [[ -e $f1out ]] && echo "$f1out  found, continuing" && continue
    [[ -e $f2out ]] && echo "$f2out  found, continuing" && continue
    echo "Trimming $f1 $f2"

    cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; cutadapt -o $f1out -p $f2out $f1 $f2 -m 3 -a 'IlluminaSmallAdapterConcatBait=GGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT' -a 'IlluminaIndexAdapter=GGAATTCTCGGGTGCCAAGGAACTCCAGTCACN{6}ATCTCGTATGCCGTCTTCTGCTTG'  -A 'IlluminaPairedEndPCRPrimer2.0=AGATCGGAAGAGCGN{6}CAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -A 'universalPrimer=GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=5' -a  'IlluminaGEX=TTTTTAATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC;min_overlap=5' -a 'IlluminaMultiplexingPCRPrimer=GGAACTCCAGTCACN{6}TCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -A 'Aseq=TGGCACCCGAGAATTCCA' -a 'Aseq=TGGCACCCGAGAATTCCA'  -a 'illuminaSmallRNAAdapter=TCGTATGCCGTCTTCTGCTTGT'"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done
