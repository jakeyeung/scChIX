#!/bin/sh
# Jake Yeung
# 2-trim_fastqs.sh
#  
# 2021-06-22

jmem='4G'
jtime='24:00:00'

inmain0="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation"
[[ ! -d $inmain0 ]] && echo "$inmain0 inmain not found, exiting" && exit 1
d1="demux_fastqs_singles"
# d2="demux_fastqs_doubles"

for d in $d1; do
    inmain=${inmain0}/${d}
    dprefix="trimmed_${d}"

    outmain="${inmain0}/${dprefix}"
    [[ ! -d $outmain ]] && mkdir $outmain

    for indir in `ls -d $inmain/E8*`; do
        bname=$(basename $indir)
        outdir=${outmain}/${bname}
        [[ ! -d $outdir ]] && mkdir $outdir

        BNAME=${outdir}/${bname}.sbatch.log
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE DBASE not found, exiting" && exit 1
        f1=$indir/demultiplexedR1.fastq.gz
        f2=$indir/demultiplexedR2.fastq.gz
        f1out=$outdir/trimmed.R1.fastq.gz
        f2out=$outdir/trimmed.R2.fastq.gz
        [[ -e $f1out ]] && echo "$f1out  found, continuing" && continue
        [[ -e $f2out ]] && echo "$f2out  found, continuing" && continue
        echo "Trimming $f1 $f2"
        echo ${BNAME}
        # exit 0
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo3; cutadapt -o $f1out -p $f2out $f1 $f2 -m 3 -a 'IlluminaSmallAdapterConcatBait=GGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT' -a 'IlluminaIndexAdapter=GGAATTCTCGGGTGCCAAGGAACTCCAGTCACN{6}ATCTCGTATGCCGTCTTCTGCTTG'  -A 'IlluminaPairedEndPCRPrimer2.0=AGATCGGAAGAGCGN{6}CAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -A 'universalPrimer=GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT;min_overlap=5' -a  'IlluminaGEX=TTTTTAATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC;min_overlap=5' -a 'IlluminaMultiplexingPCRPrimer=GGAACTCCAGTCACN{6}TCTCGTATGCCGTCTTCTGCTTG;min_overlap=5' -A 'Aseq=TGGCACCCGAGAATTCCA' -a 'Aseq=TGGCACCCGAGAATTCCA'  -a 'illuminaSmallRNAAdapter=TCGTATGCCGTCTTCTGCTTGT'"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
    done
done
