#!/bin/sh
# Jake Yeung
# 6-run.rbind_RZ_files.sh
#  
# 2021-06-28

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/processing_scripts/rbind_RZ_files.R"

# jmark=""
# jmarks="K36 K4m1 K36-K9m3 K36-K4m1 K36-K27m3 K9m3"
jmarks="K27"

for jmark in $jmarks; do
  indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/RZ_counts/${jmark}"
  infs=${indir}/*.csv
  outdir="${indir}/merged_output"
  [[ ! -d $outdir ]] && mkdir $outdir
  outf=${outdir}/${jmark}_RZ_counts_output.rds
  [[ -e $outf ]] && echo "$outf  found, continuing" && continue

  BNAME=${outdir}/${jmark}_rbind_sbatch
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $infs -outfile $outf"
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jmark} --wrap "$cmd"
done
