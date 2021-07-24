#!/bin/sh
# Jake Yeung
# 6-run.cbind_count_mats.sh
#  
# 2021-06-27

jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/cbind_count_mats.R"
# jmark="K9m3"
# jmark="K36"
# jmark="K27"

# jsuffix="counts_tables_TSStoTES_50kb_max"
jsuffix="counts_tables_topfeatures_K36_genebodies_K9m3_bins"

# jmarks="K4m1"

jmarks="K36 K9m3 K36-K9m3"

for jmark in $jmarks; do
  indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/${jsuffix}/${jmark}"
  [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
  
  outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/${jsuffix}/${jmark}/cbind_out"
  [[ ! -d $outdir ]] && mkdir $outdir
  
  outf=${outdir}/${jmark}.${jsuffix}.rds
  
  BNAME=${outdir}/${jmark}.sbatch_log
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
  infs=${indir}/*.txt
  
  # echo $infs
  
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $infs -outfile $outf --from_bed"
  echo $cmd
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jmark} --wrap "$cmd"
done
