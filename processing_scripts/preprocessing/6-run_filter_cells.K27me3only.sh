#!/bin/sh
# Jake Yeung
# 6-run_filter_cells.sh
#  
# 2021-06-27

jmem='16G'
jtime='2:00:00'

# jmarks="K9m3 K36 K36-K9m3"
# jmarks="K27 K36-K27"
# jmarks="K36-K27"
# jmark="K9m3"
# jmarks="K4m1 K36-K4m1"
# j§
# jmarks="K36 K27 K9m3 K4m1 K36-K9m3 K36-K27 K36-K4m1"
jmarks="K27"
countcutoffmin="1000"
countcutoffmax=50000
# varcutoffmin=0.1  # too strict?
varcutoffmin=0.03
varcutoffmax=1
TAcutoff=0.5
chromoskeep="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"
jsuffix="50000"

# rs="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic-functions/scripts/processing_scripts/filter_cells_by_totalcuts_TAfrac_intrachromovar.R"
rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/filter_cells_by_totalcuts_TAfrac_intrachromovar.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

for jmark in $jmarks; do

  indir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/counts_tables_${jsuffix}/${jmark}/cbind_out"
  inf="${indir}/${jmark}.countTable.binsize_${jsuffix}.rds"
  [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

  outdir=${indir}/filtered_counts_varcutoffmin_${varcutoffmin}
  [[ ! -d $outdir ]] && mkdir $outdir

 BNAME=${outdir}/${jmark}.filtered_counts_sbatch
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
  
  infilerz="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/count_tables/RZ_counts/${jmark}/merged_output/${jmark}_RZ_counts_output.rds"
  [[ ! -e $infilerz ]] && echo "$infilerz not found, exiting" && exit 1

  infcount="${inf}"

  echo $jmark

  # skip if already exists

  
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -countcutoffmin $countcutoffmin -countcutoffmax $countcutoffmax -varcutoffmin $varcutoffmin -varcutoffmax $varcutoffmax -TAcutoff $TAcutoff -infilerz $infilerz -infilecounts $infcount -names ${jmark} -outdir $outdir -chromoskeep $chromoskeep"
  # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -countcutoffmin $countcutoffmin -countcutoffmax $countcutoffmax -varcutoffmin $varcutoffmin -varcutoffmax $varcutoffmax -TAcutoff $TAcutoff -infilerz $infilerz -infilecounts $infcount -names ${jmark} -outdir $outdir -chromoskeep $chromoskeep
  # exit 0
  
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jmark} --wrap "$cmd"
  
done

# /hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/Rscript /hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic-functions/scripts/processing_scripts/filter_cells_by_totalcuts_TAfrac_intrachromovar.R -countcutoffmin 1000 500 1000 700 500 500 500 -countcutoffmax 30000 10000 50000 50000 50000 50000 50000 -varcutoffmin 0.1 0.1 0.1 0.075 0.1 0.1 0.75 -varcutoffmax 1 1 1 1 1 1 1 -TAcutoff 0.5 -infilerz /hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/embryo_bwa_B6Cmask/countTables/H3K4me1_E6_E9p5_E10_noBlRegions_MQ50_50kb.csv.gz /hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/embryo_bwa_B6Cmask/countTables/H3K4me1_E6_E9p5_E10_noBlRegions_MQ50_50kb.csv.gz /hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/embryo_bwa_B6Cmask/countTables/H3K36me3_E7_E7p5_E8_E9p5_E10_E10p5_E11p5_noBlRegions_MQ50_50kb.csv.gz /hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/embryo_bwa_B6Cmask/countTables/H3K27me3_E7-E8-E9p5-E11p5_noBlRegions_MQ50_50kb.csv.gz /hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/embryo_bwa_B6Cmask/countTables/K36me3_K27me3_E9p5_E10p5_E11p5_noBlRegions_MQ50_50kb.csv.gz /hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/embryo_bwa_B6Cmask/countTables/K36me3_K4me1_E9p5_noBlRegions_MQ50_50kb.csv.gz /hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/embryo_bwa_B6Cmask/countTables/K36me3_K9me3_E9p5_E10p5_E11p5_noBlRegions_MQ50_50kb.csv.gz -names H3K4me1 H3K36me3 H3K27me3 H3K9me3 H3K36me3-H3K27me3 H3K36me3-H3K4me1 H3K36me3-H3K9me3 -chromoskeep chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 -outdir /hpc/hub_oudenaarden/mflorescu/data/mnase/mm/mergedBAMs/embryo_bwa_B6Cmask/LDAinput/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.lessstringent.2020-04-18 --overwrite
