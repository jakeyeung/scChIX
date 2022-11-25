#!/bin/sh
# Jake Yeung
# 6-make_RZ_counts.sh
#  
# 2019-09-04

# sleep 3600

jmem='8G'
jtime='2:00:00'

inbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/raw_data/OUD7117_to_OUD7119_OUD7150/tagged_bams"
[[ ! -d $inbase ]] && echo "$inbase not found, exiting" && exit 1

outbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/raw_data/OUD7117_to_OUD7119_OUD7150/LH_RZ_counts"
[[ ! -d $outbase ]] && mkdir $outbase

# experi="BM"
jmarks="H3K4me1 H3K4me3 H3K36me3 H3K4me3-H3K36me3 H3K27me3-H3K4me1-H3K36me3"

bl="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/mm10.blacklist.copy.sorted.merged.nospikeins.nochromo.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

for jmark in $jmarks; do
  inmain=${inbase}/${jmark}
  outdir=${outbase}/${jmark}
  mkdir $outdir
  echo $jmark
  for inbam in `ls -d $inmain/*${jmark}*.bam`; do
      bname=$(basename $inbam)
      bname=${bname%.*}
  
      outf=$outdir/${bname}.LH_counts.csv
      [[ -e $outf ]] && echo "$outf found, continuing" && continue
  
      BNAME=$outdir/${bname}.LHcounts.qsub
      DBASE=$(dirname "${BNAME}")
      [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
      cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; bamToCountTable.py $inbam -sampleTags SM -featureTags lh -o $outf --dedup --r1only -blacklist $bl --filterXA -minMQ 40 --proper_pairs_only --no_softclips -max_base_edits 2 --no_indels"
      echo $cmd
      sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
  done
  # wait
done


