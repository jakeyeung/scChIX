#!/bin/sh
# Jake Yeung
# 2b-run.setup_objs_for_unmixing.sh
# Visually look for interesting topics, create umap_annotated... RData objects which 
# are input to setting up for unmixing
# 2020-02-05

# # WRAP UP
while [[ `squeue -u jyeung | grep K36_K27 |  wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/setup_objs_for_unmixing_general_from_LDA_or_countmat.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/scchix_inputs_objs"

inbase="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA"
[[ ! -d $inbase ]] && echo "$inbase not found, exiting" && exit 1

mark1="K36"

mark2vec="K27 K9m3"
# mark2vec="K27"
# mark2="K9m3"
# # mark2="K27"

# jdate="2021-07-16"
# jquant="0.2"

# jdate="2021-07-17"
# jquant="0.3"

# jdate="2021-07-20"
# jquant="manual2"

# jdate="2021-07-22"
# jquant="manual2noblood"

jdate="2021-07-23"
jquant="manual2nocenter"

jname="var_filtered"

for mark2 in $mark2vec; do

  markdbl="${mark1}-${mark2}"
  jstr="${mark1}_${mark2}_${markdbl}"
  suffix2="${jname}_${jquant}_${jstr}"
  inmain=${inbase}/${suffix2}
  
  inf1meta="${inmain}/celltyping_output_filt.${mark1}.${jdate}.rds"
  inf2meta="${inmain}/celltyping_output_filt.${mark2}.${jdate}.rds"
  infdblmeta="${inmain}/celltyping_output_filt.${markdbl}.${jdate}.rds"
  
  inf1countmat="${inmain}/countmat_output_filt.${mark1}.${jdate}.rds"
  inf2countmat="${inmain}/countmat_output_filt.${mark2}.${jdate}.rds"
  infdblcountmat="${inmain}/countmat_output_filt.${markdbl}.${jdate}.rds"
  
  inf1lda="${inmain}/lda_output_filt.${mark1}.${jdate}.rds"
  inf2lda="${inmain}/lda_output_filt.${mark2}.${jdate}.rds"
  
  [[ ! -e $inf1meta ]] && echo "$inf1meta not found, exiting" && exit 1
  [[ ! -e $inf2meta ]] && echo "$inf2meta not found, exiting" && exit 1
  [[ ! -e $infdblmeta ]] && echo "$infdblmeta not found, exiting" && exit 1
  
  [[ ! -e $inf1countmat ]] && echo "$inf1countmat not found, exiting" && exit 1
  [[ ! -e $inf2countmat ]] && echo "$inf2countmat not found, exiting" && exit 1
  [[ ! -e $infdblcountmat ]] && echo "$infdblcountmat not found, exiting" && exit 1
  
  [[ ! -e $inf1lda ]] && echo "$inf1lda not found, exiting" && exit 1
  [[ ! -e $inf2lda ]] && echo "$inf2lda not found, exiting" && exit 1
  
  outdir="${outmain}/${suffix2}"
  [[ ! -d $outdir ]] && mkdir -p $outdir
  
  outprefix="${outdir}/scchix_inputs_clstr_by_celltype_${markdbl}"
  
  BNAME=${outprefix}.qsub
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -mark1 $mark1 -mark2 $mark2 -inf1meta $inf1meta -inf2meta $inf2meta -infdblmeta $infdblmeta -inf1lda $inf1lda -inf2lda $inf2lda -infdblcountmat $infdblcountmat -outprefix $outprefix"
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=SetupObjs_${markdbl} --wrap "$cmd"
  

done


