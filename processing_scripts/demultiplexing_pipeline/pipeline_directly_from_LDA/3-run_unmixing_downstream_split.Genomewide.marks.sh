#!/bin/sh
# Jake Yeung
# 3b-run_unmixing_downstream_split_reads.sh
#  
# 2020-02-06

# echo "Sleep 2 hours"
# sleep 7200

# WRAP UP
while [[ `squeue -u jyeung | grep RunFits | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# # WRAP UP
# while [[ `qstat | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='64G'
jtime='3:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/unmixing_downstream_split_reads_script.UnevenRows.R"

mark1="K36"
# mark2vec="K27 K9m3"
# mark2vec="K27"
mark2vec="K9m3"

# mark2="K27"
# mark2="K9m3"
# jquant="0.2"
# jdate="2021-07-16"

# jquant="0.3"
# jdate="2021-07-17"

# jquant="manual"
# jdate="2021-07-19"

# jquant="manual2"
# jdate="2021-07-20"

# jquant="manual2noblood"
# jdate="2021-07-22"

# jquant="manual2nocenter"
# jdate="2021-07-23"

jquant="manual2nocenternoneu"
jdate="2021-08-03"

for mark2 in $mark2vec; do

  markdbl="${mark1}-${mark2}"  # match lda but not input in Rscript
  
  jprefix="var_filtered"
  jstr="${mark1}_${mark2}_${markdbl}"
  inputname="scchix_inputs_objs"
  outputname="scchix_outputs_objs"
  unmixname="scchix_unmixing_downstream"
  
  dname="${jprefix}_${jquant}_${jstr}"
  
  
  inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA"
  inputdir="${inmain}/${inputname}/${dname}"
  outputdir="${inmain}/${outputname}/${dname}"
  outdir="${inmain}/${unmixname}/${dname}"
  [[ ! -d $outdir ]] && mkdir $outdir
  
  [[ ! -d $inputdir ]] && echo "$inputdir not found, exiting" && exit 1
  [[ ! -d $outputdir ]] && echo "$outputdir not found, exiting" && exit 1
  
  unmixingname="scchix_inputs_clstr_by_celltype_${markdbl}.removeNA_FALSE"
  infdblinput="${inputdir}/scchix_inputs_clstr_by_celltype_${markdbl}.removeNA_FALSE.RData"
  infdbloutput="${outputdir}/unmix_scchix_inputs_clstr_by_celltype_${markdbl}.removeNA_FALSE.RData"
  
  [[ ! -e $infdblinput ]] && echo "$infdblinput not found, exiting" && exit 1
  [[ ! -e $infdbloutput ]] && echo "$infdbloutput not found, exiting" && exit 1
  
  inbase="${inmain}/objs_from_LDA"
  inmain="${inbase}/${dname}"
  
  inflda1="${inmain}/lda_output_filt.${mark1}.${jdate}.rds"
  inflda2="${inmain}/lda_output_filt.${mark2}.${jdate}.rds"
  infldadbl="${inmain}/lda_output_filt.${markdbl}.${jdate}.rds"
  
  [[ ! -e $inflda1 ]] && echo "$inflda1 not found, exiting" && exit 1
  [[ ! -e $inflda2 ]] && echo "$inflda2 not found, exiting" && exit 1
  [[ ! -e $infldadbl ]] && echo "$infldadbl not found, exiting" && exit 1
  
  infmat1="${inmain}/countmat_output_filt.${mark1}.${jdate}.rds"
  infmat2="${inmain}/countmat_output_filt.${mark2}.${jdate}.rds"
  infmatdbl="${inmain}/countmat_output_filt.${markdbl}.${jdate}.rds"
  
  [[ ! -e $infmat1 ]] && echo "$infmat1 not found, exiting" && exit 1
  [[ ! -e $infmat2 ]] && echo "$infmat2 not found, exiting" && exit 1
  [[ ! -e $infmatdbl ]] && echo "$infmatdbl not found, exiting" && exit 1
  
  outprefix="${outdir}/${unmixingname}"
  
  BNAME=$outprefix.qsub
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -mark1 $mark1 -mark2 $mark2 -inf_dbl_input $infdblinput -inf_dbl_output $infdbloutput -inf_mark1 $infmat1 -inf_mark2 $infmat2 -inf_dblmark_lda $infldadbl -outprefix $outprefix"
  echo $cmd
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=unmix_ds_${markdbl} --wrap "$cmd"
  


done


