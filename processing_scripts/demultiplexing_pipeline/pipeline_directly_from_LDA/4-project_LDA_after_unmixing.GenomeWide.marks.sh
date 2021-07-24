#!/bin/sh
# Jake Yeung
# 1-project_new_samples_on_LDA.H3K27me3.StemCells.sh
#  
# 2019-10-24

jmem='24G'
jtime='24:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# WRAP UP
while [[ `squeue -u jyeung | grep RunFits | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

sleep 600

# WRAP UP
while [[ `squeue -u jyeung | grep unmix | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# run both marks

# jmarks="K36 K27"
# jmarkdbl="K36-K27"
# jstr="K36_K27_K36-K27"

# jmarks="K36 K9m3"
jmark1="K36"
# jmark2="K9m3"
# jmark2="K27"

jmark2vec="K27 K9m3"
# jmark2vec="K27"

# jappend="TopBinsFilt"
jname="var_filtered"
# jquant="0.15"
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

jquant="manual2nocenter"
jdate="2021-07-23"

for jmark2 in $jmark2vec; do

  jmarks="${jmark1} ${jmark2}"
  jmarkdbl="${jmark1}-${jmark2}"
  jstr="${jmark1}_${jmark2}_${jmarkdbl}"
  
  
  jsuffix=${jname}_${jquant}_${jstr}
  
  ldasuffix="${jsuffix}_${jstr}"
  
  prefix="scchix_inputs_clstr_by_celltype_${jmarkdbl}.removeNA_FALSE"
  
  outmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_scchix_outputs/from_pipeline/${jsuffix}"
  [[ ! -d $outmain ]] && mkdir $outmain
  
  outdir="${outmain}/${jstr}"
  [[ ! -d $outdir ]] && mkdir $outdir
  
  for jmark in $jmarks; do
  
    BNAME=${outdir}/${jmark}_qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
    inlda="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA/${jsuffix}/lda_output_filt.${jmark}.${jdate}.rds"
    inmat="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/scchix_unmixing_downstream/${jsuffix}/${prefix}-unmixed_mat.${jmark}.rds"
    
    [[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
    [[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1
  
    outf="${outdir}/Gastru_Unmixed_DblMark.${jsuffix}.${jmark}.RData"
    [[ -e $outf ]] && echo "$outf found, exiting" && exit 1
    
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $inmat $outf --RemoveEmptyCells"
    echo $cmd
     sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=project_${jmark} --wrap "$cmd"
  done

done


