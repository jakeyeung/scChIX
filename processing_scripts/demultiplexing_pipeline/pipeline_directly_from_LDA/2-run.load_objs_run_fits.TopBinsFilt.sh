#!/bin/sh
# Jake Yeung
# 2b-run.load_objs_run_fits.from_server.sh
#  
# 2019-08-08

# echo "Sleep an hour"
# sleep 3600

# WRAP UP
while [[ `squeue -u jyeung | grep Setup |  wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='32G'
jtime='24:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/load_objs_run_fits.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

basedir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA"

jname="coords_var_filtered"
# jquant="0.2"
# jquant="0.15"
# jquant="0.3"
jquant="manual2"
jappend="TopBinsFilt"

# markdbl="K36-K9m3"
mark1="K36"
# # mark2="K27"
# mark2="K9m3"

mark2vec="K27 K9m3"

for mark2 in $mark2vec; do
  markdbl="${mark1}-${mark2}"
  jstr=${mark1}_${mark2}_${markdbl}
  jdir=${jname}_${jquant}_${jappend}_${jstr}
  
  maindir="${basedir}/scchix_inputs_objs/${jdir}"
  cd $maindir
  
  ncores=8
  jmethod="Brent"
  
  outdir="${basedir}/scchix_outputs_objs/${jdir}"
  
  [[ ! -d ${outdir} ]] && mkdir ${outdir}
  
  indat="${maindir}/scchix_inputs_clstr_by_celltype_${markdbl}.removeNA_FALSE.RData"
  [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
  
  bname=$(basename $indat)
  bname=${bname%.*}
  
  [[ ! -d $outdir ]] && mkdir -p $outdir
  outdat="${outdir}/unmix_${bname}.RData"
  
  BNAME=$outdir/$bname.qsub
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
  [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
  [[ -e $outdat ]] && echo "$outdat found, exiting" && exit 1
  [[ $ncores != [0-9]* ]] && echo "Must be integer: $ncores" && exit 1
  
  cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; cd $maindir; Rscript $rs $indat $outdat --ncores $ncores --method $jmethod"
  echo "--cpus-per-task=${ncores}"
  echo $cmd
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=${ncores} --job-name=RunFits_${markdbl} --wrap "$cmd"
done


