#!/bin/sh
# Jake Yeung
# 0-setup_lda_countmat_metas.sh
#  
# 2021-07-16

# WRAP UP
while [[ `squeue -u jyeung | grep LDA | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/cluster_make_lda_countmat_metadata_objects.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_outputs"
outmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA"

mark1="K36"
# params to tweak
# mark2="K27"

midname="var_filtered"
# mark2vec="K27 K9m3"
# mark2vec="K27"  # for noE8 we run K27y only 
mark2vec="K9m3"  # K9m3 only when filtering out bad clust


# mark2vec="K27"  # rerun 
# jquantvec="0.3"
# jdate="2021-07-15"

# jquantvec="manual2"
# jdate="2021-07-19"
# jquant="0.2"
# mark2="K9m3"

# jquantvec="manual2noblood"
# jdate="2021-07-21"

# jquantvec="manual2nocenter"
# jdate="2021-07-23"

# jquantvec="manual2nocenternoE8"
# jdate="2021-07-24"

jquantvec="manual2nocenternoneu"
jdate="2021-08-02"

for jquant in $jquantvec; do

  for mark2 in $mark2vec; do
    markdbl="${mark1}-${mark2}"
    marks="${mark1} ${mark2} ${markdbl}"
    jstr="${mark1}_${mark2}_${markdbl}"
    
    ldadname="ldaAnalysis_50000_${midname}_${jquant}_${jstr}"
    ldasuffix1="lda_outputs.countmat_var_filt.${mark1}.${jdate}.K-30.binarize.FALSE/ldaOut.countmat_var_filt.${mark1}.${jdate}.K-30.Robj"
    ldasuffix2="lda_outputs.countmat_var_filt.${mark2}.${jdate}.K-30.binarize.FALSE/ldaOut.countmat_var_filt.${mark2}.${jdate}.K-30.Robj"
    ldasuffix3="lda_outputs.countmat_var_filt.${markdbl}.${jdate}.K-30.binarize.FALSE/ldaOut.countmat_var_filt.${markdbl}.${jdate}.K-30.Robj"
    
    inlda1=${inmain}/${ldadname}/${ldasuffix1}
    inlda2=${inmain}/${ldadname}/${ldasuffix2}
    inlda3=${inmain}/${ldadname}/${ldasuffix3}
    inldas="$inlda1 $inlda2 $inlda3"

    for inldatmp in $inldas; do
      [[ ! -e $inldatmp ]] && echo "$inldatmp not found, exiting" && exit 1
    done
    
    outdir="${outmain}/${midname}_${jquant}_${jstr}"
    [[ ! -d $outdir ]] && mkdir $outdir
    
    BNAME=${outdir}/${jstr}_${midname}_${jquant}.sbatch_output
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -marks $marks -infiles $inldas -outdir ${outdir}"
    
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jstr}_${midname}_${jquant} --wrap "$cmd"
  
  done
  
done


