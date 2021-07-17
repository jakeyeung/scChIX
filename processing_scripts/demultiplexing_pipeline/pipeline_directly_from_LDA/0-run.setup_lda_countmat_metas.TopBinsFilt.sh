#!/bin/sh
# Jake Yeung
# 0-setup_lda_countmat_metas.sh
#  
# 2021-07-16

jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/cluster_make_lda_countmat_metadata_objects.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

mark1="K36"
# params to tweak
# mark2="K27"
mark2="K9m3"
jquant="0.2"
# jquant="0.15"
midname="coords_var_filtered"
jappend="TopBinsFilt"

markdbl="${mark1}-${mark2}"

marks="${mark1} ${mark2} ${markdbl}"
jstr="${mark1}_${mark2}_${markdbl}"

outmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA"
suffix="${outmain}/${midname}_${jquant}_${jappend}_${jstr}"
outdir="${outmain}/${suffix}"
[[ ! -d $outdir ]] && mkdir $outdir

inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/LDA_outputs"
ldadname="ldaAnalysis_50000_${jstr}"
ldasuffix1="lda_outputs.countmat_var_filt_${jappend}.${mark1}.2021-07-15.K-30.binarize.FALSE/ldaOut.countmat_var_filt_${jappend}.${mark1}.2021-07-15.K-30.Robj"
ldasuffix2="lda_outputs.countmat_var_filt_${jappend}.${mark2}.2021-07-15.K-30.binarize.FALSE/ldaOut.countmat_var_filt_${jappend}.${mark2}.2021-07-15.K-30.Robj"
ldasuffix3="lda_outputs.countmat_var_filt_${jappend}.${markdbl}.2021-07-15.K-30.binarize.FALSE/ldaOut.countmat_var_filt_${jappend}.${markdbl}.2021-07-15.K-30.Robj"

inlda1=${inmain}/${ldadname}/${ldasuffix1}
inlda2=${inmain}/${ldadname}/${ldasuffix2}
inlda3=${inmain}/${ldadname}/${ldasuffix3}
inldas="$inlda1 $inlda2 $inlda3"

for inlda in $inldas; do
  [[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1
done


BNAME=${outdir}/${jstr}_${midname}_${jquant}.sbatch_output
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -marks $marks -infiles $inldas -outdir ${outdir}"

sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jstr}_${midname}_${jquant} --wrap "$cmd"

