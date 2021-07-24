#!/bin/sh
# Jake Yeung
# 0-run.celltyping_from_LDA_gastru.sh
#  
# 2021-07-17

jmem='16G'
jtime='3:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/infer_celltype_from_LDA_late_gastrulation.R"

inmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/objs_from_LDA"
outmain="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/scchix_pipeline_from_LDA/celltyping_from_LDA"

jprefix="var_filtered"
# jappend="TopBinsFilt"
# jquant="0.15"
# jquant="0.2"

# jquantvec="0.15 0.2"
# jdate="2021-07-16"

# jquantvec="0.3"
# jdate="2021-07-17"

# jquantvec="manual2"
# jdate="2021-07-20"

# jquantvec="manual2noblood"
# jdate="2021-07-22"

jquantvec="manual2nocenter"
jdate="2021-07-23"

jmark1="K36"
# jmark2="K9m3"
# jmark2="K27"

jmark2vec="K27 K9m3"

for jquant in $jquantvec; do

  for jmark2 in $jmark2vec; do
  
    jmarkdbl="${jmark1}-${jmark2}"
    jstr="${jmark1}_${jmark2}_${jmarkdbl}"
    
    # dname="${jprefix}_${jquant}_${jappend}_${jstr}"
    dname="${jprefix}_${jquant}_${jstr}"
    indir="${inmain}/${dname}"
    outdir="${outmain}/${dname}"
    [[ ! -d $outdir ]] && mkdir $outdir
    
    
    indir=${inmain}/${dname}
    
    inlda1=${indir}/"lda_output_filt.${jmark1}.${jdate}.rds"
    inlda2=${indir}/"lda_output_filt.${jmark2}.${jdate}.rds"
    inldadbl=${indir}/"lda_output_filt.${jmarkdbl}.${jdate}.rds"
    
    [[ ! -e $inlda1 ]] && echo "$inlda1 not found, exiting" && exit 1
    [[ ! -e $inlda2 ]] && echo "$inlda2 not found, exiting" && exit 1
    [[ ! -e $inldadbl ]] && echo "$inldadbl not found, exiting" && exit 1
    
    markvec="${jmark1} ${jmark2} ${jmarkdbl}"
    inldavec="${inlda1} ${inlda2} ${inldadbl}"
    
    tssref="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"
    [[ ! -e $tssref ]] && echo "$tssref not found, exiting" && exit 1
    
    pbulkref="/hpc/hub_oudenaarden/jyeung/data/public_data/CaoPijuana_merged_batch_cor.2019-12-03.RData"
    [[ ! -e $pbulkref ]] && echo "$pbulkref not found, exiting" && exit 1
    
    outprefix="${outdir}/celltype_by_topic_pseudobulk.${dname}"
    
    BNAME="${outdir}/celltype_by_topic_pseudobulk.${dname}.sbatch_output"
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -ldainputvec $inldavec -marks $markvec -pseudobulk_reference ${pbulkref} -keeptop 150 -gene2tss_reference $tssref -outprefix $outprefix --calculate_var"
    
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${dname} --wrap "$cmd"
  
  done
  
done


