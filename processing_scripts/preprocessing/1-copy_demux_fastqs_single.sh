#!/bin/sh
# Jake Yeung
# copy_demux_fastqs.sh
#  
# 2021-06-22

outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/demux_fastqs_singles"
[[ ! -d $outdir ]] && mkdir $outdir

inmain0="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K4me1"
dnames0="E9p5-CB6-K4m1-190404-1 E9p5-CB6-K4m1-190404-2 E10-CB6-K4m1-192201-1 E10-CB6-K4m1-192201-2 E10-CB6-K4m1-192201-3"


# H3K36me3 and H3K4me1
# inmain1="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K36me3_H3K4me1"
inmain1="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K36me3"
dnames1="E9p5-B6C-K36m3-190404-1 E9p5-CB6-K36m3-190404-2 E10-CB6-K36m3-192201-1 E10p5-B6C-K36m3-190409-2 E11p5-CB6-K36-190618-2 E10-CB6-K36m3-192201-2 E10p5-B6C-K36m3-190409-3 E11p5-CB6-K36-190618-3 E10-CB6-K36m3-192201-3 E10p5-B6C-K36m3-190409-4 E11p5-CB6-K36-190618-4 E10p5-B6C-K36m3-190409-1 E11p5-CB6-K36-190618-1"

# inmain2="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K36me3_H3K27me3"
inmain2="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K27me3"
dname2="E9p5-B6C-K27m3-190404-1 E11p5-CB6-K27m3-190618-2 E11p5-CB6-K27m3-190618-3 E11p5-CB6-K27m3-190618-4"

inmain3="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K9me3"
dnames3="E9p5-CB6-K9m3-190409-1 E9p5-CB6-K9m3-190409-2 E10p5-B6C-K9m3-190409-1 E10p5-B6C-K9m3-190409-2 E10p5-B6C-K9m3-190409-3 E10p5-B6C-K9m3-190409-4 E11p5-CB6-K9m3-190618-1 E11p5-CB6-K9m3-190618-2 E11p5-CB6-K9m3-190618-3 E11p5-CB6-K9m3-190618-4"

cd $outdir

for d0 in $dnames0; do
    dfull=${inmain0}/${d0}
    [[ ! -d $dfull ]] && echo "$dfull not found, exiting" && exit 1
    cp -r $dfull .
done

for d1 in $dnames1; do
    dfull=${inmain1}/${d1}
    [[ ! -d $dfull ]] && echo "$dfull not found, exiting" && exit 1
    cp -r $dfull .
done

for d2 in $dnames2; do
    dfull2=${inmain2}/${d2}
    [[ ! -d $dfull2 ]] && echo "$dfull2 not found, exiting" && exit 1
    cp -r $dfull2 .
done

for d3 in $dnames3; do
    dfull3=${inmain3}/${d3}
    [[ ! -d $dfull3 ]] && echo "$dfull3 not found, exiting" && exit 1
    cp -r $dfull3 .
done
