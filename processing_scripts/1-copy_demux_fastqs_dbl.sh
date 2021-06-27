#!/bin/sh
# Jake Yeung
# copy_demux_fastqs.sh
#  
# 2021-06-22

outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/gastrulation/demux_fastqs"

# H3K36me3 and H3K4me1
inmain1="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K36me3_H3K4me1"
dnames1="E9p5-CB6-K36-K4m1-190409-1 E9p5-CB6-K36-K4m1-190409-2"

inmain2="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K36me3_H3K27me3"
dnames2="E9p5-CB6-K36-K27me3-190404-1 E10p5-B6C-K36-K27m3-190409-1 E10p5-B6C-K36-K27m3-190409-2 E11p5-CB6-K36-K27m3-190618-1 E10p5-CB6-K36-K27m3-190409-4 E11p5-CB6-K36-K27m3-190618-4 E11p5-CB6-K36-K27m3-190618-3 E9p5-CB6-K36-K27m3-190404-2 E10p5-CB6-K36-K27m3-190409-3 E11p5-CB6-K36-K27m3-190618-2"

inmain3="/hpc/archive/hub_oudenaarden/userBackups/mflorescu/mm_chic_raw_demultiplexed/H3K36me3_H3K9me3"
dnames3="E9p5-B6C-K36-K9m3-190404-2 E9p5-CB6-K36-K9m3-190409-2 E9p5-CB6-K36-K9m3-190409-1 E10p5-B6C-K36-K9m3-190409-3 E10p5-B6C-K36-K9m3-190409-2 E10p5-B6C-K36-K9m3-190409-1 E10p5-B6C-K36-K9m3-190409-4 E11p5-CB6-K36-K9m3-190618-2 E11p5-CB6-K36-K9m3-190618-1 E11p5-CB6-K36-K9m3-190618-3 E11p5-CB6-K36-K9m3-190618-4"


cd $outdir

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
