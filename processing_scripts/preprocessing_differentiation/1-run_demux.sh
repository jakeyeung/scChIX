#!/bin/sh
# Jake Yeung
# 1-run_demux.sh
# 2021-10-29

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs"
# inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data"
inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChIX_differentiation/raw_data/OUD7117_to_OUD7119_OUD7150"

jmem='16'  # no suffix for SCMO
jtime='24'
BNAME=$inmain/demux
[[ ! -d $BNAME ]] && mkdir $BNAME
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cd $inmain
. "/nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh"; conda activate scmo2022
demux.py *.fastq.gz -hd 0 --cluster -mem $jmem -time $jtime -sched slurm
