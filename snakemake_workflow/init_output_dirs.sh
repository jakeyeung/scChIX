#!/bin/sh
# Jake Yeung
# init_output_dirs.sh
#  
# 2021-09-12

# workdir="/home/hub_oudenaarden/jyeung/projects/scChIX_ES_NPC/snakemake_pipeline"
workdir=$1
[[ ! -d $workdir ]] && echo "$workdir not found, exiting" && exit 1
[[ ! -d $workdir ]] && mkdir $workdir
cd $workdir

mkdir scchix_inputs_objs scchix_outputs_objs LDA_outputs_init objs_from_LDA scchix_unmixing_downstream
