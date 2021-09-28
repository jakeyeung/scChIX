#!/bin/sh
# Jake Yeung
# init_output_dirs.sh
#  
# 2021-09-12

workdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/simulation_data/snakemake_outputs"
cd $workdir

mkdir scchix_inputs_objs scchix_outputs_objs LDA_outputs_init objs_from_LDA scchix_unmixing_downstream
