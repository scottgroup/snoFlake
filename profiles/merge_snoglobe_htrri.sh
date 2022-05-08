#!/bin/bash

set -o pipefail

ml python/3.8 &&
ml scipy-stack &&
source /home/kris98/venv38/bin/activate &&

snoglobe_dir=$SCRATCH/network_analysis/network_data/snoglobe_predictions
htrri_file=$SCRATCH/network_analysis/network_data/htrri_interactions_P_L_S.tsv
sno_chromo_file=$SCRATCH/network_analysis/profiles_data/sno8_chromo.txt
snodb=$SCRATCH/network_analysis/network_data/snoDB.tsv
outdir=$SCRATCH/network_analysis/profiles_data/snoglobe_htrri_merged

mkdir -p $outdir &&

python3 merge_snoglobe_htrri.py $snoglobe_dir $htrri_file $sno_chromo_file $snodb $outdir &&

echo 'Done'
