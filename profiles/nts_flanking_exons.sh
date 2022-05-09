#!/bin/bash

#SBATCH --account=def-scottmic
#SBATCH --mail-user=kristina.song@usherbrooke.ca
#SBATCH --mail-type=END,FAIL
#SBATCH --time=24:00:00
#SBATCH --mem=32000M
#SBATCH --cpus-per-task=16
#SBATCH --output=/home/kris98/scratch/network_analysis/logs/%j_%x.out
#SBATCH --error=/home/kris98/scratch/network_analysis/logs/%j_%x.err

#### Adapted from Gabrielle's script ####

set -o pipefail

ml python/3.8 &&
ml scipy-stack &&
ml bedtools/2.30.0 &&
source /home/kris98/venv38/bin/activate &&

scorelist=(98)
nb_list=(3)

profiles_data=$SCRATCH/network_analysis/profiles_data
snofile=$profiles_data/sno8_chromo.txt
gtf=$SCRATCH/network_analysis/network_data/hg38_Ensembl_V101_Scottlab_2020.tsv
outpath=$SCRATCH/network_analysis/profiles_output
inpath=$profiles_data/snoglobe_htrri_merged

mkdir -p $outpath &&

for score in ${scorelist[@]}; do
    echo 'score' $score;
    for n in ${nb_list[@]}; do
        echo 'min nb windows:' $n;
        python3 nts_flanking_exons.py $gtf $inpath $outpath $n $snofile $SLURM_CPUS_PER_TASK;
    done;
done &&

echo 'Done'

