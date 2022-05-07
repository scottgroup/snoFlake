#!/bin/bash

#SBATCH --account=def-scottmic
#SBATCH --mail-user=kristina.song@usherbrooke.ca
#SBATCH --mail-type=END,FAIL
#SBATCH --time=24:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --error=%u.%x-%A[%a].err
#SBATCH --mem=32000M
#SBATCH --cpus-per-task=16

set -o pipefail

ml nixpkgs/16.09 &&
ml gcc/6.4.0 &&
ml bedtools &&
ml python/3.6.3 &&
source /home/desg2718/venv36/bin/activate &&

scorelist=(98)
nb_list=(3)

profiles_data=$SCRATCH/network_analysis/profiles_data
snofile=$profiles_data/sno8_chromo.txt
gtf=$MYPROJ/human_ensembl_87_wo_dup_v3.csv


for score in ${scorelist[@]}; do
    echo 'score' $score;
    path=$profiles_data/snoglobe_htrri_merged
    cd $path
    for n in ${nb_list[@]}; do
        echo 'min nb windows:' $n;
        python3 nts_flanking_exons.py $gtf $path $n $snofile $SLURM_CPUS_PER_TASK;
    done;
done &&

echo 'Done'

