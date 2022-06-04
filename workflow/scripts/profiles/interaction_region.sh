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
ml StdEnv/2020 && 
ml viennarna/2.4.17 &&
source /home/kris98/venv38/bin/activate &&

interactions=$SCRATCH/network_analysis/profiles_data/snoglobe_htrri_merged
sno_fasta=$SCRATCH/network_analysis/profiles_data/sno8.fa
score=98
nb_windows=3
relplot_path=$SCRATCH/ViennaRNA-2.5.0/src/Utils/relplot.pl
outpath=$SCRATCH/network_analysis/profiles_output

mkdir -p $outpath &&

python3 interaction_region.py $interactions $outpath $sno_fasta $score $nb_windows $relplot_path &&


echo 'Done'