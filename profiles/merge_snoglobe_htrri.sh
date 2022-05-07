#!/bin/bash

#SBATCH --account=def-scottmic
#SBATCH --mail-user=kristina.song@usherbrooke.ca
#SBATCH --mail-type=END,FAIL
#SBATCH --time=10:00:00
#SBATCH --output=%u.%x-%A[%a].out
#SBATCH --error=%u.%x-%A[%a].err
#SBATCH --mem=32000M
#SBATCH --cpus-per-task=16

ml bedtools &&
source &&



echo 'Done'