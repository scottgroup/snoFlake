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

inpath = os.path.abspath(sys.argv[1]) # path to the interaction files
    sno_fasta = os.path.abspath(sys.argv[2]) # fasta file of snoRNA sequences
    score = sys.argv[3] # score cutoff (used for filename and fig title)
    nb_windows = int(sys.argv[4]) # nb of consecutive windows cutoff
    relplot_path = os.path.abspath(sys.argv[5]) # path to relplot utility from viennaRNA

interactions=
sno_fasta=$SCRATCH/profiles_data/sno8.fa
score=98
nb_windows=3
relplot_path=

echo 'Done'