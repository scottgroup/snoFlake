#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --time=3:00:00
#SBATCH --mail-user=kristina.song@usherbrooke.ca
#SBATCH -o ./logs/output.%j.out
#SBATCH -e ./logs/output.%j.err

set -o pipefail

module load bedtools/2.30.0 &&
module load python/3.8.10 &&
source /home/kris98/venv38/bin/activate &&

# Get target gene names for ENCODE RBP dataset

# path variables
rbp=$SCRATCH/p_3/rbp_formatted
annotation=$SCRATCH/p_3/all_transcripts.tsv
outfile=$SCRATCH/p_3/output/rbp_targets.tsv

touch $outfile
echo -e "RBP\ttarget\tinteraction" >> $outfile

# bedtools intersect between ENCODE RBP & annotated transcripts
for f in $rbp/*; do
    name=$(basename $f | sed 's/_uniq_regions.bed//g')
    bedtools intersect -a "$f" -b $SCRATCH/p_3/all_transcripts.tsv -s -wo | awk -v var="$name" '{print var"\t"$10"\tRBP"}' | uniq >> $outfile
done


echo 'Done'