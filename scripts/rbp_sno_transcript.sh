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

# Find RBPs that directly bind to snoRNA transcripts

# path variables
outfile=$SCRATCH/p_3/output/rbp_bind_to_sno_transcript.tsv
rbp=$SCRATCH/p_3/rbp_formatted

echo -e "RBP\tsnoRNA\tinteraction" >> $outfile

# bedtools intersect between ENCODE RBP data & snoRNA transcripts
for f in $rbp/*; do 
    name=$(basename $f | sed 's/_uniq_regions.bed//g')
    bedtools intersect -a "$f" -b $SCRATCH/p_3/sno_transcripts.tsv -wo -s | awk -v var="$name" '{print var"\t"$10"\trbp_sno_transcript"}' | uniq >> $outfile
done

echo 'Done'