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

# Find snoRNAs that bind to pre-mRNA trancripts that encode for RBPs

# path variables
outfile=$SCRATCH/p_3/output/sno_bind_to_rbp_transcript.tsv
sno=$SCRATCH/p_3/snoglobe_formatted

touch $outfile 
echo -e "snoRNA\tRBP\tinteraction" >> $outfile

# bedtools intersect between snoGloBe & RBP transcripts
for f in $sno/*; do 
    bedtools intersect -a "$f" -b $SCRATCH/p_3/rbp_transcripts.tsv -wo -s | awk '{print $4"\t"$10"\tsno_rbp_transcript"}' | uniq >> $outfile
done

echo 'Done'