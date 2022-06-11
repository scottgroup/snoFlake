#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --time=20:00:00
#SBATCH --mail-user=kristina.song@usherbrooke.ca
#SBATCH -o ./logs/output.%j.out
#SBATCH -e ./logs/output.%j.err

set -o pipefail

module load bedtools/2.30.0 &&
module load python/3.8.10 &&
source /home/kris98/venv38/bin/activate &&

# path variables
rbp_path=$SCRATCH/p_3/rbp_formatted 
sno_path=$SCRATCH/p_3/snoglobe_predictions
p_val_path=$SCRATCH/p_3/output/sno_rbp_overlaps_p_vals
outpath=$SCRATCH/p_3/output/sno_rbp_overlap_targets
mkdir -p $outpath

# create file to store significant interactions (used for network)
sig_p_val=$SCRATCH/p_3/output/significant_ovlp.tsv

for f in $p_val_path/*; do
    awk '($3=="0") && ($4=="100000") {print}' $f >> $sig_p_val
done

# temporarily reformat snoglobe prediction for bedtools intersect
snoglobe_temp=$SCRATCH/p_3/snoglobe_temp
mkdir -p $snoglobe_temp

for pred in $sno_path/*; do
    name=$(basename $pred | sed -e 's/pred_\(.*\).98_3.gene.tsv/\1/')
    cut -f 1-6,12 $pred > $snoglobe_temp/$name.bed
done

# run bedtools intersect for significant snoRNA-RBP pair
cat $sig_p_val | while read line; do
    sno=$(awk '{print $1}' <<< "$line")
    sno_name=$(basename $sno | sed 's/_unique_regions.bed//g')
    rbp=$(awk '{print $2}' <<< "$line")
    rbp_name=$(basename $rbp | sed 's/_uniq_regions.bed//g')
    echo -e "snoRNA\tRBP\ttarget" > $outpath/${sno_name}_${rbp_name}_overlap.tsv
    bedtools intersect -a $snoglobe_temp/$sno_name.bed -b $rbp_path/$rbp -s | awk -v var="$rbp_name" '{print $4"\t"var"\t"$7}' | sort | uniq >> $outpath/${sno_name}_${rbp_name}_overlap.tsv
done

# reformat $sig_p_val to fit other network interaction datasets
outfile=$SCRATCH/p_3/output/significant_sno_rbp_overlap.tsv
echo -e "snoRNA\tRBP\tinteraction" > $outfile
cut -f 1,2 $sig_p_val | sed 's/_unique_regions.bed//g' | sed 's/_uniq_regions.bed//g' | awk '{print $1"\t"$2"\t""sno_rbp_overlap"}' | sort | uniq >> $outfile

# remove temp files
rm $sig_p_val && rm -r $snoglobe_temp

echo 'Done'