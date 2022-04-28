#!/bin/bash
#SBATCH --time=48:00:00 # time allocated for each task of the job array
#SBATCH --mem=190000M # max required memory per task
#SBATCH --cpus-per-task=48 #number of cores to use, number of cores can be accessed through $SLURM_CPUS_PER_TASK
#SBATCH --array=[0-121] #tasks from job array to run, index can be access through $SLURM_ARRAY_TASK_ID
#SBATCH -o ./logs/output.%j.out
#SBATCH -e ./logs/output.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kristina.song@usherbrooke.ca

set -o pipefail

module load bedtools/2.30.0 &&
module load python/3.8.10 &&
source /home/kris98/venv38/bin/activate &&

rbp_path=$SCRATCH/p_3/rbp_formatted 

# read the rbplist into a bash array
readarray -t rbplist < $SCRATCH/p_3/rbp.tsv 

# use the index ($SLURM_ARRAY_TASK_ID) to get the correct RBP for this task from the bash array
rbp=${rbplist[$SLURM_ARRAY_TASK_ID]}
#print the rbpid for validation
echo $rbp

#file with rbp binding interactions for current RBP
rbpfile=$rbp_path/${rbp}_uniq_regions.bed; # need to remove redundancy in bed bfile with gene names, if a region overlaps n genes, the region will be written n times. can use cut -f1-6 file_with_geneids.bed | sort -k1,1 -k2,3n -u > unique_regions.bed

gtf=$SCRATCH/p_3/hg38_Ensembl_V101_Scottlab_2020.tsv # gtf file in tsv format
genome=$SCRATCH/p_3/chrNameLength.txt #chromosome length file

#output path
outpath=$SCRATCH/p_3/output/rbp_rbp_overlaps_p_vals/
#create ouptut paht and parent dir if needed
mkdir -p $outpath

#ouptut file
outfile=$outpath/$rbp.tsv
#create output file if needed
touch $outfile

#run the python script! you should change script path
# $SLURM_TMPDIR is a temporary directory to write temp files on the compute node
python3 $SCRATCH/p_3/scripts/sno_rbp_overlap.py $rbpfile $rbp_path $gtf $genome $outfile $SLURM_CPUS_PER_TASK $SLURM_TMPDIR;

echo 'Done'
