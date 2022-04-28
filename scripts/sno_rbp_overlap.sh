#!/bin/bash
#SBATCH --time=20:00:00 # time allocated for each task of the job array
#SBATCH --mem=190000M # max required memory per task
#SBATCH --cpus-per-task=48 #number of cores to use, number of cores can be accessed through $SLURM_CPUS_PER_TASK
#SBATCH --array=[0-252] #tasks from job array to run, index can be access through $SLURM_ARRAY_TASK_ID
#SBATCH -o ./logs/output.%j.out
#SBATCH -e ./logs/output.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kristina.song@usherbrooke.ca

set -o pipefail

module load bedtools/2.30.0 &&
module load python/3.8.10 &&
source /home/kris98/venv38/bin/activate &&

rbp_path=$SCRATCH/p_3/rbp_formatted 

# read the snolist into a bash array
readarray -t snolist < $SCRATCH/p_3/snoRNA.tsv

# use the index ($SLURM_ARRAY_TASK_ID) to get the correct snoRNA for this task from the bash array
sno=${snolist[$SLURM_ARRAY_TASK_ID]}
#print the snoid for validation
echo $sno

#file with snoglobe's predictions for current snoRNA
snofile=$SCRATCH/p_3/snoglobe_formatted/${sno}_unique_regions.bed; # need to remove redundancy in bed bfile with gene names, if a region overlaps n genes, the region will be written n times. can use cut -f1-6 file_with_geneids.bed | sort -k1,1 -k2,3n -u > unique_regions.bed

gtf=$SCRATCH/p_3/hg38_Ensembl_V101_Scottlab_2020.tsv # gtf file in csv format
genome=$SCRATCH/p_3/chrNameLength.txt #chromosome length file

#output path
outpath=$SCRATCH/p_3/output/sno_rbp_overlaps_p_vals/
#create ouptut paht and parent dir if needed
mkdir -p $outpath

#ouptut file
outfile=$outpath/$sno.tsv
#create output file if needed
touch $outfile

#run the python script! you should change script path
# $SLURM_TMPDIR is a temporary directory to write temp files on the compute node
python3 $SCRATCH/p_3/scripts/sno_rbp_overlap.py $snofile $rbp_path $gtf $genome $outfile $SLURM_CPUS_PER_TASK $SLURM_TMPDIR;

echo 'Done'
