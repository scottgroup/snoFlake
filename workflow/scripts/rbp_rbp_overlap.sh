#!/bin/bash
#SBATCH --time=48:00:00 # time allocated for each task of the job array
#SBATCH --mem=16000M # max required memory per task
#SBATCH --cpus-per-task=8 #number of cores to use, number of cores can be accessed through $SLURM_CPUS_PER_TASK
#SBATCH --array=[0-121] #tasks from job array to run, index can be access through $SLURM_ARRAY_TASK_ID
#SBATCH -o /home/kris98/scratch/network_analysis/network_output/slurmout/rbp_rbp_overlaps_p_vals/output.%j.out
#SBATCH -e /home/kris98/scratch/network_analysis/network_output/slurmout/rbp_rbp_overlaps_p_vals/output.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kristina.song@usherbrooke.ca

set -o pipefail

module load bedtools/2.30.0 &&
module load python/3.8.10 &&
source /home/kris98/venv38/bin/activate &&

data=$SCRATCH/network_analysis/network_data
rbp_path=$data/ENCODE_RBP_formatted

# copy all input data to $SLURM_TMPDIR
cp -r $rbp_path $SLURM_TMPDIR/
cp $data/rbp.tsv $SLURM_TMPDIR/
cp $data/hg38_Ensembl_V101_Scottlab_2020.tsv $SLURM_TMPDIR/
cp $data/chrNameLength.txt $SLURM_TMPDIR/
cp $SCRATCH/network_analysis/sno_RBP_network/scripts/compute_overlaps.py $SLURM_TMPDIR/

# read the rbplist into a bash array
# rbp name in the 1st column
readarray -t rbplist < <(awk '(NR>1) {print $1}' $SLURM_TMPDIR/rbp.tsv)

# use the index ($SLURM_ARRAY_TASK_ID) to get the correct RBP for this task from the bash array
rbp=${rbplist[$SLURM_ARRAY_TASK_ID]}
# print the rbpid for validation
echo $rbp

#file with rbp binding interactions for current RBP
rbpfile=$rbp_path/${rbp}_uniq_regions.bed;

gtf=$SLURM_TMPDIR/hg38_Ensembl_V101_Scottlab_2020.tsv # gtf file in tsv format
genome=$SLURM_TMPDIR/chrNameLength.txt # chromosome length file
tmp_rbp_path=$SLURM_TMPDIR/ENCODE_RBP_formatted # rbp path in $SLURM_TMPDIR/

# output path in $SLURM_TMPDIR
tmp_outpath=$SLURM_TMPDIR/rbp_rbp_overlaps_p_vals/
# create ouptut path
mkdir -p $tmp_outpath

# ouptut file
tmp_outfile=$tmp_outpath/${rbp}_rbp_overlap.tsv
# create output file if needed
touch $tmp_outfile

# go to $SLURM_TMPDIR/
cd $SLURM_TMPDIR/

# run the python script! you should change script path
# $SLURM_TMPDIR is a temporary directory to write temp files on the compute node
python3 compute_overlaps.py $rbpfile $tmp_rbp_path $gtf $genome $tmp_outfile $SLURM_CPUS_PER_TASK $SLURM_TMPDIR rbp;

# copy output files in $SLURM_TMPDIR to original output directory
outpath=$SCRATCH/network_analysis/network_output/rbp_rbp_overlaps_p_vals/
mkdir -p $outpath
cp $tmp_outfile $outpath/

echo 'Done'