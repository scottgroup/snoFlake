#!/bin/bash
#SBATCH --time=20:00:00 # time allocated for each task of the job array
#SBATCH --mem=16000M # max required memory per task
#SBATCH --cpus-per-task=8 #number of cores to use, number of cores can be accessed through $SLURM_CPUS_PER_TASK
#SBATCH --array=[0] #tasks from job array to run, index can be access through $SLURM_ARRAY_TASK_ID
#SBATCH -o /home/kris98/scratch/network_analysis/network_output/slurmout/sno_sno_overlaps_p_vals/output.%j.out
#SBATCH -e /home/kris98/scratch/network_analysis/network_output/slurmout/sno_sno_overlaps_p_vals/output.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kristina.song@usherbrooke.ca

set -o pipefail

module load bedtools/2.30.0 &&
module load python/3.8.10 &&
source /home/kris98/venv/bin/activate &&

data=$SCRATCH/network_analysis/network_data
sno_path=$data/snoglobe_formatted 

# copy all input data to $SLURM_TMPDIR
cp -r $sno_path $SLURM_TMPDIR/
cp $data/snoRNA_rerun.tsv $SLURM_TMPDIR/
cp $data/hg38_Ensembl_V101_Scottlab_2020.tsv $SLURM_TMPDIR/
cp $data/chrNameLength.txt $SLURM_TMPDIR/
cp $SCRATCH/network_analysis/sno_RBP_network/scripts/compute_overlaps.py $SLURM_TMPDIR/

# read the snolist into a bash array
readarray -t snolist < <(awk '(NR>1) {print $1}' $SLURM_TMPDIR/snoRNA_rerun.tsv)

# use the index ($SLURM_ARRAY_TASK_ID) to get the correct snoRNA for this task from the bash array
sno=${snolist[$SLURM_ARRAY_TASK_ID]}
# print the snoid for validation
echo $sno

#file with snoglobe's predictions for current snoRNA
snofile=$SLURM_TMPDIR/snoglobe_formatted/${sno}_uniq_regions.bed; 

gtf=$SLURM_TMPDIR/hg38_Ensembl_V101_Scottlab_2020.tsv # gtf file in tsv format
genome=$SLURM_TMPDIR/chrNameLength.txt # chromosome length file
tmp_sno_path=$SLURM_TMPDIR/snoglobe_formatted # sno path in $SLURM_TMPDIR/

# output path in $SLURM_TMPDIR
tmp_outpath=$SLURM_TMPDIR/sno_sno_overlaps_p_vals/
# create ouptut path
mkdir -p $tmp_outpath

# ouptut file
tmp_outfile=$tmp_outpath/${sno}_sno_overlap.tsv
# create output file if needed
touch $tmp_outfile

# go to $SLURM_TMPDIR/
cd $SLURM_TMPDIR/

# run the python script! you should change script path
# $SLURM_TMPDIR is a temporary directory to write temp files on the compute node
python3 compute_overlaps.py $snofile $tmp_sno_path $gtf $genome $tmp_outfile $SLURM_CPUS_PER_TASK $SLURM_TMPDIR sno;

# copy output files in $SLURM_TMPDIR to original output directory
outpath=$SCRATCH/network_analysis/network_output/sno_sno_overlaps_p_vals/
mkdir -p $outpath
cp $tmp_outfile $outpath/

echo 'Done'