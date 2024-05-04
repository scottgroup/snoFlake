#!/bin/bash

# Run compute_overlaps.py and access $SLURM_TMPDIR for I/O operations

# store current working directory path
workdir=$(pwd)

# copy all input data to $SLURM_TMPDIR
cp -r $2 $SLURM_TMPDIR/
cp $1 $3 $4 workflow/scripts/compute_overlaps.py $SLURM_TMPDIR/

# initialize input variables
in_file=$(basename "$1")
target_path=$(basename "$2")
gtf=$(basename "$3")
genome=$(basename "$4")
tmp_outfile=$(basename "$5")
outfile=$workdir/$5

cd $SLURM_TMPDIR/ &&
python3 compute_overlaps.py $in_file $target_path $gtf $genome $tmp_outfile $SLURM_CPUS_PER_TASK $SLURM_TMPDIR $6

# copy output file in $SLURM_TMPDIR to original output directory
cp $tmp_outfile $outfile