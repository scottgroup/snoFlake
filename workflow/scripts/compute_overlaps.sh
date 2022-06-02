#!/bin/bash

# Run compute_overlaps.py and access $SLURM_TMPDIR for I/O operations

# copy all input data to $SLURM_TMPDIR
cp -r $2 $SLURM_TMPDIR/
cp $1 $3 $4 scripts/compute_overlaps.py $SLURM_TMPDIR/

# initialize input variables
in_file=$SLURM_TMPDIR/$(basename "$1")
target_path=$SLURM_TMPDIR/$(basename "$2")
gtf=$SLURM_TMPDIR/$(basename "$3")
genome=$SLURM_TMPDIR/$(basename "$4")
tmp_outfile=$SLURM_TMPDIR/$(basename "$5")
outfile=$5
ovlp_type=$6

cd $SLURM_TMPDIR &&
python3 compute_overlaps.py $in_file $target_path $gtf $genome $tmp_outfile $SLURM_CPUS_PER_TASK $SLURM_TMPDIR $ovlp_type

# copy output file in $SLURM_TMPDIR to original output directory
cp $tmp_outfile $outfile