#!/usr/bin/env python3

"""
Adapted from https://bitbucket.org/snakemake/snakemake/issues/28/clustering-jobs-with-snakemake
Launch with :
snakemake -j 999 --use-conda --immediate-submit --notemp --cluster-config cluster.json --cluster 'python3 scripts/slurmSubmit.py {dependencies}'
"""

import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

cmdline = "sbatch "

# direct slurm output and error file to a specific directory
slurm_outdir = '../network_output/slurmout'
if not os.path.exists(slurm_outdir):
    os.makedirs(slurm_outdir)
    
cmdline += " --output=../network_output/slurmout/%x_%j.out --error=../network_output/slurmout/%x_%j.err "

for param, val in job_properties['cluster'].items():
    cmdline += "--{param} {val} ".format(param=param, val=val)

# Set up dependencies
dependencies = set(sys.argv[1:-1])
if dependencies:
    cmdline += " --dependency=afterok:{} ".format(":".join(dependencies))

# Adding the actual job
cmdline += jobscript

# remove the leading and trailing white space for the submitted jobid
cmdline += r" | awk '{print substr($NF, 0, length($NF))}'"

sys.stdout.write(cmdline)

os.system(cmdline)