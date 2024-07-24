# **snoFlake**: snoRNA-RBP Functional Interaction Network Models

A pipeline to gather snoRNA and RBP interactions to perform a sno-RBP interaction network analysis.

**Author:** [Kristina Sungeun Song](mailto:kristina.song@usherbrooke.ca)


## Running the Snakemake workflow
For a dry-run of this Snakemake workflow, simply run the following code from `CRSSANT/`.
```
snakemake -n
```
To run this Snakemake workflow, simply run the following code from `CRSSANT/`.
```
snakemake --profile profile_slurm

OR

snakemake --profile profile_local
```

**Types of snoRNA-RBP Interactions:**
ADD FIGURE

## Step 0: Download and prepare environment
snoFlake is supported on Linux (tested on Ubuntu).

## Installing Snakemake
Please follow the [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install the Snakemake workflow management tool. We recommend using `Conda/Mamba` to install Snakemake.

This Snakemake workflow has been tested with `v7.32.4`.

## Step 1: Get input datasets from various sources

## Step 2: Compute snoRNA and RBP interactions

## Step 3: Build snoFlake on Cytoscape