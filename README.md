# **snoFlake**: snoRNA Functional Interaction Network Model

![snoFlake logo](images/snoFlake_logo.svg)

A Snakemake pipeline to gather snoRNA–RNA-binding protein (RBP) interactions from multiple sources and construct a snoRNA–RBP interaction network for visualization and analysis in Cytoscape.

**Author:** [Kristina Sungeun Song](mailto:kristina.song@usherbrooke.ca)  

---

## Table of Contents

- [Overview](#overview)
- [Installation and Environment Setup](#installation-and-environment-setup)
- [Running the Snakemake Workflow](#running-the-snakemake-workflow)
- [Step 0: Configure the Pipeline](#step-0-configure-the-pipeline)
- [Step 1: Get Input Datasets](#step-1-get-input-datasets)
- [Step 2: Compute snoRNA–RBP Interactions](#step-2-compute-snorna-rbp-interactions)
- [Step 3: Build snoFlake Network in Cytoscape](#step-3-build-snoflake-network-in-cytoscape)
- [Citation](#citation)

---

## Overview

snoFlake is a reproducible bioinformatics pipeline that integrates multiple types of snoRNA–RBP interaction evidence to build a comprehensive interaction network. The pipeline collects and processes data from eCLIP-seq experiments (ENCODE), computationally predicted snoRNA-RNA interactions (snoGloBe), and other interaction sources, computes interaction scores, and outputs network files ready for visualization in Cytoscape.

**Types of snoRNA–RBP Interactions integrated by snoFlake:**

> ⚠️ **TODO:** Add interaction type figure here.

---

## Installation and Environment Setup and Environment Setup

snoFlake is supported on **Linux** (tested on Ubuntu).

### 1. Clone the repository

```bash
git clone https://github.com/scottgroup/snoFlake.git
cd snoFlake
```

### 2. Install Snakemake

We recommend installing Snakemake via [Conda/Mamba](https://github.com/conda-forge/miniforge). Please follow the [official Snakemake installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.32.4
conda activate snakemake
```

> This workflow has been tested with Snakemake **v7.32.4**.

### 3. Install rule-level dependencies

Each rule in the workflow uses a dedicated Conda environment defined in `workflow/envs/`. Snakemake will automatically install these when you use the `--use-conda` flag (recommended).

Alternatively, required Python packages include:

- `pandas`
- `numpy`
- `pybedtools`
- `networkx`

---

## Running the Snakemake Workflow

First, perform a **dry run** to verify that the workflow is correctly configured and all input files are accessible:

```bash
# Run from the root snoFlake/ directory
snakemake -n
```

To execute the workflow, choose the appropriate profile for your system:

```bash
# On a SLURM cluster
snakemake --profile profile_slurm

# On a local machine
snakemake --profile profile_local
```

To enable automatic Conda environment management per rule:

```bash
snakemake --profile profile_slurm --use-conda
```

---

## Step 0: Configure the Pipeline

Before running the workflow, edit the main configuration file at `config/config.yaml`. Key parameters to set include:

| Parameter | Description |
|---|---|
| `output_dir` | Path to the desired output directory |
| `eclip_metadata` | Path to ENCODE eCLIP metadata table |

---

## Step 1: Get Input Datasets

This step downloads and prepares all input data required for computing snoRNA–RBP interactions.

### 1a. snoRNA annotations — snoDB

snoRNA annotations and snoRNA–protein interaction data are sourced from [snoDB 2.0](https://bioinfo-scottgroup.med.usherbrooke.ca/snoDB/), the Scott Lab's curated human snoRNA database. This includes snoRNA genomic coordinates, snoRNA class (C/D box or H/ACA box), host gene information, and pre-existing snoRNA–RBP interaction evidence extracted from ENCODE eCLIP data.

### 1b. RBP eCLIP data — ENCODE

eCLIP-seq data for RNA-binding proteins are downloaded from the [ENCODE portal](https://www.encodeproject.org/). The pipeline uses fold-enrichment over size-matched input controls to quantify RBP binding at snoRNA loci across cell lines (e.g., K562, HepG2). Datasets are filtered based on ENCODE quality metrics.

### 1c. Additional interaction sources

> ⚠️ **TODO:** snoRNA–target interaction predictions from snoGloBe

### Expected inputs summary

| Source | Type | Format |
|---|---|---|
| snoDB 2.0 | snoRNA annotations & known interactions | BED / TSV |
| ENCODE eCLIP | RBP–snoRNA binding evidence | BAM / BED / TSV |
| *(Additional sources)* | *(Interaction type)* | *(Format)* |

---

## Step 2: Compute snoRNA–RBP Interactions

This step processes the input datasets and computes interaction evidence scores for each snoRNA–RBP pair.

![snoFlake edges](images/snoFlake_edges.svg)

### 2a. Quantify RBP binding at snoRNA loci

eCLIP fold-enrichment values for each RBP are computed at snoRNA loci defined by snoDB coordinates. A snoRNA–RBP interaction is retained if the RBP shows significant enrichment (fold-enrichment ≥ threshold, as defined in `config.yaml`) over the size-matched input in at least one cell line.

### 2b. Integrate interaction evidence

Interaction evidence from multiple sources is merged and scored per snoRNA–RBP pair. Each edge in the final network is supported by one or more evidence types:

| Interaction Type | Evidence Source | Description |
|---|---|---|
| eCLIP-based binding | ENCODE eCLIP | RBP physically binds the snoRNA locus |
| *(Add additional types)* | *(Source)* | *(Description)* |

> ⚠️ **TODO:** Complete the interaction type table to match the figure referenced in the Overview section.

### 2c. Output interaction table

The pipeline produces a scored, tab-separated interaction table with one row per snoRNA–RBP pair, which is used as input for network construction in Step 3.

> ⚠️ **TODO:** Specify the exact output file name(s) and column definitions.

---

## Step 3: Build snoFlake Network in Cytoscape

The final interaction table produced by Step 2 is used to build and visualize the snoRNA–RBP interaction network in [Cytoscape](https://cytoscape.org/).

### 3a. Import the network

1. Open Cytoscape (tested with version ≥ **3.9**).
2. Go to **File → Import → Network from File** and select the interaction table output from Step 2.
3. Configure the import dialog:
   - **Source node column:** `snoRNA_id`
   - **Target node column:** `RBP_name`
   - **Interaction type column:** `interaction_type`
   - **Edge weight column:** `score` *(or equivalent)*

> ⚠️ **TODO:** Confirm exact column names from the Step 2 output file and update accordingly.

### 3b. Apply network style

A Cytoscape style file is provided in `resources/` to reproduce the network visualization used in the manuscript.

1. Go to **File → Import → Styles from File** and load `resources/snoFlake_style.xml` *(or equivalent)*.
2. Apply the imported style from the **Style** panel.

> ⚠️ **TODO:** Confirm the style file name and path in `resources/`.

### 3c. Network layout and analysis

We recommend the **Degree-sorted Circular Layout** for visualizing large snoRNA–RBP networks. Node and edge attributes exported from the pipeline can be used to color or size nodes by snoRNA class, RBP function, or interaction score.

---

## Citation

If you use snoFlake in your research, please cite:

> ⚠️ **TODO:** Replace with final citation once the manuscript is published in *Cell Genomics*.

```
[Author list]. snoFlake: [Full manuscript title].
[Journal]. [Year]. DOI: [DOI]
```

---

*For questions or issues, please open a GitHub Issue or contact [kristina.song@usherbrooke.ca](mailto:kristina.song@usherbrooke.ca).*
