A pipeline to gather snoRNA and RBP interactions to perform a network analysis.

Dependencies:
snakemake
bedtools
pandas

Author: Kristina Sungeun Song kristina.song@usherbrooke.ca

inputs:
- RBP data downloaded from ENCODE
    - Target category: select 'RNA binding protein'
- snoGloBe predictions against all protein coding biotypes
- snoDB: all snoRNAs that have protein coding host genes
- list of RBPs and snoRNAs

Types of snoRNA-RBP Interactions:
- snoRNAS that are embedded in RBP host genes
- RBPs that directly bind to snoRNAs
- RBP-RBP interactions obtained from STRING
- snoRNAS that bind to RBP mRNA transcripts
- Overlaps of snoRNA and RBP target sites
- High-throughput RNA-RNA interactions
- snoGloBe predictions

- Overlap of target sites:
    - snoRNA-snoRNA
    - RBP-RBP
    - snoRNA-RBP
- Direct binding interactions:
- Interactions from external database:
    - STRING