#!/bin/bash

# Filter physically binding STRING interactions

# path variables
string=$SCRATCH/p_3/STRING.9606.protein.actions.v11.0.txt.gz
out=$SCRATCH/p_3/output/rbp_rbp_STRING.tsv
rbp=$SCRATCH/p_3/rbp_protein_id.tsv

# filter by interaction type and combined score
gunzip -c $string | awk -F '\t' '($3 == "binding") && ($7 > 900) {print}' | cut -f 1,2 | uniq > $SCRATCH/p_3/STRING_temp.tsv

# call python script
python3 $SCRATCH/p_3/scripts/filter_STRING.py $rbp $SCRATCH/p_3/STRING_temp.tsv $out

# delete temp file
rm $SCRATCH/p_3/STRING_temp.tsv 