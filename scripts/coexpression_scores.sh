#!/bin/bash
#PBS -l walltime=6:00:00,nodes=1:ppn=1,mem=60gb
#PBS -o /home/hirschc1/della028/projects/stem-rust_networks/analysis/msi_dump
#PBS -e /home/hirschc1/della028/projects/stem-rust_networks/analysis/msi_dump
#PBS -V
#PBS -N coexpression_scores_${NAME}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# start virtual environment
source activate camoco

# go to project folder
cd ~/projects/stem-rust_networks

# get number of clusters that contains S genes orthologs
clusters=$(sed 1d analysis/clusters/s_gene_orthologs_clusters.${NAME}.txt | cut -f 6 | sort -n | uniq | tr "\n" "," | sed '$ s/.$//')
# get coexpression zscores (edges) between genes of these clusters
~/.conda/envs/camoco/bin/python scripts/coexpression_scores.py ${NAME} ${clusters} analysis/clusters/orthologs/coexpr_scores_all_genes.${NAME}

# stop virtual environment
source deactivate
