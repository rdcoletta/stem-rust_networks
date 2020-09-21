#!/bin/bash
#PBS -l walltime=12:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/stem-rust_networks/analysis/msi_dump
#PBS -e /home/hirschc1/della028/projects/stem-rust_networks/analysis/msi_dump
#PBS -V
#PBS -N cluster_enrichment_${NAME}_${GO}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# if no optional variable requested, set it blank
if [[ -z ${CLUSTERS} ]]; then
  CLUSTERS=""
else
  CLUSTERS=$(echo ${CLUSTERS} | tr "-" ",")
  CLUSTERS="--clusters ${CLUSTERS}"
fi

# start virtual environment
source activate camoco

# go to project folder
cd ~/projects/stem-rust_networks
# perform enrichment analysis
~/.conda/envs/camoco/bin/python scripts/ClusterEnrichment.py --cob ${NAME} --go ${GO} ${CLUSTERS} > ${OUT}

# stop virtual environment
source deactivate
