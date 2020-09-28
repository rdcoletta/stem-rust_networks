#!/bin/bash
#PBS -l walltime=48:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/stem-rust_networks/analysis/msi_dump
#PBS -e /home/hirschc1/della028/projects/stem-rust_networks/analysis/msi_dump
#PBS -V
#PBS -N network_health
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# start virtual environment
source activate camoco

# go to project folder
cd ~/projects/stem-rust_networks
# generate health statistics
camoco health --out ${OUT} --refgen ${REF} --go ${GO} ${NAME}

# stop virtual environment
source deactivate
