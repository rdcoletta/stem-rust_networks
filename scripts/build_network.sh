#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/stem-rust_networks/analysis/msi_dump
#PBS -e /home/hirschc1/della028/projects/stem-rust_networks/analysis/msi_dump
#PBS -V
#PBS -N build_network_${NAME}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# start virtual environment
source activate camoco

# go to project folder
cd ~/projects/stem-rust_networks
# build network
camoco build-cob --rawtype RNASEQ --sep , ${OPT} ${IN} ${NAME} ${DESC} ${REF}
# get summary stats
camoco health --out ${HEALTH}/${NAME} --refgen ${REF} ${NAME}

# stop virtual environment
source deactivate
