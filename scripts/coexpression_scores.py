#!/usr/bin/python3

import camoco as co
import pandas as pd
import numpy as np
import argparse as ap

# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script extracts co-expression scores between genes of a
             Camoco cluster (i.e. subnetwork).''')
# add positional arguments
parser.add_argument("network_name", type=str,
                    help="name of Camoco network")
parser.add_argument("cluster_list", type=str,
                    help="comma-separated list of Camoco clusters")
parser.add_argument("output_prefix", type=str,
                    help="prefix of output name")
# pass arguments into variables
args = parser.parse_args()
network_name = args.network_name
cluster_list = args.cluster_list
output_prefix = args.output_prefix


# transform network into a COB object
COB = co.COB(network_name)

# get all gene ids in camoco network
genes_network = COB.genes()

# transform clusters into integers
cluster_list = [int(cluster) for cluster in cluster_list.split(",")]

for cluster in cluster_list:

    print("Retrieving cluster", str(cluster))
    # get all gene names in cluster
    genes_cluster = COB.clusters.loc[COB.clusters["cluster"] == cluster]
    # names are actually row names
    genes_cluster = list(genes_cluster.index)

    # the gene ids from networks are actually a class with different attributes
    # (e.g. name, chromosome, start/end position, etc.), but the gene ids in the
    # cluster shows only the name of the gene

    # get full gene ids for all genes in cluster
    genes_cluster = [gene for gene in genes_network if gene.name in genes_cluster]

    # get data frame with co-expression scores
    scores_df = COB.subnetwork(genes_cluster, sig_only = False, names_as_cols = True)
    # keep only important columns
    scores_df = scores_df[["gene_a", "gene_b", "score"]]

    # export data frame data to tsv
    output_name = output_prefix + ".cluster_" + str(cluster) + ".txt"
    scores_df.to_csv(output_name, sep = "\t", index = False)
