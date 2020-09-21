# !/usr/bin/python3
"""This script is to be used in conjunction with the Camoco software.
It will perform go enrichments for every MCL cluster found within a network.

Written by Jean-Michel Michno (Meesh)
- Source: https://github.com/MeeshCompBio/Wheat_Brachy_Coex/blob/master/ClusterEnrichment.py

Slightly modified by Rafael Della Coletta
"""

# even though only a few packages are called in this script
# all of them are needed for Camoco software built in function
import sys
import os
import pandas as pd
import numpy as np
import camoco as co
import scipy as sp
import collections
import itertools
from itertools import chain
import getopt


def usage():
    print("\nThis is the usage function\n")
    print("python3 ClusterEnrichment.py -c <COB> -g <GOnt> -s <clusters>")
    print("-c or --cob       The name of your co-expression network")
    print("-g or --go        The name of your GO network object")
    print("-s or --clusters  The clusters to perform enrichment (comma-separated list)\n")

try:
    opts, args = getopt.getopt(sys.argv[1:], "c:g:s:h", ["cob=", "go=", "clusters=", "help"])
except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-c", "--cob"):
        cob = arg
    elif opt in ("-g", "--go"):
        go = arg
    elif opt in ("-s", "--clusters"):
        cluster_list = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"


# check if '--clusters' option was required
try: cluster_list
except NameError: cluster_list = None

# Set the network and GO object
cob = co.COB(cob)
go = co.GOnt(go)

TotGenes = cob.num_genes()

# change from np.ndarray to pd dataframe
Clusters = pd.DataFrame(cob.clusters)

if cluster_list is None:
    # if no specific clusters were required, get all cluster from network
    cluster_list = list(set(Clusters["cluster"]))
else:
    # otherwise, keep only requested clusters
    cluster_list = [int(cluster) for cluster in cluster_list.split(",")]

# Make a ordered dictionary to each key
# is a cluster and each value is the
# genes in that cluster
ClustDict = collections.OrderedDict()
for index, row in Clusters.iterrows():
    if row[0] in ClustDict.keys():
        ClustDict[row[0]].append(index)
    else:
        ClustDict[row[0]] = [index]

# write out the header
header = ("Cluster", "Term_info", "Loci", "Source", "pval(BF)",
          "Terms_tested", "Num_in_Common", "Source_term_size",
          "Target_term_size", "Num_of_terms", "Number_Universe")
sys.stdout.write("\t".join(header)+"\n")

# go through each cluster and set of genes
for key, value in ClustDict.items():
    ClusterID = key
    # perform enrichment only for cluster of interest
    if ClusterID in cluster_list:
        # convert gene list to IDs in network
        genes = (cob.refgen.from_ids(value))
        # run the built in enrichment funciton
        sigterms = go.enrichment(locus_list=genes, num_universe=TotGenes)
        if len(sigterms) == 0:
            final = [key, "No enrichment", "NA",
                     "NA", "NA", "NA",
                     "NA", "NA", "NA"]
            sys.stdout.write("\t".join(map(str, final))+"\n")
        else:
            # Pull out the relevant information
            for i in range(len(sigterms)):
                pval = str(sigterms[i].attrs["pval"])
                loci = (",".join(map(str, sigterms[i].loci.intersection(genes))))
                source = sigterms[i].attrs["hyper"]["source"]
                terms_tested = sigterms[i].attrs["hyper"]["terms_tested"]
                num_common = sigterms[i].attrs["hyper"]["num_common"]
                num_universe = sigterms[i].attrs["hyper"]["num_universe"]
                source_term_size = sigterms[i].attrs["hyper"]["source_term_size"]
                target_term_size = sigterms[i].attrs["hyper"]["target_term_size"]
                num_terms = sigterms[i].attrs["hyper"]["num_terms"]
                # num_sampled = sigterms[i].attrs["hyper"]["num_sampled"]
                final = itertools.chain([key], [sigterms[i]],
                                        [loci], [source],
                                        [pval], [terms_tested],
                                        [num_common], [source_term_size],
                                        [target_term_size],
                                        [num_terms], [num_universe])
                # write to STD OUT
                sys.stdout.write(("\t".join(map(str, final)))+"\n")
