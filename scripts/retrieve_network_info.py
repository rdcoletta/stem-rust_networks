#!/usr/bin/python3

'''
created by Rafael Della Coletta
2019-09-05
'''

import camoco as co
import pandas as pd
import numpy as np
import argparse as ap

# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script extracts both gene expression and cluster information
             from a specific network of Camoco's database.''')
# add positional arguments
parser.add_argument("network_name", type=str,
                    help="name of Camoco network")
parser.add_argument("out_dir", type=str,
                    help="path to store output")
# pass arguments into variables
args = parser.parse_args()
network_name = args.network_name
out_dir = args.out_dir

# transform network into a COB object
COB = co.COB(network_name)
# assign expression and cluster data to a variable
expression = COB.expr()
clusters = COB.clusters
# export expression and cluster data to csv
expression.to_csv(out_dir + "/network_expr." + network_name + ".csv")
clusters.to_csv(out_dir + "/network_clusters." + network_name + ".csv")
