library(data.table)


#### functions ----

usage <- function() {
  cat("
description: get the number of Camoco cluster containing a susceptibility gene

usage: Rscript S_genes_clusters.R [file_clusters] [file_orthologs] [species] [output_name]

positional arguments:
  file_clusters           file with gene ids and cluster information from Camoco
  file_orthologs          file with known susceptibility (S) genes and their orthologs
  species                 species to associate ortholog ids with cluster info ('brachy' or 'wheat')
  output_name             name of output file

optional argument:
  --help                  show this helpful message

"
  )
}



#### command line options ----

# get arguments
args <- commandArgs(trailingOnly = TRUE)

# assert to have the correct optional arguments
if ("--help" %in% args) usage() & q(save = "no")

if (length(args) != 4) stop(usage(), "missing positional argument(s)")



#### assign cluster number to orthologs ----

file_clusters <- args[1]
file_ortho <- args[2]
species <- args[3]
output_name <- args[4]
# file_clusters <- "analysis/clusters/network_clusters.Brachy_stem_rust_Meesh.csv"
# file_ortho <- "data/s_gene_orthologs.corrected-IDs.txt"
# species <- "brachy"
# output_name <- "analysis/clusters/s_gene_orthologs_clusters.Brachy_stem_rust_Meesh.txt"

# load info
info_clusters <- fread(file_clusters, header = TRUE, data.table = FALSE)
info_ortho <- fread(file_ortho, header = TRUE, data.table = FALSE)

# filter orthologs info to have only species of interest
species_col <- grep(species, colnames(info_ortho))
info_ortho <- info_ortho[, c(1:4, species_col)]
# remove NAs
info_ortho <- info_ortho[which(!is.na(info_ortho[, 5])), ]
# make gene ids upper case for compatibility with Camoco
info_ortho[, 5] <- toupper(info_ortho[, 5])

# keep clusters that have the species orthologs
info_clusters <- info_clusters[which(info_clusters[, 1] %in% info_ortho[, 5]), ]

# create new column to keep cluster number
info_ortho$cluster <- NA
# assign cluster number to orthologs
for (row in 1:NROW(info_clusters)) {

  gene <- info_clusters[row, 1]
  cluster <- info_clusters[row, 2]

  info_ortho[which(info_ortho[, 5] == gene), "cluster"] <- cluster

}

# order by cluster
info_ortho <- info_ortho[order(info_ortho$cluster, info_ortho$orthogroup, info_ortho$S_gene), ]

# remove genes not in a cluster
info_ortho <- info_ortho[which(!is.na(info_ortho$cluster)), ]

fwrite(info_ortho, output_name, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
