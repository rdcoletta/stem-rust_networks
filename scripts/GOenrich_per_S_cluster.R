library(data.table)


#### functions ----

usage <- function() {
  cat("
description: get enriched GO terms in clusters with susceptibility gene orthologs.

usage: Rscript GOenrich_per_S_cluster.R [file_orthologs] [file_go] [network] [output_name]

positional arguments:
  file_orthologs          file with known susceptibility (S) genes and their orthologs
  file_go                 file with go enrichment per cluster
  network                 name of the network
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



#### for each cluster with orthologs, what GOs are enriched? ----

file_ortho <- args[1]
file_go <- args[2]
network <- args[3]
output_name <- args[4]
# file_ortho <- "analysis/clusters/s_gene_orthologs_clusters.SR_brachy_1.txt"
# file_go <- "analysis/go_enrichment/go_slim.SR_brachy_1.all_clusters.txt"
# network <- "SR_brachy_1"
# output_name <- "analysis/go_enrichment/GOenrich_S_clusters.SR_brachy_1.txt"



# load matrix with orthologs information
info_ortho <- fread(file_ortho, header = TRUE, data.table = FALSE)
# open file with enriched go terms
info_go <- fread(file_go, header = TRUE, data.table = FALSE, fill = TRUE)

# get clusters with susceptibility genes
clusters <- unique(info_ortho[, "cluster"])

# subset by clusters with known S genes
info_go_clusters <- subset(info_go, Cluster %in% clusters)
# keep only important GO information
info_go_clusters <- apply(info_go_clusters[, 1:2], MARGIN = 1, function(go) {

  # go descriptions can be very long and not very friendly to be parsed programatically
  # that's why i'm adding extra steps to retrieve go term, name and description

  go_term <- unlist(strsplit(go[2], split = "Term: "))
  go_term <- unlist(strsplit(go_term, split = ","))[1]

  go_name <- unlist(strsplit(go[2], split = "Name: "))
  go_name <- unlist(strsplit(rev(go_name)[1], split = ","))[1]

  go_desc <- unlist(strsplit(go[2], split = "\""))[2]

  return(data.frame(cluster = go[1], go_term, go_name, go_desc))

})
info_go_clusters <- do.call(rbind, info_go_clusters)

# remove empty spaces on cluster column
info_go_clusters$cluster <- gsub(" ", "", info_go_clusters$cluster)

# add two new empty columns in data frame
info_go_clusters$network_genes <- NA
info_go_clusters$S_genes <- NA

for (cluster in clusters) {

  # get which genes from the network in are in the cluster
  network_genes <- unique(info_ortho[which(info_ortho[, "cluster"] == cluster), 5])
  network_genes <- paste0(network_genes, collapse = ",")
  info_go_clusters[which(info_go_clusters[, "cluster"] == cluster), "network_genes"] <- network_genes

  # get which S gene orthologs are in the cluster
  S_genes <- unique(info_ortho[which(info_ortho[, "cluster"] == cluster), 1])
  S_genes <- paste0(S_genes, collapse = ",")
  info_go_clusters[which(info_go_clusters[, "cluster"] == cluster), "S_genes"] <- S_genes

}

# add name of network
info_go_clusters <- data.frame(network = network, info_go_clusters)

# write output
fwrite(info_go_clusters, output_name, sep = "\t", quote = FALSE, na = NA, row.names = FALSE)
