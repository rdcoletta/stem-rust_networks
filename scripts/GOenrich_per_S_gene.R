library(data.table)


#### functions ----

usage <- function() {
  cat("
description: add enriched GO terms to susceptibility gene orthologs.

usage: Rscript GOenrich_per_S_gene.R [file_orthologs] [go_folder] [network] [output_name]

positional arguments:
  file_orthologs          file with known susceptibility (S) genes and their orthologs
  go_folder               folder with go enrichment files
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



#### get S genes with enriched GO terms ----

file_ortho <- args[1]
go_folder <- args[2]
network <- args[3]
output_name <- args[4]
# file_ortho <- "analysis/clusters/s_gene_orthologs_clusters.SR_brachy_1.txt"
# go_folder <- "analysis/go_enrichment/slim"
# network <- "SR_brachy_1"
# output_name <- "analysis/go_enrichment/GOenrich_S_genes.SR_brachy_1.txt"


# load matrix with orthologs information
info_ortho <- fread(file_ortho, header = TRUE, data.table = FALSE)

# create empty df to store resuts
df_ortho_go <- data.frame(stringsAsFactors = FALSE)

# get names of files containing GO enrichment per gene
go_files <- list.files(go_folder, pattern = network, full.names = TRUE)

for (file in go_files) {

  # get gene id -- which comes after network name
  gene <- unlist(strsplit(file, split = ".", fixed = TRUE))
  gene <- gene[grep(network, gene) + 1]

  # get onthology information for gene
  gene_ortho <- info_ortho[which(info_ortho[, 5] == gene), 1:3]
  colnames(gene_ortho) <- c("S_gene_name", "S_gene_species", "S_gene_id")
  # get cluster number of that gene
  cluster <- info_ortho[which(info_ortho[, 5] == gene), "cluster"]

  # open file with enriched go terms for that gene
  gene_go <- fread(file, header = TRUE, data.table = FALSE)

  for (go in 1:NROW(gene_go)) {

    go_info <- gene_go[go, "Term_info"]

    # go descriptions can be very long and not very friendly to be parsed programatically
    # that's why i'm adding extra steps to retrieve go term, name and description

    go_term <- unlist(strsplit(go_info, split = "Term: "))
    go_term <- unlist(strsplit(go_term, split = ","))[1]

    go_name <- unlist(strsplit(go_info, split = "Name: "))
    go_name <- unlist(strsplit(rev(go_name)[1], split = ","))[1]

    go_desc <- unlist(strsplit(go_info, split = "\""))[2]

    # add information to data frame
    df_ortho_go <- rbind(df_ortho_go, data.frame(network = network,
                                                 cluster = cluster,
                                                 gene = gene,
                                                 GO_term = go_term,
                                                 GO_name = go_name,
                                                 GO_desc = go_desc,
                                                 gene_ortho,
                                                 stringsAsFactors = FALSE))
  }
}

# sort df by cluster number (lower the number, more genes in the cluster)
df_ortho_go <- df_ortho_go[order(df_ortho_go$cluster), ]

# write output
fwrite(df_ortho_go, output_name, sep = "\t", quote = FALSE, na = NA, row.names = FALSE)
