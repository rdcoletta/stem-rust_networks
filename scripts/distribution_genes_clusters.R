library(data.table)
library(ggplot2)


#### functions ----

usage <- function() {
  cat("
description: plot distribution of genes per network cluster

usage: Rscript distribution_genes_clusters.R [cluster_info] [network_name] [outfile_name]

positional arguments:
  cluster_info            CSV file with information of which genes belong to which Camoco cluster
  network_name            name of network (to go in the title of the plot)
  outfile_name            name of output file

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

if (length(args) != 3) stop(usage(), "missing positional argument(s)")




#### plot pca ----

cluster_file <- args[1]
network_name <- args[2]
outfile_name <- args[3]
# cluster_file <- "analysis/clusters/network_clusters.SR_brachy_2.csv"
# network_name <- "SR_brachy_2"
# outfile_name <- "analysis/clusters/distribution_clusters.SR_brachy_2.png"

# load data
cluster_info <- fread(cluster_file, header = TRUE, data.table = FALSE)

# plot
plot <- ggplot(cluster_info, aes(x = cluster)) +
  geom_histogram(binwidth = 1) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 1000)) +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  labs(title = paste0("First 100 clusters of ", network_name),
       subtitle = "(y-axis cropped at 1,000 genes)",
       x = "Cluster",
       y = "Number of genes")
  
ggsave(filename = outfile_name, plot = plot, device = "png")
