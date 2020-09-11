library(data.table)
library(ggplot2)
library(ggrepel)


#### functions ----

usage <- function() {
  cat("
description: plot PCA of gene expression data

usage: Rscript pca_expr_data.R [expression_data] [outfile_name]

positional arguments:
  expression_data         input file with expression data
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

if (length(args) < 2) stop(usage(), "missing positional argument(s)")




#### plot pca ----

expr_file <- args[1]
outfile_name <- args[2]
# expr_file <- "data/brachy_counts_raw.txt"
# outfile_name <- "data/pca_brachy_raw.png"
# expr_file <- "data/brachy_counts_fpkm.txt"
# outfile_name <- "data/pca_brachy_fpkm.png"

# load data
expr <- fread(expr_file, header = TRUE, data.table = FALSE)

# in PCA, we need samples as rows and genes as columns
expr <- data.frame(t(expr), stringsAsFactors = FALSE)

# adjust column names
colnames(expr) <- expr[1, ]
expr <- expr[-1, ]

# convert df to numeric after transposition
expr <- data.frame(apply(expr, MARGIN = c(1,2), as.numeric), stringsAsFactors = FALSE)

# get pcs
pca <- prcomp(expr)
pca <- as.data.frame(pca$x)

# add sample ids for ploting
sample_ids <- lapply(rownames(pca), function(name) unlist(strsplit(name, split = "_")))
sample_ids <- data.frame(do.call(rbind, sample_ids), stringsAsFactors = FALSE)
colnames(sample_ids) <- c("genotype", "dpi", "group", "rep")

sample_ids <- cbind(sample_ids, dpi_rep = paste(sample_ids[, 2], sample_ids[, 4], sep = "_"))

pca <- cbind(pca, sample_ids)

# plot
plot <- ggplot(data = pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = group, shape = rep), size = 3) +
  geom_label_repel(aes(label = dpi), label.size = 0) +
  theme(panel.background = element_blank()) +
  theme(legend.position="top", legend.title=element_blank(), legend.key=element_blank()) +
  theme(axis.line = element_line(colour = "black", size=.25)) +
  scale_color_manual(values = c("#1b9e77", "#d95f02"))

ggsave(filename = outfile_name, plot = plot, device = "png")
