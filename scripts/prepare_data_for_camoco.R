library(data.table)
library(ggplot2)


#### functions ----

usage <- function() {
  cat("
description: adjust expression dataset to match the matrix format required by Camoco
             for building the networks

usage: Rscript prepare_data_for_camoco.R [expression_data] [output_name] [...]

positional arguments:
  expression_data         input file with expression data
  output_name             name of output file (user needs to include '.csv' extension if filtering
                          or '.png' if plotting distribution)

optional argument:
  --help                  show this helpful message
  --plot-cv               only plot distribution of genes CV (no filtering happens)
  --filter-cv=VALUE       filter out genes with CV < 0.1 (default)


"
  )
}

getArgValue <- function(arg) {

  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(as.numeric(arg[2]))

}


#### command line options ----

# set default
filter_cv <- NULL
plot_cv <- FALSE

# get arguments
args <- commandArgs(trailingOnly = TRUE)

# assert to have the correct optional arguments
if ("--help" %in% args) usage() & q(save = "no")

if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {

  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--filter-cv", "--plot-cv")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")

  # change default based on the argument provided
  if (any(grepl("--filter-cv", opt_args_requested))) filter_cv <- getArgValue(opt_args[grep("--filter-cv", opt_args)])
  if (any(grepl("--plot-cv", opt_args_requested))) plot_cv <- getArgValue(opt_args[grep("--plot-cv", opt_args)])

}



#### calculate CV ----

expr_file <- args[1]
output_name <- args[2]
# expr_file <- "data/brachy_counts_fpkm.txt"
# output_name <- "data/brachy_counts_fpkm.csv"
# output_name <- "data/brachy_counts_fpkm.cv0-1.csv"

# load data
expr <- fread(expr_file, header = TRUE, data.table = FALSE)

# format data
rownames(expr) <- expr[, 1]
expr <- expr[, -1]

# calculate cv
genes_cv <- apply(expr,  MARGIN = 1, function(gene) sd(gene) / mean(gene))


if (plot_cv) {

  # transform NA into zero so those genes also appear in the plot
  genes_cv[which(is.na(genes_cv))] <- 0
  
  cat(sum(genes_cv == 0), " genes (", round((sum(genes_cv == 0) / length(genes_cv)) * 100, digits = 2),
      "%) not expressed in any sample\n", sep = "")
  for (threshold in seq(0.1, 1.0, 0.1)) {
    cat(sum(genes_cv < threshold), " genes (", round((sum(genes_cv < threshold) / length(genes_cv)) * 100, digits = 2),
        "%) with CV < ", threshold, "\n", sep = "")
  }
  
  # plot distribution
  plot <- ggplot(data.frame(cv = genes_cv), aes(x = cv)) +
    geom_histogram(binwidth = 0.1) +
    scale_x_continuous(breaks = seq(0, ceiling(max(genes_cv)), 0.5)) +
    coord_cartesian(xlim = c(0, ceiling(max(genes_cv))))

  ggsave(filename = output_name, plot = plot, device = "png")

  # exit after plotting
  q(save = "no")

} else if (!is.null(filter_cv)) {

  # filter data by cv
  genes_above_threshold <- which(genes_cv >= filter_cv)

  n_genes_removed <- NROW(expr) - length(genes_above_threshold)
  cat("Removed ", n_genes_removed, " genes (", round((n_genes_removed / NROW(expr) * 100), digits = 2),
      "%) with CV < ", filter_cv, "\n", sep = "")

  expr <- expr[genes_above_threshold, ]

}

# write output
fwrite(expr, output_name, quote = FALSE, sep = ",", row.names = TRUE)

# remove first comma from first row (header) for compatibility with Camoco
delete_first_comma <- paste0("sed -i '1 s/^,//' ", output_name)
system(delete_first_comma)
