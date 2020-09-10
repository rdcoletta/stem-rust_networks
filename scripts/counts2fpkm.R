library(data.table)
library(doParallel)


#### functions ----

usage <- function() {
  cat("
description: generate FPKM values from raw counts of HTseq

usage: Rscript counts2fpkm.R [raw_counts_file] [annotation_file] [output_name] [...]

positional arguments:
  raw_counts_file         input file with raw counts from HTseq
  annotation_file         gff3 file with annotation for each gene
  output_name             name of output file

optional argument:
  --help                  show this helpful message
  --cores=VALUE           number of cores to use for parallel computing (default = 1)

"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(as.integer(arg[2]))
  
}


#### command line options ----

# set default
cores <- 1

# get arguments
args <- commandArgs(trailingOnly = TRUE)

# assert to have the correct optional arguments
if ("--help" %in% args) usage() & q(save = "no")

if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {
  
  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--cores")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")
  
  # change default based on the argument provided
  if (any(grepl("--cores", opt_args_requested))) cores <- getArgValue(opt_args[grep("--cores", opt_args)])
  
}


#### raw counts to fpkm ----

raw_counts_file <- args[1]
annotation_file <- args[2]
output_name <- args[3]
# raw_counts_file <- "data/brachy_counts_raw.txt"
# annotation_file <- "data/BdistachyonBd21_3_537_v1.2.gene_exons.no-header.gff3"
# output_name <- "data/brachy_counts_fpkm.txt"
# raw_counts_file <- "data/wheat_counts_raw.txt"
# annotation_file <- "data/IWGSC_v1.1_HC_20170706.no-header.gff3"
# output_name <- "data/wheat_counts_fpkm.txt"
# cores <- 5

# load files
counts <- fread(raw_counts_file, header = TRUE, data.table = FALSE)
annotation <- fread(annotation_file, header = FALSE, data.table = FALSE)
colnames(annotation) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# remove extra htseq-related rows
counts <- counts[grep("^_", counts[, 1], perl = TRUE, invert = TRUE), ]

# note: total gene length should be exons + introns (not just exons), because we counted
# any read landing within the start/stop as counting for that particular gene

# get gene name and length
annotation <- annotation[annotation$type == "gene", ]
annotation$geneid <- lapply(annotation[, "attributes"], function(att) {
  
  id <- unlist(strsplit(att, split = ";"))
  id <- id[grep("ID=", id)]
  id <- unlist(strsplit(id, split = "="))[2]
  
  return(id)
  
})
annotation$geneid <- do.call(c, annotation$geneid)
annotation[, "length"] <- as.numeric(annotation[, "end"]) - as.numeric(annotation[, "start"])

# make sure that 'counts' and 'annotation' have the same genes
if (NROW(counts) == NROW(annotation) & all(counts$geneid %in% annotation$geneid) & all(annotation$geneid %in% counts$geneid)) {
  
  # make both datasets with the same gene order
  annotation <- annotation[order(match(annotation[, "geneid"], counts[, "geneid"])), ]

} else {
  
  stop("Expression data and annotation have different genes")
  
}

# calculate fpkm per sample
registerDoParallel(cores)
fpkm <- foreach(sample = 2:NCOL(counts), .combine = cbind) %dopar% {
  
  # number of reads mapped to all protein-coding genes
  reads_all <- as.numeric(sum(counts[, sample]))
  
  fpkm_sample <- c()
  for (gene in 1:NROW(counts)) {
    
    # convert raw count of gene into fpkm
    reads_gene <- as.numeric(counts[gene, sample])
    length_gene <- as.numeric(annotation[gene, "length"])
    fpkm_gene <- (reads_gene * 10^9) / (length_gene * reads_all)
    
    fpkm_sample <- c(fpkm_sample, fpkm_gene)
    
  }
  
  # return fpkm of all genes for a sample
  fpkm_sample
  
}
stopImplicitCluster()

# add gene ids and sample names
fpkm <- data.frame(fpkm, stringsAsFactors = FALSE)
fpkm <- data.frame(geneid = counts[, "geneid"], fpkm)
colnames(fpkm)[2:NCOL(fpkm)] <- colnames(counts)[2:NCOL(counts)]

# write final files
fwrite(fpkm, file = output_name, sep = "\t", na = NA, row.names = FALSE, quote = FALSE)
