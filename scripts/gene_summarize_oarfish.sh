#!/usr/bin/env Rscript

#gene_summarize.R - Summarize oarfish transcript quantifications to gene level

#Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript gene_summarize.R [--gene-names] [--ignore-tx-version] <gtf_file> <quant_file1> [quant_file2 ...]\n")
  cat("   or: Rscript gene_summarize.R [--gene-names] [--ignore-tx-version] --tx2gene <tx2gene.txt> <quant_file1> [quant_file2 ...]\n")
  cat("\nOptions:\n")
  cat("  --gene-names         Use gene names instead of gene IDs in output\n")
  cat("  --ignore-tx-version  Strip version suffixes from transcript IDs\n")
  cat("                       Handles both directions:\n")
  cat("                         - quant has versions, tx2gene does not: stripped via ignoreTxVersion in tximport\n")
  cat("                         - tx2gene has versions, quant does not: stripped from tx2gene directly\n")
  cat("\nOutput: gene_counts.tsv and gene_tpm.tsv in current directory\n")
  quit(status = 1)
}

#Check for --gene-names option
use_gene_names <- FALSE
if ("--gene-names" %in% args) {
  use_gene_names <- TRUE
  args <- args[args != "--gene-names"]
}

#Check for --ignore-tx-version option
ignore_tx_version <- FALSE
if ("--ignore-tx-version" %in% args) {
  ignore_tx_version <- TRUE
  args <- args[args != "--ignore-tx-version"]
}

cat(sprintf("ignore_tx_version = %s\n", ignore_tx_version))

#Load required libraries
suppressPackageStartupMessages({
  library(tximport)
  library(txdbmaker)
  library(GenomicFeatures)
  library(AnnotationDbi)
})

#Check if using tx2gene file or gtf
if (args[1] == "--tx2gene") {
  tx2gene <- read.table(args[2], header = FALSE, col.names = c("TXNAME","GENEID"))
  quant_files <- args[3:length(args)]
  
  if (use_gene_names) {
    cat("Warning: --gene-names option requires GTF file, not tx2gene file. Using gene IDs.\n")
    use_gene_names <- FALSE
  }
} else {
  gtf_file <- args[1]
  quant_files <- args[2:length(args)]
  
  cat("Creating transcript database from GTF...\n")
  suppressWarnings({
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_file)
  })
  
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
  tx2gene <- tx2gene[, c("TXNAME", "GENEID")]
  
  if (use_gene_names) {
    cat("Extracting gene names from GTF...\n")
    gtf_lines <- readLines(gtf_file)
    gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
    
    gene_id_pattern   <- 'gene_id "([^"]+)"'
    gene_name_pattern <- 'gene_name "([^"]+)"'
    
    gene_ids   <- regmatches(gtf_lines, regexpr(gene_id_pattern,   gtf_lines))
    gene_ids   <- gsub('gene_id "|"', '', gene_ids)
    
    gene_names <- regmatches(gtf_lines, regexpr(gene_name_pattern, gtf_lines))
    gene_names <- gsub('gene_name "|"', '', gene_names)
    
    id2name <- data.frame(GENEID = gene_ids, GENENAME = gene_names, stringsAsFactors = FALSE)
    id2name <- unique(id2name[id2name$GENEID != "" & id2name$GENENAME != "", ])
    
    cat("Found", nrow(id2name), "gene ID to name mappings\n")
  }
}

#Detect whether quant files have version suffixes on transcript IDs
first_quant <- read.table(quant_files[1], header = TRUE, nrow = 100)
quant_ids <- as.character(first_quant[[1]])
quant_has_version <- any(grepl("\\.[0-9]+$", quant_ids))
cat(sprintf("quant_has_version = %s\n", quant_has_version))

#Strip GeneIDs versions when needed
if (ignore_tx_version) {
  # Strip versions from both sides unconditionally
  cat("--ignore-tx-version set: stripping version suffixes from both tx2gene and quant file IDs...\n")
  before <- nrow(tx2gene)
  tx2gene$TXNAME <- sub("\\.[0-9]+$", "", tx2gene$TXNAME)
  tx2gene$GENEID <- sub("\\.[0-9]+$", "", tx2gene$GENEID)
  tx2gene <- unique(tx2gene)
  cat(sprintf("tx2gene: %d -> %d rows after deduplication\n", before, nrow(tx2gene)))
  tximport_ignore_version <- TRUE
  
} else if (quant_has_version) {
  # Quant IDs have versions: keep tx2gene as-is to match
  cat("Quant IDs have versions — keeping tx2gene as-is\n")
  tximport_ignore_version <- FALSE
  
} else {
  # Quant IDs have no versions: strip versions from tx2gene to match
  cat("Quant IDs have no versions — stripping version suffixes from tx2gene to match...\n")
  before <- nrow(tx2gene)
  tx2gene$TXNAME <- sub("\\.[0-9]+$", "", tx2gene$TXNAME)
  tx2gene$GENEID <- sub("\\.[0-9]+$", "", tx2gene$GENEID)
  tx2gene <- unique(tx2gene)
  cat(sprintf("tx2gene: %d -> %d rows after deduplication\n", before, nrow(tx2gene)))
  tximport_ignore_version <- FALSE
}

#Set sample names from file names
names(quant_files) <- gsub(".*/", "", gsub("\\..*$", "", quant_files))
cat("Importing", length(quant_files), "quantification file(s)...\n")

# ignoreTxVersion = TRUE only when flag set AND versions are on the quant file side
txi <- tximport(quant_files, type = "oarfish", tx2gene = tx2gene, ignoreTxVersion = ignore_tx_version)

#Write gene-level counts (rounded to integers)
cat("Writing gene-level counts to gene_counts.tsv...\n")
gene_counts <- as.data.frame(round(txi$counts))
gene_counts$gene_id <- rownames(gene_counts)

if (use_gene_names) {
  matched_names <- id2name$GENENAME[match(gene_counts$gene_id, id2name$GENEID)]
  gene_counts$gene_id <- ifelse(is.na(matched_names), gene_counts$gene_id, matched_names)
  cat("Converted", sum(!is.na(matched_names)), "gene IDs to names\n")
}
gene_counts <- gene_counts[, c("gene_id", setdiff(names(gene_counts), "gene_id"))]
colnames(gene_counts)[1] <- if(use_gene_names) "gene_name" else "gene_id"
write.table(gene_counts, "gene_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Write gene-level TPMs (not rounded, keep as floating point)
cat("Writing gene-level TPM to gene_tpm.tsv...\n")
gene_tpm <- as.data.frame(txi$abundance)
gene_tpm$gene_id <- rownames(gene_tpm)

if (use_gene_names) {
  matched_names <- id2name$GENENAME[match(gene_tpm$gene_id, id2name$GENEID)]
  gene_tpm$gene_id <- ifelse(is.na(matched_names), gene_tpm$gene_id, matched_names)
}
gene_tpm <- gene_tpm[, c("gene_id", setdiff(names(gene_tpm), "gene_id"))]
colnames(gene_tpm)[1] <- if(use_gene_names) "gene_name" else "gene_id"
write.table(gene_tpm, "gene_tpm.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done! Output files: gene_counts.tsv, gene_tpm.tsv\n")