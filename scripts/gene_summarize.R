#!/usr/bin/env Rscript

#gene_summarize.R - Summarize oarfish transcript quantifications to gene level

#Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat ("Usage: Rscript gene_summarize.R [--gene-names] <gtf_file> <quant_file1> [quant_file2 ...]\n")
  cat ("   or: Rscript gene_summarize.R [--gene-names] --tx2gene <tx2gene.txt> <quant_file1> [quant_file2 ...]\n")
  cat ("\nOptions:\n")
  cat ("  --gene-names    Use gene names instead of gene IDs in output\n")
  cat ("\nOutput: gene_counts.tsv and gene_tpm.tsv in current directory\n")
  quit(status = 1)
}

#Check for gene names option
use_gene_names <- FALSE
if (args[1] == "--gene-names") {
  use_gene_names <- TRUE
  args <- args[-1]
}

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
  
  # Detect if transcript IDs have version numbers
  sample_txnames <- head(tx2gene$TXNAME, 100)
  has_version <- any(grepl("\\.[0-9]+$", sample_txnames))
  ignore_tx_version <- !has_version
  
  if (has_version) {
    cat("Detected transcript IDs WITH version numbers (e.g., ENST00000516118.1)\n")
    cat("Setting ignoreTxVersion = FALSE\n")
  } else {
    cat("Detected transcript IDs WITHOUT version numbers (e.g., ENST00000516118)\n")
    cat("Setting ignoreTxVersion = TRUE\n")
  }
  
  if (use_gene_names) {
    cat("Warning: --gene-names option requires GTF file, not tx2gene file. Using gene IDs.\n")
    use_gene_names <- FALSE
  }
} else {
  gtf_file <- args[1]
  quant_files <- args[2:length(args)]
  
  cat("Creating transcript databases from GTF...\n")
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_file)
  
  # Get transcript to gene mapping
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
  tx2gene <- tx2gene[, c("TXNAME", "GENEID")]
  
  # Detect if transcript IDs have version numbers
  sample_txnames <- head(tx2gene$TXNAME, 100)
  has_version <- any(grepl("\\.[0-9]+$", sample_txnames))
  ignore_tx_version <- !has_version
  
  if (has_version) {
    cat("Detected transcript IDs WITH version numbers (e.g., ENST00000516118.1)\n")
    cat("Setting ignoreTxVersion = FALSE\n")
  } else {
    cat("Detected transcript IDs WITHOUT version numbers (e.g., ENST00000516118)\n")
    cat("Setting ignoreTxVersion = TRUE\n")
  }
  
  # If gene names requested, extract gene name mapping from GTF
  if (use_gene_names) {
    cat("Extracting gene names from GTF...\n")
    # Read GTF to get gene_id to gene_name mapping
    gtf_lines <- readLines(gtf_file)
    gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
    
    gene_id_pattern <- 'gene_id "([^"]+)"'
    gene_name_pattern <- 'gene_name "([^"]+)"'
    
    gene_ids <- regmatches(gtf_lines, regexpr(gene_id_pattern, gtf_lines))
    gene_ids <- gsub('gene_id "|"', '', gene_ids)
    
    gene_names <- regmatches(gtf_lines, regexpr(gene_name_pattern, gtf_lines))
    gene_names <- gsub('gene_name "|"', '', gene_names)
    
    id2name <- data.frame(GENEID = gene_ids, GENENAME = gene_names, stringsAsFactors = FALSE)
    id2name <- unique(id2name[id2name$GENEID != "" & id2name$GENENAME != "", ])
    
    cat("Found", nrow(id2name), "gene ID to name mappings\n")
  }
}

#Set sample names from files names
names(quant_files) <- gsub(".*/", "", gsub("\\..*$", "", quant_files))
cat("Importing", length(quant_files), "quantification file(s)...\n")

# Use ignoreTxVersion based on detected transcript ID format
txi <- tximport(quant_files, type = "oarfish", tx2gene = tx2gene, ignoreTxVersion = ignore_tx_version)

#Write gene-level counts (rounded to integers)
cat("Writing gene-level counts to gene_counts.tsv...\n")
gene_counts <- as.data.frame(round(txi$counts))  # Round to nearest integer
gene_counts$gene_id <- rownames(gene_counts)

# Replace gene IDs with gene names if requested
if (use_gene_names) {
  matched_names <- id2name$GENENAME[match(gene_counts$gene_id, id2name$GENEID)]
  # Keep gene ID if no name is found
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

# Replace gene IDs with gene names if requested
if (use_gene_names) {
  matched_names <- id2name$GENENAME[match(gene_tpm$gene_id, id2name$GENEID)]
  # Keep gene ID if no name is found
  gene_tpm$gene_id <- ifelse(is.na(matched_names), gene_tpm$gene_id, matched_names)
}
gene_tpm <- gene_tpm[, c("gene_id", setdiff(names(gene_tpm), "gene_id"))]
colnames(gene_tpm)[1] <- if(use_gene_names) "gene_name" else "gene_id"
write.table(gene_tpm, "gene_tpm.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("Done! Output files: gene_counts.tsv, gene_tpm.tsv\n")