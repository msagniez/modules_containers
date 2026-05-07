#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat ("Usage: Rscript AMLmapR_predict.R </full/path/to/input/quant_file.csv> </full/path/to/output/directory> \n")
  cat ("\nOutput: Prediction_results.tsv in output directory\n")
  quit(status = 1)
}


quant <- args[1]
outdir <- args[2]
curdir <- getwd()

#install Libraries
suppressPackageStartupMessages({
  library(caret)
  library(kernlab)
  library(AMLmapR)
  library(tools)
})


#Load matrix
quant_df <- read.table(quant, sep=",", header = T, row.names = 1)
df_matrix <- as.matrix(quant_df[,colnames(quant_df)])
storage.mode(df_matrix) <- "integer"

#Predict using AMLmapR
predictions <- predict_AML_clusters(df_matrix)

#Save results
out_name <- paste(outdir,"/","AMLmapR_predictions.csv",sep="",collapse=NULL)
write.table(predictions, out_name, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)