#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat ("Usage: Rscript MDALL_predict.R </full/path/to/input/quant_files/directory> </full/path/to/output/directory>\n")
  cat ("\nOutput: Prediction_results.tsv in output directory\n")
  quit(status = 1)
}

quant <- args[1]
outdir <- args[2]
curdir <- getwd()

#install Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(Rphenograph)
  library(SummarizedExperiment)
  library(umap)
  library(ggrepel)
  library(Rphenograph)
  library(MDALL)
  library(tools)
})

#Initialize the list of samples to process
file_list <- list.files(path = quant, 
                        pattern = "\\.tsv$", 
                        full.names = TRUE) # full.names=TRUE returns the full path

# Initialize results matrix
df_results <- data.frame(
  id = character(),
  PhenoGraph_pred = character(),
  PhenoGraph_predScore = numeric(),
  PhenoGraph_predLabel = character(),
  svm_pred = character(),
  svm_predScore = numeric(),
  svm_predLabel = character(),
  stringsAsFactors = FALSE
)

# Keep track of failed samples
failed_samples <- c()

#Parse through the list of samples
for (file_path in file_list) {
  
  # Wrap entire processing in tryCatch
  tryCatch({
    
    df_count=read_input(file_path,delimiter = "\t",header = F)
    samp_name <- file_path_sans_ext(basename(file_path))
    print(paste("Processing file:", samp_name))
    
    #Normalize with reference data
    df_vst=get_vst_values(obj_in = obj_234_HTSeq, df_count = df_count)
    
    #Get normalized expression values for feature genes
    df_feateure_exp=get_geneExpression(df_vst = df_vst, genes = c("CDX2","CRLF2","NUTM1"))
    
    #Get gene expression box plot
    obj_boxplot=obj_merge(obj_in = obj_1821, df_in = df_vst, assay_name_in = "vst")
    boxplot_name <- paste(outdir,"/",samp_name,"_expression-boxplot.pdf",sep="",collapse=NULL)
    pdf(file = boxplot_name, width = 8, height = 6)
    print(draw_BoxPlot(obj_in = obj_boxplot, group.by = "diag_raw1", features = "CRLF2", highlightLevel = "TestSample", plot_title = "Box Plot of Expression"))
    dev.off()
    
    #Inputation
    df_vst_i=f_imputation(obj_ref = obj_234_HTSeq, df_in = df_vst)
    
    #Add testing sample to reference daTaset for subtype prediction
    obj_=obj_merge(obj_in = obj_1821, df_in = df_vst_i, assay_name_in = "vst")
    
    #Draw UMAP plot
    obj_=run_umap(obj_in = obj_, out_label = "umap", n_neighbors = 10, variable_n = c(1058), feature_panel = "keyFeatures")
    umap_name <- paste(outdir,"/",samp_name,"_umap.pdf",sep="",collapse=NULL)
    pdf(file = umap_name, width = 7, height = 6)
    print(draw_DimPlot(obj_, group.by = "diag_raw", reduction = "umap", highlightLevel = "TestSample"))
    dev.off()
    
    #Run PhenoGraph clustering and SVM prediction
    df_out_phenograph=get_PhenoGraphPreds(obj_in = obj_, feature_panel = "keyFeatures", SampleLevel = "TestSample",
                                          neighbor_k = 10,
                                          variable_n_list = c(100,1000,1058)
    )
    
    df_out_svm=get_SVMPreds(models_svm, df_in = df_vst_i)
    
    df_pred=bind_rows(df_out_phenograph, df_out_svm) %>% mutate(N=sprintf("%04d",featureN))
    
    heatmap_name <- paste(outdir,"/",samp_name,"_prediction-heatmap.pdf",sep="",collapse=NULL)
    pdf(file = heatmap_name, width = 8, height = 6)
    print(gg_tilePlot(df_in = df_pred, x = "N", y = "method", var_col = "pred", x_tick_label_var = "featureN", title = "Prediction Heatmap"))
    dev.off()
    
    #Get GEP prediction results
    df_gep_pred=get_GEP_pred(
      count_matrix = df_count,
      featureN_PG = c(100,200,300,400,500,600,700,800,900,1000,1058)
    )
    
    # Save Sample results
    new_row <- data.frame(
      id = samp_name,
      PhenoGraph_pred = df_gep_pred$phenoGraph_pred,
      PhenoGraph_predScore = df_gep_pred$phenoGraph_predScore,
      PhenoGraph_predLabel = df_gep_pred$phenoGraph_predLabel,
      svm_pred = df_gep_pred$svm_pred,
      svm_predScore = df_gep_pred$svm_predScore,
      svm_predLabel = df_gep_pred$svm_predLabel,
      stringsAsFactors = FALSE
    )
    
    # Add the new row to df_listing
    df_results <- rbind(df_results, new_row)
    
    print(paste("Successfully processed:", samp_name))
    
  }, error = function(e) {
    # Error handler - this runs if any error occurs
    samp_name <- file_path_sans_ext(basename(file_path))
    failed_samples <<- c(failed_samples, samp_name)
    cat("ERROR processing", samp_name, ":", conditionMessage(e), "\n")
    cat("Skipping to next sample...\n\n")
    
    # Close any open PDF devices in case error occurred during plotting
    if (dev.cur() > 1) {
      dev.off()
    }
  })
  
}

#Save results (without row names, with column names)
out_name <- paste(outdir,"/","MD-ALL_predictions.csv",sep="",collapse=NULL)
write.table(df_results, out_name, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

# Report on failed samples
if (length(failed_samples) > 0) {
  cat("\n=== PROCESSING SUMMARY ===\n")
  cat("Successfully processed:", nrow(df_results), "samples\n")
  cat("Failed samples (", length(failed_samples), "):\n", sep="")
  cat(paste("  -", failed_samples, collapse="\n"), "\n")
  
  # Optionally save failed samples list
  failed_name <- paste(outdir,"/","failed_samples.txt",sep="",collapse=NULL)
  writeLines(failed_samples, failed_name)
  cat("\nFailed samples list saved to:", failed_name, "\n")
} else {
  cat("\nAll samples processed successfully!\n")
}