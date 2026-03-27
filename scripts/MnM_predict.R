#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat ("Usage: Rscript MnM_predict.R <quant_file> </full/path/to/output/directory> \n")
  cat ("\nOutput: Prediction_results.tsv in output directory\n")
  quit(status = 1)
}

SAMPLE_NAME <- gsub("\\..*$", "", args[1])
quant <- args[1]
outdir <- args[2]
curdir <- getwd()

#install Libraries
suppressPackageStartupMessages({
  library("tidyverse")
  library("MnM")
})

#Load test set
countDataTest <- read.csv(quant, header = T, row.names = 1, sep = "\t") %>% as.matrix()

#Load pre-trained models
modelsMinority <- readRDS("./resources/createdModelsMinority.rds")
modelsMajority <- readRDS("./resources/createdModelsMajority.rds")

#Run MnM predictions
predictionsTestMinority <- newPredictionsMinority(createdModelsMinority = modelsMinority,
                                                  countDataNew = countDataTest,
                                                  outputDir = outdir
)
predictionsTestMajority <- newPredictionsMajority(createdModelsMajority = modelsMajority,
                                                  countDataNew = countDataTest,
                                                  outputDir = outdir
                                                  
) 


#Generate subtype-level predictions
predictionsMMTestList <- integrateMM(minority = predictionsTestMinority,
                                     majority = predictionsTestMajority,
                                     subtype = T)

predictionsMMTest <- predictionsMMTestList$predictionsMMFinal

finalpred <- paste(outdir,"/MnM_subtype-predictions.csv",sep="",collapse=NULL)
write.csv(predictionsMMTest, file = finalpred)

#Generate Type level predictions
predictionsMMTestList <- integrateMM(minority = predictionsTestMinority,
                                     majority = predictionsTestMajority,
                                     subtype = F)

predictionsMMTest <- predictionsMMTestList$predictionsMMFinal

finalpred <- paste(outdir,"/MnM_lineage-predictions.csv",sep="",collapse=NULL)
write.csv(predictionsMMTest, file = finalpred)