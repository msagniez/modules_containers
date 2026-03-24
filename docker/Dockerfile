# Use a specific R version base image from the Rocker project
FROM rocker/r-ver:4.2.0 

# Set the working directory inside the container
WORKDIR /opt

# Copy your local R script and other files into the container's working directory
COPY mdall /opt/mdall
WORKDIR /opt/mdall

# Install any R packages required by your script
RUN install.packages(c("devtools", "BiocManager", "dplyr", "stringr", "Seurat", "ggplot2", "ggrepel", "cowplot", "umap"), repos='cran.rstudio.com')
RUN BiocManager::install(c("DESeq2", "SingleR", "SummarizedExperiment"))
RUN devtools::install_github("JinmiaoChenLab/Rphenograph")
RUN devtools::install_github("gu-lab20/MD-ALL")

ENV PATH="/opt/mdall/bin:${PATH}"