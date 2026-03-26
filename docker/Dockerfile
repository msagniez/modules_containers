# Use a specific R version base image from the Rocker project
FROM rocker/r-ver:4.5.0 

RUN apt-get update && apt-get install -y \
    cmake \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libwebp-dev \
    libglpk-dev \
    pandoc \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Install R packages layer by layer to respect dependency order
RUN Rscript -e "install.packages(c('BiocManager','devtools','remotes','dplyr','stringr','umap','ggrepel'), repos = 'https://cran.rstudio.com/')"
RUN Rscript -e "BiocManager::install(version = '3.22')"
RUN Rscript -e "BiocManager::install(c('SummarizedExperiment','DESeq2','SingleR'))"
RUN Rscript -e "pak::pak('JinmiaoChenLab/Rphenograph')"
RUN Rscript -e "pak::pak('gu-lab20/MD-ALL')"

# Copy the prediction script into the container
COPY scripts/MDALL_predict.R /opt/scripts/MDALL_predict.R
RUN chmod +x /opt/scripts/MDALL_predict.R

# Set the working directory
WORKDIR /opt