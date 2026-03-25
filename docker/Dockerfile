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
RUN Rscript -e "install.packages('BiocManager', repos = 'https://cran.rstudio.com/')"
RUN Rscript -e "BiocManager::install(version = '3.22')"
RUN Rscript -e "BiocManager::install(c('tximport', 'txdbmaker', 'GenomicFeatures', 'AnnotationDbi'))"

# Copy the prediction script into the container
COPY scripts/gene_summarize.R /opt/scripts/gene_summarize.R
RUN chmod +x /opt/scripts/gene_summarize.R

# Set the working directory
WORKDIR /opt