# Use a specific R version base image from the Rocker project
FROM rocker/r-ver:4.2.0

# Install system dependencies required by devtools and R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN Rscript -e "install.packages(c('devtools', 'tools'), repos='https://cran.rstudio.com/')"
RUN Rscript -e "devtools::install_github('ThomasBeder/ALLCatchR')"
RUN Rscript -e "devtools::install_github('ThomasBeder/ALLCatchR2')"
RUN Rscript -e "devtools::install_github('ThomasBeder/ALLCatchRbcrabl1')"

# Copy the prediction script into the container
COPY scripts/ALLcatchR_predict.R /scripts/ALLcatchR_predict.R
RUN chmod +x /scripts/ALLcatchR_predict.R

# Set the working directory
WORKDIR /opt