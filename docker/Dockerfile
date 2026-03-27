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
    wget \
    && rm -rf /var/lib/apt/lists/*


WORKDIR /opt/models/
RUN wget "https://zenodo.org/records/14167359/files/createdModelsMajority.rds?download=1"
RUN wget "https://zenodo.org/records/14167359/files/createdModelsMinority.rds?download=1"

WORKDIR /opt

# Install R packages layer by layer to respect dependency order
RUN Rscript -e "install.packages(c('tidyverse','remotes'), repos = 'https://cran.rstudio.com/')"
RUN Rscript -e "remotes::install_github('princessmaximacenter/MnM', dependencies = T)"

# Copy the prediction script into the container
COPY scripts/MnM_predict.R /opt/scripts/MnM_predict.R
RUN chmod +x /opt/scripts/MnM_predict.R

# Set the working directory
WORKDIR /opt