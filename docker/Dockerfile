FROM rocker/r-ver:4.2.0

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
RUN Rscript -e "install.packages('curl')"
RUN Rscript -e "install.packages(c('fs', 'httr2', 'credentials'), repos='https://cran.rstudio.com/')"
RUN Rscript -e "install.packages(c('gert', 'gh', 'ragg'), repos='https://cran.rstudio.com/')"
RUN Rscript -e "install.packages(c('bslib', 'rmarkdown', 'shiny'), repos='https://cran.rstudio.com/')"
RUN Rscript -e "install.packages(c('pkgload', 'rcmdcheck', 'usethis', 'miniUI', 'pkgdown'), repos='https://cran.rstudio.com/')"
RUN Rscript -e "install.packages(c('roxygen2', 'testthat', 'profvis', 'rversions', 'urlchecker'), repos='https://cran.rstudio.com/')"
RUN Rscript -e "install.packages(c('devtools', 'remotes'), repos='https://cran.rstudio.com/')"
RUN Rscript -e "install.packages(c('BiocManager', 'dplyr'), repos='https://cran.rstudio.com/')"
RUN Rscript -e "BiocManager::install(c('singscore'))"


# Install ALLCatchR packages from GitHub
RUN Rscript -e "devtools::install_github('ThomasBeder/ALLCatchR')"
RUN Rscript -e "devtools::install_github('ThomasBeder/ALLCatchR2')"
RUN Rscript -e "devtools::install_github('ThomasBeder/ALLCatchR_bcrabl1')"

# Copy the prediction script into the container
COPY scripts/ALLcatchR_predict.R /opt/scripts/ALLcatchR_predict.R
RUN chmod +x /opt/scripts/ALLcatchR_predict.R

# Set the working directory
WORKDIR /opt