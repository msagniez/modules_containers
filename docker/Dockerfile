FROM rocker/r-ver:4.5.0

RUN apt-get update && apt-get install -y wget --no-install-recommends \
    build-essential \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

#Install R libraries
RUN Rscript -e "install.packages(c('caret','kernlab', 'AMLmapR', 'tools'), repos = 'https://cran.rstudio.com/')"

WORKDIR /opt/

#Copy the prediction script into the container
COPY scripts/AMLmapR_predict.R /opt/scripts/AMLmapR_predict.R
RUN chmod +x /opt/scripts/AMLmapR_predict.R
COPY scripts/format_matrix.py /opt/scripts/format_matrix.py
RUN chmod +x /opt/scripts/format_matrix.py
COPY scripts/example_matrix_format.csv /opt/scripts/example_matrix_format.csv