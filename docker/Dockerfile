FROM rocker/tidyverse:4.4.0

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Install python libraries to run format_matrix.py ; call with "/opt/venv/bin/python3 /opt/scripts/format_matrix.py ..."
RUN python3 -m venv /opt/venv
RUN /opt/venv/bin/pip install pandas numpy

# Install R libraries
RUN R -e "install.packages(c('devtools', 'caret'), repos='https://cloud.r-project.org')"
RUN R -e "devtools::install_github('jeppeseverens/AMLmapR')"

WORKDIR /opt/

# Copy the prediction script into the container
COPY scripts/AMLmapR_predict.R /opt/scripts/AMLmapR_predict.R
RUN chmod +x /opt/scripts/AMLmapR_predict.R
COPY scripts/format_matrix.py /opt/scripts/format_matrix.py
RUN chmod +x /opt/scripts/format_matrix.py
COPY scripts/example_matrix_format.csv /opt/scripts/example_matrix_format.csv