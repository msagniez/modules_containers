FROM rocker/r-ver:4.4
 
# --- System deps -----------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-venv \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libpng-dev \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*
 
# --- Python deps -----------------------------------------------------------
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir pandas==2.2.*
 
# --- R deps ----------------------------------------------------------------
RUN Rscript -e "install.packages('ggplot2', repos='https://cloud.r-project.org')"


WORKDIR /opt/

# Copy the prediction script into the container
COPY scripts/build_classification_summary.py  /opt/scripts/build_classification_summary.py
RUN chmod +x /opt/scripts/build_classification_summary.py
COPY scripts/plot_classification_summary.R  /opt/scripts/plot_classification_summary.R
RUN chmod +x /opt/scripts/plot_classification_summary.R
COPY scripts/run_classif_summary.sh  /opt/scripts/run_classif_summary.sh
RUN chmod +x /opt/scripts/run_classif_summary.sh

ENTRYPOINT ["/opt/scripts/run_classif_summary.sh"]
CMD ["Classif-summary"]