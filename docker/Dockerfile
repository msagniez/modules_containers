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


WORKDIR /opt/

# Copy the prediction script into the container
COPY scripts/build_classif_summary.py  /opt/scripts/build_classif_summary.py
RUN chmod +x /opt/scripts/build_classif_summary.py