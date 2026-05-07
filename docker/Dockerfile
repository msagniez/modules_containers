FROM condaforge/mambaforge:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    make \
    gcc \
    g++ \
    build-essential \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone the repository
RUN git clone https://github.com/Maggi-Chen/FusionSeeker.git /opt/fusionseeker

# Install conda dependencies into base environment
RUN mamba install -c bioconda -c conda-forge \
    minimap2 \
    pysam \
    samtools \
    -y

# Build bsalign and install to /usr/local/bin
RUN git clone https://github.com/ruanjue/bsalign.git /opt/bsalign
RUN cd /opt/bsalign && make
RUN cp /opt/bsalign/bsalign /usr/local/bin/

WORKDIR /opt/fusionseeker