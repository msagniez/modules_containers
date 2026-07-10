FROM anaconda/miniconda:latest

WORKDIR /opt

RUN conda install -y -c bioconda oarfish && \
    conda clean -a -y

ENV PATH="/opt/conda/bin:${PATH}"

RUN oarfish -h