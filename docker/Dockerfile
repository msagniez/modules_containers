FROM anaconda/miniconda:latest

WORKDIR /opt

RUN conda install -c bioconda oarfish && \
    conda clean --all

ENV PATH="/opt/conda/envs/base/bin:${PATH}"

RUN oarfish -h