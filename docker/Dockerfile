FROM docker.io/mambaorg/micromamba:latest

WORKDIR /opt

RUN micromamba install -y -n base -c bioconda -c conda-forge oarfish && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"

RUN oarfish -h