
FROM mambaorg/micromamba:latest

RUN micromamba install -y -n base -c bioconda -c conda-forge oarfish && \
    micromamba clean --all --yes

WORKDIR /opt