# tximport tools

This repository provides:

1\. The tximport Docker container definition.
2\. A script to summarize gene-level counts and TPMs from Oarfish transcript-level quantifications.

------------------------------------------------------------------------

## Prerequisites

To use this tool, you must have the following software installed on your system:

-   [Docker](https://www.docker.com/) OR [Singularity](https://sylabs.io/singularity/) OR [Apptainer](https://apptainer.org/)

------------------------------------------------------------------------

## Installation

``` bash
docker build -t tximport:latest .
```

To run gene summarization, use the scripts located in the /opt/scripts directory.


------------------------------------------------------------------------

## License

This project is licensed under the following licenses:
[GNU GENERAL PUBLIC LICENSE](https://bioconductor.org/packages/release/bioc/html/tximport.html)

------------------------------------------------------------------------

## Contact

## CI/CD

[![Build Status](https://github.com/msagniez/modules_containers/.github/workflows/build-push.yml/badge.svg?branch=)](https://github.com/msagniez/modules_containers/.github/workflows/build-push.yml?query=branch%3A)