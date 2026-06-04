# SUMMARY

This repository provides scripts to parse classification results directories to build the final summary figure.

------------------------------------------------------------------------

## Prerequisites

To use this tool, you must have the following software installed on your system:

-   [Docker](https://www.docker.com/) OR [Singularity](https://sylabs.io/singularity/) OR [Apptainer](https://apptainer.org/)

------------------------------------------------------------------------

## Installation

``` bash
docker build -t summary:latest .
```

To run the summary, use the scripts located in the /opt/scripts directory.


------------------------------------------------------------------------

## Contact

## CI/CD

[![Build Status](https://github.com/msagniez/modules_containers/.github/workflows/build-push.yml/badge.svg?branch=)](https://github.com/msagniez/modules_containers/.github/workflows/build-push.yml?query=branch%3A)