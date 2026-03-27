# MnM Classifier

This repository provides:

1\. The MnM Docker container definition.
2\. A script to run MnM predictions from gene-level counts.

------------------------------------------------------------------------

## Prerequisites

To use this tool, you must have the following software installed on your system:

-   [Docker](https://www.docker.com/) OR [Singularity](https://sylabs.io/singularity/) OR [Apptainer](https://apptainer.org/)

------------------------------------------------------------------------

## Installation

``` bash
docker build -t mnm:latest .
```

To run MnM, use the scripts located in the /opt/scripts directory.


------------------------------------------------------------------------

## License

This project is licensed under the following licenses:
[GNU GENERAL PUBLIC LICENSE](https://github.com/princessmaximacenter/MnM/blob/main/LICENSE.md)

------------------------------------------------------------------------

## Contact

## CI/CD

[![Build Status](https://github.com/msagniez/modules_containers/.github/workflows/build-push.yml/badge.svg?branch=)](https://github.com/msagniez/modules_containers/.github/workflows/build-push.yml?query=branch%3A)