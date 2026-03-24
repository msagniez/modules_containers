# ALLcatchR Classifiers

This repository provides:

1\. The ALLcatchR Docker container definition.

------------------------------------------------------------------------

## Prerequisites

To use this classifier, you must have the following software installed on your system:

-   [Docker](https://www.docker.com/) OR [Singularity](https://sylabs.io/singularity/) OR [Apptainer](https://apptainer.org/)

------------------------------------------------------------------------

## Installation

``` bash
docker build -t allcatchr:latest .
```

To run the full set of ALLcatchR classifiers (lineage, B-ALL, T-ALL, BCR-ABL1), use the scripts located in the /scripts directory.


------------------------------------------------------------------------

## License

This project is licensed under the following licenses:
[ALLcatchR MIT LICENSE](https://github.com/ThomasBeder/ALLCatchR/blob/main/LICENSE.txt)
[ALLcatchR2 MIT LICENSE](https://github.com/ThomasBeder/ALLCatchR2/blob/main/LICENSE)

------------------------------------------------------------------------

## Contact

## CI/CD

[![Build Status](https://github.com/msagniez/modules_containers/.github/workflows/build-push.yml/badge.svg?branch=)](https://github.com/msagniez/modules_containers/.github/workflows/build-push.yml?query=branch%3A)