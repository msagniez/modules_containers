# Use a specific R version base image from the Rocker project
FROM rocker/r-ver:4.2.0 

# Set the working directory inside the container
WORKDIR /opt

# Copy your local R script and other files into the container's working directory
COPY allcatchr /opt/allcatchr
WORKDIR /opt/allcatchr

# Install any R packages required by your script
RUN install.packages(c("devtools","tools"), repos='cran.rstudio.com')
RUN devtools::install_github("ALLCatchR")
RUN devtools::install_github("ALLCatchR2")
RUN devtools::install_github("ALLCatchRbcrabl1")

ENV PATH="/opt/allcatchr/bin:${PATH}"