FROM condaforge/mambaforge:latest

# Clone the repository
RUN git clone https://github.com/jwanglab/fusilli.git /opt/fusilli

# Create conda environment and install dependencies
RUN conda create -n fusilli python -y
RUN conda run -n fusilli conda install pip -y
RUN conda run -n fusilli pip install -r /opt/fusilli/requirements.txt

WORKDIR /opt/fusilli

# Set the env's bin as default so you can call fusilli directly
ENV PATH="/opt/conda/envs/fusilli/bin:$PATH"