FROM condaforge/mambaforge:latest

# Clone the repository
RUN git clone https://github.com/jwanglab/fusilli.git /opt/fusilli

# Install dependencies into base environment
RUN pip install --ignore-installed -r /opt/fusilli/requirements.txt

WORKDIR /opt/fusilli