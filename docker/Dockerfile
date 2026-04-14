FROM python:3.9-slim

RUN apt-get update && apt-get install -y wget --no-install-recommends \
    build-essential \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

#Install python dependencies
RUN pip install --upgrade pip
RUN pip install pandas numpy seaborn matplotlib scikit-learn

#Download classifier models and test data
WORKDIR /opt/models/
RUN wget -O features.tsv "https://nextcloud.computecanada.ca/index.php/s/Ax6fakA2aCSrj3e/download"
RUN wget -O classifier.pkl "https://nextcloud.computecanada.ca/index.php/s/58Eb4gMGtiBYBWN/download"

WORKDIR /opt/test_data/
RUN wget -O sample_example.tsv "https://nextcloud.computecanada.ca/index.php/s/sP3kmHJ87cDgPXN/download"

WORKDIR /opt/

#Copy the prediction script into the container
COPY scripts/Maylis_predict.py /opt/Maylis_predict.py
RUN chmod +x /opt/Maylis_predict.py