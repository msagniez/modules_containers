FROM python:3.10-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/

# ── Clone AttentionAML ────────────────────────────────────────────────────────
RUN git clone https://github.com/wan-mlab/AttentionAML.git

# ── Install Python dependencies ───────────────────────────────────────────────
# Core scientific stack expected by Wan Lab tools
RUN pip install --no-cache-dir \
    numpy \
    pandas \
    scikit-learn \
    scipy \
    torch \
    tqdm

# ── Install the package itself ────────────────────────────────────────────────
RUN pip install --no-cache-dir ./AttentionAML


WORKDIR /opt/AttentionAML

# Copy the prediction script into the container
COPY scripts/AttentionAML_predict.py /opt/AttentionAML/AttentionAML_predict.py
RUN chmod +x /opt/AttentionAML/AttentionAML_predict.py
COPY scripts/format_matrix.py /opt/AttentionAML/format_matrix.py
RUN chmod +x /opt/AttentionAML/format_matrix.py
COPY scripts/TPM_test_format.csv /opt/AttentionAML/TPM_test_format.csv