FROM debian:bookworm-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl && \
    rm -rf /var/lib/apt/lists/*

ARG OARFISH_VERSION=v0.10.0

RUN curl --proto '=https' --tlsv1.2 -LsSf \
        "https://github.com/COMBINE-lab/oarfish/releases/download/${OARFISH_VERSION}/oarfish-installer.sh" \
    | sh && \
    mv /root/.cargo/bin/oarfish /usr/local/bin/oarfish && \
    chmod +x /usr/local/bin/oarfish && \
    oarfish --version

WORKDIR /opt