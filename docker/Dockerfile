FROM rockylinux/rockylinux:9.1

# Update system and install dependencies
RUN dnf -y update && \
    dnf -y install \
        curl \
        ca-certificates && \
    dnf clean all

# Install oarfish (prebuilt binary release)
ARG OARFISH_VERSION=v0.10.0

RUN curl --proto '=https' --tlsv1.2 -LsSf \
        "https://github.com/COMBINE-lab/oarfish/releases/download/${OARFISH_VERSION}/oarfish-installer.sh" \
    | sh && \
    mv /root/.cargo/bin/oarfish /usr/local/bin/oarfish && \
    chmod +x /usr/local/bin/oarfish && \
    oarfish --version

# Set the working directory
WORKDIR /opt