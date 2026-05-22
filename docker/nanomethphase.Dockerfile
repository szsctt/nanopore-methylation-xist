# Custom multi-stage image for NanoMethPhase runtime.
# Build dependencies stay in builder stage.

FROM ubuntu:20.04 AS builder
ARG DEBIAN_FRONTEND=noninteractive
ARG NANOMETHPHASE_REF=master
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    git \
    python3 \
    python3-dev \
    python3-venv \
    python3-pip \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

RUN pip install --no-cache-dir --upgrade "pip<24" "setuptools<60" "wheel<0.38"
RUN pip install --no-cache-dir "numpy==1.18.5"
RUN pip install --no-cache-dir \
    "git+https://github.com/vahidAK/NanoMethPhase.git@${NANOMETHPHASE_REF}" \
    whatshap


FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    python3 \
    samtools \
    tabix \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

# Install DSS for nanomethphase dma subcommand.
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
 && Rscript -e "BiocManager::install('DSS', ask=FALSE, update=FALSE)"

RUN printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' 'exec /opt/venv/bin/python3 -m nanomethphase "$@"' > /usr/local/bin/nanomethphase \
 && chmod +x /usr/local/bin/nanomethphase \
 && ln -sf /usr/local/bin/nanomethphase /usr/local/bin/nanomethphase.py

CMD ["/bin/bash"]
