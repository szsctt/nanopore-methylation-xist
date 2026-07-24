# Custom multi-stage image for NanoMethPhase runtime.
# Build dependencies stay in the builder stage.

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
ARG DSS_VERSION=2.58.0
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    bzip2 \
    curl \
    python3 \
    samtools \
    tabix \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/venv /opt/venv
ENV MAMBA_ROOT_PREFIX="/opt/mamba"
ENV PATH="/opt/mamba/envs/dss/bin:/opt/venv/bin:${PATH}"

# Install an isolated R environment with a pinned DSS version for DMA.
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C /usr/local/bin --strip-components=1 bin/micromamba

RUN micromamba create -y -p /opt/mamba/envs/dss -c conda-forge -c bioconda \
      r-base=4.5 \
      bioconductor-dss=${DSS_VERSION} \
      r-sys \
      r-r.utils \
      bioconductor-rsamtools \
      bioconductor-genomicranges \
      bioconductor-iranges \
      bioconductor-s4vectors \
      bioconductor-biobase \
 && micromamba clean -a -y

RUN Rscript -e "v <- as.character(packageVersion('DSS')); if (v != '${DSS_VERSION}') stop(paste('Expected DSS', '${DSS_VERSION}', 'got', v)); cat('DSS version:', v, '\n')"

RUN printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' 'exec /opt/venv/bin/python3 -m nanomethphase "$@"' > /usr/local/bin/nanomethphase \
 && chmod +x /usr/local/bin/nanomethphase \
 && ln -sf /usr/local/bin/nanomethphase /usr/local/bin/nanomethphase.py

CMD ["/bin/bash"]
