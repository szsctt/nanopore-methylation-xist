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
# Ubuntu 20.04 ships R 3.6; current CRAN locfit requires newer R, so we pin
# an archive release when needed.
RUN Rscript -e "cran <- c('BiocManager','sys','R.utils','scales'); miss <- cran[!vapply(cran, requireNamespace, logical(1), quietly=TRUE)]; if (length(miss)) install.packages(miss, repos='https://cloud.r-project.org')" || true

RUN Rscript -e "if (!requireNamespace('locfit', quietly=TRUE)) { if (getRversion() < '4.1.0') { install.packages('https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz', repos=NULL, type='source') } else { install.packages('locfit', repos='https://cloud.r-project.org') } }"

RUN set -eux; \
    apt-get update; \
    for pkg in r-cran-biocmanager r-cran-sys r-cran-r.utils r-cran-scales r-cran-locfit r-bioc-bsgenome r-bioc-rtracklayer; do \
      apt-get install -y --no-install-recommends "$pkg" || echo "APT package not available: $pkg"; \
    done; \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e "cran <- c('BiocManager','sys','R.utils','scales'); miss <- cran[!vapply(cran, requireNamespace, logical(1), quietly=TRUE)]; if (length(miss)) install.packages(miss, repos='https://cloud.r-project.org')" \
 && Rscript -e "if (!requireNamespace('locfit', quietly=TRUE)) { if (getRversion() < '4.1.0') { install.packages('https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz', repos=NULL, type='source') } else { install.packages('locfit', repos='https://cloud.r-project.org') } }" \
 && Rscript -e "bioc <- c('rtracklayer','BSgenome','DSS'); miss <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly=TRUE)]; if (length(miss)) BiocManager::install(miss, ask=FALSE, update=FALSE)" \
 && Rscript -e "pkgs <- c('sys','R.utils','DSS','Rsamtools','GenomicRanges','IRanges','S4Vectors','Biobase'); miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if (length(miss)) stop(paste('Missing R packages:', paste(miss, collapse=', '))); cat('R lib paths:', paste(.libPaths(), collapse=':'), '\n')"

RUN printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' 'exec /opt/venv/bin/python3 -m nanomethphase "$@"' > /usr/local/bin/nanomethphase \
 && chmod +x /usr/local/bin/nanomethphase \
 && ln -sf /usr/local/bin/nanomethphase /usr/local/bin/nanomethphase.py

CMD ["/bin/bash"]
