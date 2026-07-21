FROM ubuntu:18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG CRAN_SNAPSHOT=2020-04-01
ARG READR_CRAN_SNAPSHOT=2021-09-28
ARG BIOC_VERSION=3.10
ARG NUCDYN_REF=master
ARG R_INSTALL_CPUS=2

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    CRAN_REPO=https://packagemanager.posit.co/cran/${CRAN_SNAPSHOT} \
    READR_CRAN_REPO=https://packagemanager.posit.co/cran/${READR_CRAN_SNAPSHOT} \
    BIOC_VERSION=${BIOC_VERSION} \
    R_INSTALL_CPUS=${R_INSTALL_CPUS}

# System utilities, R package compilation tools, and common native libraries.
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
        apt-transport-https \
        build-essential \
        ca-certificates \
        curl \
        gfortran \
        git \
        libbz2-dev \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libjpeg-dev \
        liblzma-dev \
        libpcre3-dev \
        libpng-dev \
        libssl-dev \
        libtiff5-dev \
        libxml2-dev \
        software-properties-common \
        tcsh \
        vim \
        wget \
        zlib1g-dev; \
    rm -rf /var/lib/apt/lists/*

# Install the CRAN build of R 3.6.3 for Ubuntu 18.04 (bionic).
RUN set -eux; \
    wget -qO /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
        https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc; \
    echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
        > /etc/apt/sources.list.d/cran-r.list; \
    apt-get update; \
    apt-get install -y --no-install-recommends r-base r-base-dev; \
    R --version; \
    Rscript --vanilla -e 'stopifnot(getRversion() >= "3.6", getRversion() < "3.7")'; \
    rm -rf /var/lib/apt/lists/* /tmp/*

# Make the historical CRAN snapshot the default for every R invocation.
RUN printf '%s\n' \
    'local({' \
    '  repo <- Sys.getenv("CRAN_REPO", unset = "https://packagemanager.posit.co/cran/2020-04-01")' \
    '  options(repos = c(CRAN = repo))' \
    '})' \
    >> /usr/lib/R/etc/Rprofile.site

WORKDIR /home/NucleosomeDynamics

# Install all CRAN and Bioconductor dependencies directly in the image.
RUN printf '%s\n' \
    'cran_repo <- Sys.getenv("CRAN_REPO")' \
    'bioc_version <- Sys.getenv("BIOC_VERSION", unset = "3.10")' \
    'ncpus <- as.integer(Sys.getenv("R_INSTALL_CPUS", unset = "2"))' \
    'options(repos = c(CRAN = cran_repo), Ncpus = ncpus)' \
    'install.packages("BiocManager")' \
    'cran_packages <- c("getopt", "plyr", "purrr", "dplyr", "ggplot2", "doParallel", "magrittr", "Rcpp", "knitr", "rmarkdown")' \
    'install.packages(cran_packages)' \
    'repos <- suppressWarnings(BiocManager::repositories(version = bioc_version))' \
    'repos["CRAN"] <- cran_repo' \
    'options(repos = repos)' \
    'bioc_packages <- c("Biobase", "BiocGenerics", "Biostrings", "GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools", "S4Vectors", "ShortRead", "BSgenome", "plyranges", "nucleR")' \
    'BiocManager::install(bioc_packages, version = bioc_version, ask = FALSE, update = FALSE, Ncpus = ncpus)' \
    'required <- unique(c(cran_packages, bioc_packages))' \
    'missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]' \
    'if (length(missing)) stop("Failed to install required R packages: ", paste(missing, collapse = ", "))' \
    'cat("R:", as.character(getRversion()), "\n")' \
    'cat("Bioconductor:", as.character(BiocManager::version()), "\n")' \
    'for (pkg in c("dplyr", "ggplot2", "plyranges", "nucleR")) cat(pkg, as.character(packageVersion(pkg)), "\n")' \
    | Rscript --vanilla -


# Pin rlang and readr to versions required by the downstream analysis code.
# This is done after the general dependency installation so these versions are
# the ones visible when NucDyn and the final runtime environment are validated.
USER root

# Remove any abandoned package-installation locks.
RUN rm -rf /usr/local/lib/R/site-library/00LOCK-*

# Install a version of rlang that exports inject(). Use a separate R process so
# an older rlang namespace cannot remain loaded during verification.
RUN R --vanilla -e '\
    install.packages( \
        "https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.11.tar.gz", \
        repos = NULL, \
        type = "source" \
    ); \
    cat("rlang version: "); \
    print(packageVersion("rlang")); \
    stopifnot(packageVersion("rlang") == package_version("0.4.11")); \
    stopifnot("inject" %in% getNamespaceExports("rlang")) \
'

# Install readr and its historical dependencies from the 2021-09-28 snapshot.
RUN R --vanilla -e '\
    options(repos = c(CRAN = Sys.getenv("READR_CRAN_REPO"))); \
    install.packages( \
        "readr", \
        dependencies = c("Depends", "Imports", "LinkingTo") \
    ); \
    cat("readr version: "); \
    print(packageVersion("readr")); \
    stopifnot(packageVersion("readr") == package_version("2.0.2")); \
    stopifnot("show_col_types" %in% names(formals(readr::read_delim))) \
'

# Verify the pinned packages in a fresh R process.
RUN R --vanilla -e '\
    cat("R: ", R.version.string, "\n"); \
    cat("rlang: ", as.character(packageVersion("rlang")), "\n"); \
    cat("readr: ", as.character(packageVersion("readr")), "\n"); \
    cat( \
        "rlang::inject exported: ", \
        "inject" %in% getNamespaceExports("rlang"), \
        "\n" \
    ); \
    cat( \
        "show_col_types supported: ", \
        "show_col_types" %in% names(formals(readr::read_delim)), \
        "\n" \
    ); \
    stopifnot(packageVersion("rlang") == package_version("0.4.11")); \
    stopifnot(packageVersion("readr") == package_version("2.0.2")); \
    stopifnot("inject" %in% getNamespaceExports("rlang")); \
    stopifnot("show_col_types" %in% names(formals(readr::read_delim))) \
'
# Install NucDyn. For a fully reproducible build, override NUCDYN_REF with a
# Git commit SHA rather than leaving it as "master".
RUN set -eux; \
    git clone https://github.com/nucleosome-dynamics/NucDyn.git NucDyn; \
    cd NucDyn; \
    git checkout "${NUCDYN_REF}"; \
    cd ..; \
    R CMD INSTALL NucDyn; \
    Rscript --vanilla -e 'stopifnot(requireNamespace("NucDyn", quietly = TRUE)); cat("Installed NucDyn:", as.character(packageVersion("NucDyn")), "\n")'; \
    rm -rf NucDyn/.git /tmp/*

# Copy the Nucleosome Dynamics CLI from the Docker build context.
COPY nucleosome_dynamics/ /home/NucleosomeDynamics/NucleosomeDynamics/

# Install UCSC command-line utilities required by the CLI.
RUN set -eux; \
    mkdir -p /home/NucleosomeDynamics/NucleosomeDynamics/wig_utils; \
    wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig \
        -O /home/NucleosomeDynamics/NucleosomeDynamics/wig_utils/wigToBigWig; \
    wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig \
        -O /home/NucleosomeDynamics/NucleosomeDynamics/wig_utils/bigWigToWig; \
    wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes \
        -O /home/NucleosomeDynamics/NucleosomeDynamics/wig_utils/fetchChromSizes; \
    chmod 0755 /home/NucleosomeDynamics/NucleosomeDynamics/wig_utils/*; \
    test -s /home/NucleosomeDynamics/NucleosomeDynamics/wig_utils/wigToBigWig; \
    test -s /home/NucleosomeDynamics/NucleosomeDynamics/wig_utils/bigWigToWig; \
    test -s /home/NucleosomeDynamics/NucleosomeDynamics/wig_utils/fetchChromSizes

ENV PATH=/home/NucleosomeDynamics/NucleosomeDynamics/wig_utils:${PATH}

COPY runNuclDyn /home/NucleosomeDynamics/runNuclDyn
COPY help/ /home/NucleosomeDynamics/help/

# Validate the final R environment and prepare the non-root runtime user.
RUN set -eux; \
    chmod 0755 /home/NucleosomeDynamics/runNuclDyn; \
    printf '%s\n' \
        'required <- c("NucDyn", "nucleR", "getopt", "BSgenome", "plyranges", "plyr", "purrr", "dplyr", "readr", "ggplot2")' \
        'missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]' \
        'if (length(missing)) stop("Final image is missing R packages: ", paste(missing, collapse = ", "))' \
        'stopifnot(packageVersion("rlang") == package_version("0.4.11"))' \
        'stopifnot(packageVersion("readr") == package_version("2.0.2"))' \
        'stopifnot("inject" %in% getNamespaceExports("rlang"))' \
        'stopifnot("show_col_types" %in% names(formals(readr::read_delim)))' \
        | Rscript --vanilla -; \
    groupadd --gid 99990 nucdyn; \
    useradd --uid 99990 --gid 99990 --home-dir /home/NucleosomeDynamics \
        --no-create-home --shell /bin/bash nucdyn; \
    chown -R nucdyn:nucdyn /home/NucleosomeDynamics

USER nucdyn:nucdyn
WORKDIR /home/NucleosomeDynamics
ENTRYPOINT ["/home/NucleosomeDynamics/runNuclDyn"]

LABEL base.image="ubuntu:18.04" \
      software="NucleosomeDynamics" \
      software.version="R 3.6.3; Bioconductor 3.10; rlang 0.4.11; readr 2.0.2; NucDyn; nucleR" \
      version="v1.0" \
      about.summary="R-based package for nucleosome positioning analysis using MNase-seq data" \
      about.home="https://mmb.irbbarcelona.org/NucleosomeDynamics" \
      about.license="Apache-2.0" \
      maintainer="Josep Ll. Gelpi <gelpi@ub.edu>"
