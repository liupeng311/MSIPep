FROM condaforge/mambaforge:24.7.1-0

LABEL description="MSIPep immunopeptidomics-transcriptomics pipeline"

ENV DEBIAN_FRONTEND=noninteractive \
    MSIPEP_ROOT=/opt/msipep \
    MSIPEP_DATA=/data \
    MSIPEP_REF=/ref \
    MSIPEP_SOFTWARE=/software

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget unzip curl ca-certificates procps bash \
    && rm -rf /var/lib/apt/lists/*

WORKDIR ${MSIPEP_ROOT}
COPY . ${MSIPEP_ROOT}/

RUN chmod +x docker/install-conda-envs.sh docker/entrypoint.sh \
    && bash docker/install-conda-envs.sh

RUN mkdir -p /data /ref /software
ENV PATH="${MSIPEP_ROOT}:${PATH}"

ENTRYPOINT ["/opt/msipep/docker/entrypoint.sh"]
CMD ["help"]
