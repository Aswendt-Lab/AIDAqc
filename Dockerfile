# Multi-architecture AIDAqc Docker image
# linux/amd64 -> aidaqc-intel.yaml
# linux/arm64 -> aidaqc-arm64.yaml

FROM condaforge/miniforge3:latest

ARG TARGETARCH

COPY aidaqc-arm64.yaml /opt/aidaqc-arm64.yaml
COPY aidaqc-intel.yaml /opt/aidaqc-intel.yaml

SHELL ["/bin/bash", "-lc"]

RUN if [ "${TARGETARCH}" = "arm64" ]; then \
        echo "Building AIDAqc ARM64 environment"; \
        mamba env create -n aidaqc -f /opt/aidaqc-arm64.yaml; \
    elif [ "${TARGETARCH}" = "amd64" ]; then \
        echo "Building AIDAqc AMD64 environment"; \
        mamba env create -n aidaqc -f /opt/aidaqc-intel.yaml; \
    else \
        echo "Unsupported architecture: ${TARGETARCH}"; \
        exit 1; \
    fi && \
    conda clean -afy

ENV PATH=/opt/conda/envs/aidaqc/bin:$PATH

RUN useradd -m -s /bin/bash aida

WORKDIR /app

COPY --chown=aida:aida . /app

USER aida

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "aidaqc", "python", "/app/scripts/ParsingData.py"]