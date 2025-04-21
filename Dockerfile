# Use the official Python 3.6 image as the base
FROM python:3.6-slim

# Set environment variables
ENV PATH /opt/conda/bin:$PATH

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

FROM continuumio/miniconda:latest

COPY aidaqc.yaml /opt/aidaqc.yaml

SHELL ["/bin/bash", "-l", "-c"]
# Create the conda environment
RUN conda create -n aidaqc python=3.6 && \
    conda env update -n aidaqc --file /opt/aidaqc.yaml

RUN useradd -ms /bin/bash aida
USER aida
SHELL ["conda", "run", "-n", "aidaqc", "/bin/bash", "-c"]

# Activate the environment and ensure it's activated
RUN echo "source activate aidaqc" > ~/.bashrc
ENV PATH "$(dirname (dirname (which conda)))"/envs/aidaqc/bin:$PATH
RUN /bin/bash -c "source activate aidaqc"

# Set the working directory
WORKDIR /app

# Copy the rest of the application code to the container
COPY . /app

ENTRYPOINT ["conda", "run", "-n", "aidaqc", "python", "/app/scripts/ParsingData.py"]
