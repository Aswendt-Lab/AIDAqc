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

# Create the conda environment
RUN conda env create -n aidaqc python=3.6 && \
    conda env install -n aidaqc --file /opt/aidaqc.yaml

# Activate the environment and ensure it's activated
RUN echo "source activate aidaqc" > ~/.bashrc
ENV PATH "$(dirname (which conda))"/envs/aidaqc/bin:$PATH

# Set the working directory
WORKDIR /app

# Copy the rest of the application code to the container
COPY . /app

# Specify the default command to run the application
CMD ["python", "ParsingData.py"]
