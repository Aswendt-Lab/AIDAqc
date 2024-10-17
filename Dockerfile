# Use the official Python 3.6 image as the base
FROM python:3.6-slim

# Set environment variables
ENV PATH /opt/conda/bin:$PATH

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda for python3.6
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh \
    && /bin/bash /opt/miniconda.sh -b -p /opt/conda \
    && rm /opt/miniconda.sh

COPY aidaqc.yaml /opt/aidaqc.yaml

# Create the conda environment
RUN /opt/conda/bin/conda env create -f /opt/aidaqc.yaml 

# Activate the environment and ensure it's activated
RUN echo "source activate aidaqc" > ~/.bashrc
ENV PATH /opt/conda/envs/aidaqc/bin:$PATH

# Set the working directory
WORKDIR /app

# Copy the rest of the application code to the container
COPY . /app

# Specify the default command to run the application
CMD ["python", "ParsingData.py"]
