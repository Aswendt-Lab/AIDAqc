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
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

COPY aidaqc.yml /tmp/aidaqc.yml

# Create the conda environment
RUN conda env create -n aidaqc python=3.6 -f /tmp/aidaqc.yml

# Activate the environment and ensure it's activated
RUN echo "source activate aidaqc" > ~/.bashrc
ENV PATH /opt/conda/envs/aidaqc/bin:$PATH

# Set the working directory
WORKDIR /app

# Copy the rest of the application code to the container
COPY . /app

# Specify the default command to run the application
CMD ["python", "ParsingData.py"]
