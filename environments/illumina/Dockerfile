FROM continuumio/miniconda3:latest
LABEL authors="Matt Bull" \
      description="Docker image containing all requirements for an Illumina mpxv pipeline"

COPY environments/extras.yml /extras.yml
COPY environments/illumina/environment.yml /environment.yml
RUN apt-get update && apt-get install -y curl g++ git make procps && apt-get clean -y
RUN /opt/conda/bin/conda env create -f /environment.yml 
RUN /opt/conda/bin/conda env update -f /extras.yml -n artic-mpxv-illumina && /opt/conda/bin/conda clean -a
ENV PATH=/opt/conda/envs/artic-mpxv-illumina/bin:$PATH
