FROM docker.io/library/ubuntu:22.04
RUN apt update -y && apt install -y python3-hisat2=2.2.1-3
COPY hisat/hisat-genotype /opt/hisat-genotype
RUN ln -s /usr/bin/python3 /usr/bin/python
ENV PATH=/opt/hisat-genotype:$PATH
ENV PYTHONPATH=/opt/hisat-genotype/hisatgenotype_modules
