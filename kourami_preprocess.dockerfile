FROM docker.io/continuumio/miniconda3:4.12.0
RUN conda install -y -c conda-forge -c bioconda bwakit=0.7.17 bamUtil=1.0.14 samtools=1.15.1
