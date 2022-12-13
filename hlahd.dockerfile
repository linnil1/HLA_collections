FROM docker.io/continuumio/miniconda3:4.12.0
RUN conda install -y -c bioconda -c conda-forge bowtie2=2.5.0 gxx
COPY hlahd/hlahd/bin/* /usr/local/bin/
