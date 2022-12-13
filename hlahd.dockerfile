FROM docker.io/continuumio/miniconda3:4.12.0
RUN conda install -y -c bioconda -c conda-forge bowtie2=2.5.0 gxx
ARG folder=hlahd
COPY ${folder}/hlahd/bin/* /usr/local/bin/
COPY ${folder}/hlahd/src/* /usr/local/src/
COPY ${folder}/hlahd/compile.sh /usr/local/
RUN cd /usr/local/ && bash compile.sh
