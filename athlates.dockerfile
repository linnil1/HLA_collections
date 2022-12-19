FROM docker.io/continuumio/miniconda3:4.12.0
RUN conda install -c bioconda -c conda-forge -y perl samtools=1.14 novoalign=3.09.04 bamtools=2.5.2 gxx gcc=12.2.0

RUN git clone https://github.com/cliu32/athlates /usr/local/athlates \
    && cd /usr/local/athlates \
    && git checkout 40c7334 \
    && chmod +x *.pl \
    && cp *.pl /usr/local/bin

ARG folder=athlates
COPY ${folder}/athlates/src /usr/local/src
# see makefile
RUN g++ -O3 -I/opt/conda/include/bamtools/ /usr/local/src/*.cpp -o /usr/local/bin/typing -lpthread -lbamtools -lz
RUN mkdir -p /path/to/athlates/bin/ \
    && ln -s /usr/local/bin/typing /path/to/athlates/bin/typing

