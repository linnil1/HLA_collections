FROM docker.io/continuumio/miniconda3:4.12.0
RUN conda install -c bioconda -c conda-forge bowtie=1.3.1 perl-bioperl=1.7.8 perl-list-moreutils=0.430

RUN git clone https://github.com/FNaveed786/hlaforest.git \
    && cd hlaforest \
    && git checkout 75edb46 \
    && chmod +x scripts/*.sh \
    && cp -r scripts/* /usr/local/bin \
    && mkdir /root/hla/hlaforest/scripts/ -p \
    && ln -s /usr/local/bin/config.sh.sh /root/hla/hlaforest/scripts/ \
    && rm ../hlaforest -rf

# linnil1 bug fix
RUN sed -i \
        -e 's/$BOWTIE_INDEX/-x\ $BOWTIE_INDEX/g' \
        -e 's/FILTERED_PREFIX\ /FILTERED_PREFIX.sam /g' \
        -e 's/_1/_1.sam/g' \
        -e 's/_2/_2.sam/g' \
        /usr/local/bin/build-forest-from-fastq.sh
