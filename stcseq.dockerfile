FROM docker.io/continuumio/miniconda3:4.12.0
RUN conda install -c conda-forge -c bioconda -y bowtie=1.3.1 bowtie2=2.2.5 bbmap=38.18 r-base=4.2.0

RUN cd /usr/local/src/ \
    && wget https://ngdc.cncb.ac.cn/biocode/tools/7068/releases/v1.0/file/download?filename=STC-Seq.tar.gz -O STC-Seq.tar.gz \
    && tar -vxf STC-Seq.tar.gz \
    && rm STC-Seq.tar.gz \
    && mv STC-Seq/bin/* /usr/local/bin/ \
    && sed -i -e '/^rm/d' \
              -e '/bin\/bash/a set -e' \
              -e 's/exon-70bp/$2\/exon-70bp/g' \
              -e 's/connect_exon_70bp/$2\/connect_exon_70bp/g' \
              -e 's/hla.exon.fasta/$2\/hla.exon.fasta/g' \
              -e 's/nullAllele_list.txt/$2\/nullAllele_list.txt/g' \
              -e 's/cigar.txt/$2\/cigar.txt/g' \
              -e 's/-Xmx500m//g' \
              STC-Seq/script/step_*.sh \
    && mv STC-Seq/script/* /usr/local/bin/ \
    && chmod +x /usr/local/bin/*
