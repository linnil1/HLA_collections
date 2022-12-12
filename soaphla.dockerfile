FROM docker.io/continuumio/miniconda3:4.12.0
RUN conda install -y -c adefelicibus -c bioconda -c conda-forge soap-hla=1.0.0 samtools
# bug https://github.com/adefelicibus/soap-hla/issues/2
RUN wget https://raw.githubusercontent.com/adefelicibus/soap-hla/407812d90acc37ae54bd1543d6984171636d01ce/src/MHC_autopipeline.pl -O /opt/conda/bin/MHC_autopipeline \
    && chmod +x /opt/conda/bin/MHC_autopipeline
ENV CONDA_PREFIX=/opt/conda
