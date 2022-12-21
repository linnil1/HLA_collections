FROM quay.io/biocontainers/seq2hla:2.2--2
RUN mv /usr/local/share/seq2hla-2.2-2/*dbmhc /usr/local/share/seq2hla-2.2-2/references \
    && ln -s references/HLAI.dbmhc /usr/local/share/seq2hla-2.2-2/HLAI.dbmhc \
    && ln -s references/HLAII.dbmhc /usr/local/share/seq2hla-2.2-2/HLAII.dbmhc \

