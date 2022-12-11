FROM alpine
RUN wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.7.5/graphtyper -O /usr/local/bin/graphtyper && \
    chmod a+x /usr/local/bin/graphtyper
