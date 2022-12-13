FROM alpine
RUN wget https://github.com/SyntekabioTools/HLAscan/releases/download/v2.1.4/hla_scan_r_v2.1.4 -O /usr/local/bin/hlascan \
    && chmod +x /usr/local/bin/hlascan
FROM ubuntu:22.04
COPY --from=0 /usr/local/bin/hlascan /usr/local/bin/hlascan
