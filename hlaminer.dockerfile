FROM alpine
RUN wget https://github.com/bcgsc/HLAminer/releases/download/v1.4/HLAminer_1-4.tar.gz \
    && tar -vxf HLAminer_1-4.tar.gz \
    && mv HLAminer-1.4/HLAminer_v1.4/bin/* /usr/local/bin/ \
    && rm HLA* -rf
FROM ubuntu:22.04
RUN apt update && apt install -y perl
COPY --from=0 /usr/local/bin/* /usr/local/bin/
