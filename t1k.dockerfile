# Remove this when quay.io/biocontainers/t1k:v1.0.1 is ready
FROM ubuntu:20.04
RUN apt update && \
    apt install -y  g++-7 make zlib1g-dev curl git \
    && ln -s /usr/bin/g++-7 /usr/bin/g++  \
    && ln -s /usr/bin/gcc-7 /usr/bin/gcc

ARG version
RUN git clone https://github.com/mourisl/T1K /opt/T1K/ \
    && cd /opt/T1K/ \
    && git checkout ${version} \
    && make \
    && chmod +x *.pl
ENV PATH=/opt/T1K:$PATH
