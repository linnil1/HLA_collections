FROM docker.io/library/openjdk:8-jdk
RUN wget https://nagasakilab.csml.org/hla/HLAVBSeq.jar -P /opt \
    && wget https://nagasakilab.csml.org/hla/parse_result.pl -P /opt
