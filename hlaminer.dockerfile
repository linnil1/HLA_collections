FROM alpine
RUN apk add perl
COPY hlaminer/HLAminer-1.4/HLAminer_v1.4/bin/* /usr/local/bin/
