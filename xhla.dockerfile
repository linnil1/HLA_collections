FROM docker.io/humanlongevity/hla
RUN apt update -y && apt install -y muscle parallel
ENTRYPOINT []
CMD bash
