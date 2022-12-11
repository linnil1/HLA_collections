FROM docker.io/biocontainers/biocontainers:v1.2.0_cv2
RUN conda install -y -c conda-forge bwakit=0.7.17 bamUtil=1.0.14 samtools=1.15.1
