FROM docker.io/continuumio/miniconda3:4.12.0
RUN conda install -c defaults -c conda-forge -c bioconda -y gxx=12.2.0 make boost=1.80 r-base=4.2.2 blat=35 unzip

RUN cd /usr/local/src \
    && wget https://www.ikmb.uni-kiel.de/sites/default/files/downloads/hlassign/pilot_source_code.zip \
    && unzip pilot_source_code.zip \
    && cd pilot_source_code \
    && unzip NextCallHLA.zip \
    && unzip getHLAMultiAlignExonWise.zip \
    && unzip prefilterHlaReads.zip \
    && sed -i -e 's/-lboost_regex//g' -e 's/g++/g++ -I\/opt\/conda\/include/g' */nbproject/Makefile-Debug.mk \
    && sed -i '/void erase(iterator i){return the_alignments.mapMultipleAlignment.erase(i);}/ s/return//'  NextCallHLA/CImgtDb.h \
    && cd getHLAMultiAlignExonWise && make -j4 && cd .. \
    && cd NextCallHLA && make -j4 && cd .. \
    && cd prefilterOnTarget && make -j4 && cd .. \
    && cp */dist/Debug/GNU-Linux-x86/* /usr/local/bin
