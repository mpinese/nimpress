FROM nimlang/nim:onbuild

MAINTAINER Mark Pinese <MPinese@ccia.org.au>

LABEL \
    description="Docker image for the nimpress polygenic score calculator"

RUN \
    apt-get update -y && \
    apt-get install -y \
        apt-utils \
        autoconf \
        build-essential \
        libbz2-dev \
        liblzma-dev

ENV HTSLIB_INSTALL_DIR=/opt/htslib

# htslib dependencies, adapted from https://github.com/brentp/hts-nim/blob/master/Dockerfile
RUN \
    mkdir -p /usr/local/include && \
    git clone --depth 1 https://github.com/ebiggers/libdeflate.git && \
    cd libdeflate && make -j4 CFLAGS="-fPIC -O3" install && \
    cd .. && rm -rf libdeflate && \
    git clone https://github.com/cloudflare/zlib cloudflare-zlib && \
    cd cloudflare-zlib && ./configure && make install && \
    cd .. && rm -rf cloudflare-zlib

# htslib, adapted from https://github.com/brentp/hts-nim/blob/master/Dockerfile
RUN \
    git clone https://github.com/samtools/htslib && \
    cd htslib && git checkout 1.10.2 && autoheader && autoconf && \
    ./configure --disable-s3 --disable-libcurl --with-libdeflate --prefix=$HTSLIB_INSTALL_DIR && \
    make -j4 CFLAGS="-fPIC -O3" install && \
    cd .. && rm -rf htslib && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

ENTRYPOINT ["./nimpress"]
