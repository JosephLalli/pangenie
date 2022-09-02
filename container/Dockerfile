# PanGenie docker file.
# In this Dockerfile PanGenie is installed via git cloning.

FROM ubuntu:20.04 as builder

RUN apt-get update && \
    apt-get install --assume-yes software-properties-common && \
    add-apt-repository universe && \
    apt-get update && \
    apt-get install --no-install-recommends -qqy \
    build-essential \
    git \
    zlib1g-dev \
    libjellyfish-2.0-dev \
    cmake \
    pkg-config

RUN mkdir /metadata && echo `dpkg -l | grep jellyfish | tr -s " " | cut -d " " -f 2,3` > /metadata/jellyfish.lib.version

RUN mkdir /repos && cd /repos && \
    git clone https://github.com/eblerjana/pangenie.git && \
    cd pangenie && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 4 && \
    cd .. && \
    echo `git rev-parse --short HEAD` > /metadata/pangenie.git.version

FROM jlalli/samtools-with-libdeflate:1.16.0.1.13.0 as runner
ENV LD_LIBRARY_PATH=/usr/local/lib
COPY --from=builder /repos/pangenie/build/src/PanGenie /usr/local/bin/pangenie
COPY --from=builder /repos/pangenie/build/src/PanGenie-graph /usr/local/bin/PanGenie-graph
COPY --from=builder /repos/pangenie/build/src/libPanGenieLib.so /usr/local/lib/libPanGenieLib.so
COPY --from=builder /metadata/pangenie.git.version /pangenie.git.version
COPY --from=builder /usr/lib/x86_64-linux-gnu/libjellyfish-2.0* /usr/lib/x86_64-linux-gnu/
RUN apt-get update -y && apt-get install -y libcurl4-openssl-dev jellyfish pigz

ENV LC_ALL=C