FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    curl \
    git \
    pkg-config \
    unzip \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libjansson-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libtbb-dev \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libhdf5-dev \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*
    


RUN git clone https://github.com/jzhang-dev/GFAse-MaxCut-solver && \
  cd GFAse-MaxCut-solver && \
  mkdir build && \
  cd build && \
  cmake .. && \
  make -j 32 

ENV PATH="/GFAse-MaxCut-solver/build/:${PATH}"





