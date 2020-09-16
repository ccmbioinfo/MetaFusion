
#https://hub.docker.com/r/jsschrepping/velocyto-docker/dockerfile

FROM continuumio/miniconda2:4.7.12
LABEL maintainer="michael.j.apostolides@gmail.com"
ENV DEBIAN_FRONTEND noninteractive
# INSTALL BASE REQUIREMENTS
RUN apt-get update && \
    apt-get install -y git sudo apache2 python2.7 
RUN apt-get install -y python-pip 

# Update and install system libraries
RUN apt-get update --yes && \
    apt-get install --yes --no-install-recommends  \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    gcc \
    mono-mcs \
    build-essential \
    cmake \
    less \
    libbamtools-dev \
    libboost-dev \
    libboost-iostreams-dev \
    libboost-log-dev \
    libboost-system-dev \
    libboost-test-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libz-dev \
    nano \
    pandoc && \
    rm -rf /var/lib/apt/lists/*



# Update conda
RUN conda update -n base -c defaults conda

# Set channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

## Install conda packages
RUN conda install -y R=3.5.1 && \
    conda install -c bioconda -y bedtools==2.27

RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages("BiocManager", repos = "http://cran.us.r-project.org"); BiocManager::install("RBGL")' 

RUN pip install pybedtools==0.8.1
RUN pip install pandas==0.24.2

#PERL STUFF
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm Set::IntervalTree JSON::XS 

