
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
#RUN conda install -y numpy=1.16.1 \
#    	  	 scipy=1.2.0 \
#		     cython=0.29.4 \
#		     numba=0.41.0 \
#		     h5py=2.9.0 \
#		     click=7.0 \
#		     R=3.5.1 \
#		     rpy2=2.9.4 \
#		     r-devtools=2.0.1  && \
#    conda install -c bioconda -y samtools=1.9 && \
#    conda install -c bioconda -y bedtools==2.27 
RUN conda install -y R=3.5.1 && \
    conda install -c bioconda -y bedtools==2.27

RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages("BiocManager", repos = "http://cran.us.r-project.org"); BiocManager::install("RBGL")' 

RUN pip install pybedtools==0.8.1
RUN pip install pandas==0.24.2

#PERL STUFF
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm Set::IntervalTree JSON::XS 
#perl -MCPAN -e shell
#install Set::IntervalTree
#install JSON::XS 

#FusionAnnotator - BUILD CTAT LIB
#RUN git clone https://github.com/FusionAnnotator/CTAT_HumanFusionLib.git
#RUN perl /CTAT_HumanFusionLib/build_fusion_lib.pl > /CTAT_HumanFusionLib/fusion_lib.dat



#COPY /Users/mapostolides/MetaFusion/scripts /home 
#SCRAP
#RUN \
#  R -e 'chooseCRANmirror(ind=52); install.packages(c("devtools", "Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown", "withr"))'
#RUN pip install pysam==0.15.2 \
#    		loompy==2.0.17 \
#		velocyto==0.17.17
#RUN pip install pysam==0.15.2 
#RUN pip install bedtools==2.27
