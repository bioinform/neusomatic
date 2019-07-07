FROM ubuntu:16.04
		
ENV NEUSOMATIC_VERSION 0.2.1
ENV ZLIB_VERSION 1.2.11
ENV NUMPY_VERSION 1.15.4
ENV SCIPY_VERSION 1.2.0
ENV IMAGEIO_VERSION 2.5.0
ENV PYTORCH_VERSION 1.1.0
ENV TORCHVISION_VERSION 0.3.0
ENV CUDATOOLKIT_VERSION 9.0
ENV CMAKE_VERSION 3.13.2
ENV PYBEDTOOLS_VERSION 0.8.0
ENV PYSAM_VERSION 0.15.2
ENV SAMTOOLS_VERSION 1.9
ENV TABIX_VERSION 0.2.6
ENV BEDTOOLS_VERSION 2.27.1
ENV BIOPYTHON_VERSION 1.72
ENV GCC_VERSION 5

RUN apt-get update && apt-get install -y --fix-missing \
				build-essential zlib1g-dev curl less vim bzip2
RUN apt-get install -y --fix-missing git wget

RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
ENV LD_LIBRARY_PATH=/miniconda/lib:${LD_LIBRARY_PATH}
RUN conda update -y conda


RUN conda install -y zlib=${ZLIB_VERSION} numpy=${NUMPY_VERSION} scipy=${SCIPY_VERSION} \
					 imageio=${IMAGEIO_VERSION} && conda clean -a
RUN conda install -y cmake=${CMAKE_VERSION} -c conda-forge && conda clean -a
RUN conda install -y pysam=${PYSAM_VERSION} pybedtools=${PYBEDTOOLS_VERSION} \
					 samtools=${SAMTOOLS_VERSION} tabix=${TABIX_VERSION} \
					 bedtools=${BEDTOOLS_VERSION} \
					 biopython=${BIOPYTHON_VERSION} -c bioconda && conda clean -a
RUN conda install -y pytorch=${PYTORCH_VERSION} \
					 torchvision=${TORCHVISION_VERSION} \
					 cudatoolkit=${CUDATOOLKIT_VERSION} -c pytorch && conda clean -a

RUN apt-get install -y --fix-missing gcc-${GCC_VERSION} g++-${GCC_VERSION}

ADD https://github.com/bioinform/neusomatic/archive/v${NEUSOMATIC_VERSION}.tar.gz /opt/v${NEUSOMATIC_VERSION}.tar.gz 
RUN cd /opt/ && tar -xzvf v${NEUSOMATIC_VERSION}.tar.gz && mv neusomatic-${NEUSOMATIC_VERSION} neusomatic && rm /opt/v${NEUSOMATIC_VERSION}.tar.gz
RUN cd /opt/neusomatic/ && ./build.sh 
ENV PATH=/opt/neusomatic/neusomatic/bin:/opt/neusomatic/neusomatic/python/:${PATH}
