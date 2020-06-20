FROM ubuntu:18.04
		
ENV NEUSOMATIC_VERSION 0.3.0
ENV ZLIB_VERSION 1.2.11
ENV NUMPY_VERSION 1.18.1
ENV SCIPY_VERSION 1.4.1
ENV IMAGEIO_VERSION 2.8.0
ENV PILLOW_VERSION 7.1.2
ENV PYTORCH_VERSION 1.4.0
ENV TORCHVISION_VERSION 0.5.0
ENV CUDATOOLKIT_VERSION 9.2
ENV CMAKE_VERSION 3.14.0
ENV PYSAM_VERSION 0.15.3
ENV SAMTOOLS_VERSION 1.9
ENV TABIX_VERSION 0.2.6
ENV BEDTOOLS_VERSION 2.29.2
ENV BIOPYTHON_VERSION 1.76
ENV FISHER_VERSION 0.1.9
ENV GCC_VERSION 5

RUN apt-get update && apt-get install -y --fix-missing \
				build-essential zlib1g-dev curl less vim bzip2
RUN apt-get install -y --fix-missing git wget

RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
ENV LD_LIBRARY_PATH=/miniconda/lib:${LD_LIBRARY_PATH}
RUN conda update -y conda


RUN conda install -y zlib=${ZLIB_VERSION} numpy=${NUMPY_VERSION} scipy=${SCIPY_VERSION} \
					 pillow=${PILLOW_VERSION} cmake=${CMAKE_VERSION} imageio=${IMAGEIO_VERSION} && conda clean -a
RUN conda install -y fisher=${FISHER_VERSION} -c conda-forge && conda clean -a
RUN conda install -y pysam=${PYSAM_VERSION} \
					 samtools=${SAMTOOLS_VERSION} tabix=${TABIX_VERSION} \
					 bedtools=${BEDTOOLS_VERSION} \
					 biopython=${BIOPYTHON_VERSION} -c bioconda && conda clean -a
RUN conda install -y pytorch=${PYTORCH_VERSION} \
					 torchvision=${TORCHVISION_VERSION} \
					 cudatoolkit=${CUDATOOLKIT_VERSION} -c pytorch && conda clean -a

RUN apt-get install -y --fix-missing gcc-${GCC_VERSION} g++-${GCC_VERSION}


ADD https://github.com/bioinform/neusomatic/archive/extended_standalone.tar.gz /opt/extended_standalone.tar.gz 
RUN cd /opt/ && tar -xzvf extended_standalone.tar.gz && mv neusomatic-extended_standalone neusomatic && rm /opt/extended_standalone.tar.gz
RUN cd /opt/neusomatic/ && ./build.sh 
ENV PATH=/opt/neusomatic/neusomatic/bin:/opt/neusomatic/neusomatic/python/:${PATH}
