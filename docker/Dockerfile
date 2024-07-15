FROM ubuntu:24.04

ENV LC_ALL=C

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get update \
	&& apt-get install -yq build-essential wget default-jdk \
	&& wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ./miniconda.sh \
	&& /bin/bash ./miniconda.sh -b -p /opt/conda \
	&& rm -f ./miniconda.sh \
	&& /opt/conda/bin/conda init bash
ENV PATH=/opt/conda/bin:$PATH
RUN conda install -y conda-forge::r-base conda-forge::python conda-forge::r-stringr conda-forge::r-devtools bioconda::bioconductor-singlecellexperiment bioconda::bioconductor-dropletutils conda-forge::r-seurat conda-forge::r-hdf5r conda-forge::r-seuratobject conda-forge::r-matrix conda-forge::scanpy bioconda::pegasuspy bioconda::pegasusio bioconda::bioconductor-zellkonverter \
	&& R -e 'devtools::install_github("JSB-UCLA/scDEED")'
