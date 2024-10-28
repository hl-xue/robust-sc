FROM debian:12.6

ENV LC_ALL=C

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt update \
	&& apt-get update \
	&& apt install -y r-base \
	&& apt-get -yq install build-essential default-jdk curl \
	&& apt -y install git python3 python3-pip python3-venv

RUN cd /usr \
	&& python3 -m venv venv \
	&& . venv/bin/activate \
	&& pip3 install scipy==1.12 \ 
	&& git clone https://github.com/ayshwaryas/ddqc.git \
	&& cd ddqc/ \
	&& python3 -m pip install . \
	&& cd \
	&& rm -rf /usr/ddqc/ \
	&& pip3 install cython pegasuspy[forceatlas,rpy2,torch,scvi] scrublet scanpy liana "decoupler>=0.1.4" mofax muon omnipath marsilea mofapy2 scarches --no-cache \
	&& python3 -m pip cache purge 
 
RUN R -e "install.packages(c('BiocManager', 'RCurl', 'Seurat', 'stringr', 'dplyr'))" \
	&& R -e "BiocManager::install(c('DESeq2', 'edgeR', 'SingleCellExperiment', 'clusterProfiler') , dependencies = TRUE, ask = FALSE)"

ENV PATH="/usr/venv/bin":$PATH
