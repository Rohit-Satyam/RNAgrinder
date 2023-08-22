FROM ubuntu:20.04

######################
# Prerequisites
#######################
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -yq \
build-essential \
zlib1g-dev \
apt-utils \
# gnupg requirement
dirmngr \
gnupg \
# curl requirement
ca-certificates \
curl \
cmake \
vim \
procps \
tabix \
pkg-config \
whiptail \
# install script requirements
sudo \
wget \
locales \
liblzma-dev \
libbz2-dev \
git \
libtinfo6 \
libcurl4-gnutls-dev 

# Install miniconda
RUN curl -LO "https://repo.anaconda.com/miniconda/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh"
RUN chmod +x Miniconda3-py39_22.11.1-1-Linux-x86_64.sh
RUN bash Miniconda3-py39_22.11.1-1-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-py39_22.11.1-1-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda init

## Installing svelter
RUN mkdir /bin/wham && git clone --recursive  https://github.com/zeeev/wham.git /bin/wham &&\
cd /bin/wham/src/bamtools/ && \
mkdir lib && \
git checkout master && \
cmake -DCMAKE_INSTALL_PREFIX=. && \
make && \
cp src/libbamtools.a lib/ && \
cd /bin/wham && \
make
RUN conda install -y -c conda-forge -n base mamba


#RUN git clone https://github.com/leklab/haplocheckCLI.git && cd haplocheckCLI &&  \
#jar cvfe haplocheckCLI.jar haplocheck_contam * && chmod 777 *
#ENV PATH "$PATH:haplocheckCLI"
#RUN echo 'export PATH="haplocheckCLI:${PATH}"' >> /root/.bashrc
#RUN echo 'alias haplocheckCLI="java -jar haplocheckCLI/haplocheckCLI.jar"' >> /root/.bashrc
#RUN source $HOME/.bashrc
## Adding the pipeline

ADD haplocheckCLI haplocheckCLI/
#ADD bin bin/
COPY bin/Post-variant-call-pipeline.Rmd bin/
COPY bin/collectwgsmetrics.R bin/
COPY bin/coverageAtEveryBase.R bin/
COPY bin/getContaminationAndFilter.sh bin/
COPY bin/header bin/
COPY bin/maftools.R bin/
COPY bin/multiqc_config.yaml bin/
ADD modules modules/
ADD test test/
ADD nextflow.config $HOME
ADD main.nf $HOME
ADD download_gatk_resources.sh $HOME
ADD bin/vcf2maf bin/vcf2maf
RUN wget https://github.com/dellytools/delly/releases/download/v1.1.6/delly_v1.1.6_linux_x86_64bit -P bin/

RUN mamba install -y -c conda-forge -c bioconda -n base bioconductor-maftools \
bioconductor-biocstyle \ 
r-dygraphs \
r-r.utils \
bioconductor-genomicranges \
bioconductor-rtracklayer \
bioconductor-bsgenome \
gnuplot  

RUN mamba create -y -c conda-forge -c bioconda -n whamg wham breakdancer=1.4.5
RUN R -e "install.packages(c('BiocManager','NMF'),dependencies=TRUE, repos='http://cran.rstudio.com/');  if (!library(BiocManager, logical.return=T)) quit(status=10)"
RUN R -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38"); if (!library(BSgenome.Hsapiens.UCSC.hg38, logical.return=T)) quit(status=10)'
RUN mamba install -c biobuilds mrfast

## Some separate environment for tools with conflicting dependencies
RUN mamba create -y -n breakseq python=2.7 numpy && \
mamba run -n breakseq pip install https://github.com/bioinform/breakseq2/archive/2.2.tar.gz

RUN mamba install -y -c bioconda -c conda-forge -n breakseq lumpy-sv \
manta 

## Installing svelter
RUN git clone https://github.com/mills-lab/svelter.git /bin/svelter && \
cd /bin/svelter && \
mamba run -n breakseq python setup.py install && \
sed -i-e '351,353d;849,851d;1825,1828d;1946,1949d;6381,6384d;10935,10938d' svelter_sv/svelter.py && \
sed -i "s/if not chrom_single in chromos:/chromos=chrom_single.split(',')/g" svelter_sv/svelter.py

RUN git clone https://github.com/BilkentCompGen/tardis.git --recursive /bin/tardis && \
cd /bin/tardis && \
make libs && \
make

RUN mamba create -y -n cnvnator -c conda-forge -c bioconda cnvnator

RUN mamba install -y -n base -c conda-forge unzip
RUN wget https://github.com/SFGLab/ConsensuSV-core/archive/refs/tags/1.7.zip -P bin/ && \
cd bin/ && \
unzip 1.7.zip && \
rm 1.7.zip && \
mv ConsensuSV-core-1.7 ConsensuSV-core && \
cd ConsensuSV-core && \
unzip ALL_Illumina_Integrate_20170206.zip


RUN mamba install -y -c conda-forge -c bioconda -n base  biopython \
bwa=0.7.17 \
bcftools=1.16 \
covtobed=1.3.5 \
cnvkit=0.9.10 \
gsutil=5.17 \
fastp=0.23.2 \
fastqc=0.11.9 \
fonttools=4.38.0 \
gatk4=4.3.0.0 \
ghostscript=9.54.0 \
matplotlib=3.6.2 \
mosdepth=0.3.3 \
multiqc=1.13 \
nextflow=21.10.6 \
ensembl-vep=108.2 \
picard=2.27.4 \
samplot=1.3.0 \
r-ggplot2 \
r-gplots \
r-gsalib \
r-argparse \
r-curl \
openjdk==8.0.332=h166bdaf_0

RUN mamba install -y -c conda-forge -c bioconda -n base samtools
RUN pip install vatools==5.0.1
RUN pip install scikit-learn==1.2.2
