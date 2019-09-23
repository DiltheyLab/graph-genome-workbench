################## BASE IMAGE ######################
FROM biocontainers/biocontainers:v1.0.0_cv4


################## METADATA ######################
LABEL base_image="biocontainers:v1.0.0_cv4"
LABEL version="4"
LABEL software="graphtyper"
LABEL software.version="2.0beta"
LABEL about.summary="Genotyping software. Utilitzing a reference genome and known variants in a graph structure (a pangenome reference)"
LABEL about.home="https://github.com/bioinformatics-centre/BayesTyper"
LABEL about.documentation="https://github.com/bioinformatics-centre/BayesTyper"
LABEL about.license_file="https://github.com/bioinformatics-centre/BayesTyper"
LABEL about.license="SPDX:GPL-3.0"
LABEL about.tags="Genomics"
LABEL extra.binaries="/usr/bin/bayestyper"

################## MAINTAINER ######################
MAINTAINER Torsten Houwaart <houwaart@hhu.du>

################## INSTALLATION ######################
USER root
WORKDIR /app
RUN wget https://github.com/bioinformatics-centre/BayesTyper/releases/download/v1.5/bayesTyper_v1.5_linux_x86_64.tar.gz
RUN tar -xzf bayesTyper_v1.5_linux_x86_64.tar.gz 
RUN mv bayesTyper_v1.5_linux_x86_64 /usr/bin/bayestyper
RUN wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.0-beta/graphtyper
RUN chmod a+x graphtyper &&  mv graphtyper /usr/bin/

USER biodocker


################## RUN ######################
#CMD ./graphtyper --help
