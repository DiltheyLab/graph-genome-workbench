################## BASE IMAGE ######################
FROM biocontainers/biocontainers:v1.0.0_cv4


################## METADATA ######################
LABEL base_image="biocontainers:v1.0.0_cv4"\
      version="1"\
      about.summary="Genotyping software. Utilitzing a reference genome and known variants in a graph structure (a pangenome reference)"\
      about.home="https://github.com/DiltheyLab/graph-genome-workbench"\
      about.documentation="https://github.com/DiltheyLab/graph-genome-workbench"\
      about.license_file="https://github.com/DiltheyLab/graph-genome-workbench"\
      about.license="SPDX:GPL-3.0"\
      about.tags="Genomics"\
      extra.binaries="/usr/bin/graphtyper"\
      extra.binaries="/usr/bin/bayesTyper"\
      extra.binaries="/usr/bin/bayesTyperTools"

################## MAINTAINER ######################
MAINTAINER Torsten Houwaart <houwaart@hhu.du>

################## INSTALLATION ######################
USER root
WORKDIR /app
RUN wget https://github.com/bioinformatics-centre/BayesTyper/releases/download/v1.5/bayesTyper_v1.5_linux_x86_64.tar.gz
RUN tar -xzf bayesTyper_v1.5_linux_x86_64.tar.gz 
RUN wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.0-beta/graphtyper
RUN mv bayesTyper_v1.5_linux_x86_64/bin/bayesTyper /usr/bin/ && chmod a+x /usr/bin/bayesTyper
RUN mv bayesTyper_v1.5_linux_x86_64/bin/bayesTyperTools /usr/bin/ && chmod a+x /usr/bin/bayesTyperTools
RUN chmod a+x graphtyper &&  mv graphtyper /usr/bin/

USER biodocker


################## RUN ######################
#CMD ./graphtyper --help
