################## BASE IMAGE ######################
FROM biocontainers/biocontainers:v1.0.0_cv4


################## METADATA ######################
LABEL base_image="biocontainers:v1.0.0_cv4"
LABEL version="4"
LABEL software="graphtyper"
LABEL software.version="2.0beta"
LABEL about.summary="Genotyping software. Utilitzing a reference genome and known variants in a graph structure (a pangenome reference)"
LABEL about.home="https://github.com/DecodeGenetics/graphtyper"
LABEL about.documentation="https://github.com/DecodeGenetics/graphtyper"
LABEL about.license_file="https://github.com/DecodeGenetics/graphtyper"
LABEL about.license="SPDX:GPL-3.0"
LABEL about.tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Torsten Houwaart <houwaart@hhu.du>

################## INSTALLATION ######################
USER root
WORKDIR /app
RUN wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.0-beta/graphtyper
RUN chmod a+x graphtyper
USER biodocker


################## RUN ######################
CMD ./graphtyper --help
