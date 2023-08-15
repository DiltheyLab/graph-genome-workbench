#!/bin/bash

snakemake -s download-data.smk --cores 27 
snakemake -s genotyping.smk --cores 27 -n
snakemake -s evaluation.smk --cores 27 -n
