# Downloading data, variant calling and graph construction for paper experiments

* downloading data: `` snakemake -s download-data.smk ``
* variant calling and graph construction for leave-one-out experiments: `` snakemake -s create-data.smk --configfile config.json --use-conda ``
* variant calling and graph construction for HLA experiments (inclusing MHC region): `` snakemake -s create-data.smk --configfile config-no-alt.json --use-conda ``

Note: Paths to assemblies must be set in config files prior to running.
