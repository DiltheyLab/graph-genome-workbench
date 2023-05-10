This is a demonstration of how to use the PanGenie tool in order to evaluate the resulting variant-call output of PanGenie. The pipeline is applied to a sample (i.e. NA24385/HG002).
The pipeline is written in Snakemake and has three steps:

- `download-data.smk`: This downloads all the necessary data. This includes: short-reads of our sample, a reference genome, a pangenome graph in VCF format (used to represent variants of already known haplotypes) and a benchmark dataset in VCF format (used to evaluate the resulting genotyped data).

- `genotyping.smk`: This applies the PanGenie tool to the short-reads of our sample, the reference genome and the pangenome graph. As output, it produces a VCF file containing genotypes for the variants provided in the input pangenome graph.

- `evaluation`: With the help of VCFeval, this evaluates the resulting genotyped VCF file for our sample by comparing it to the corresponding assembly of our sample in VCF format. As output, it yields the measured performance between both VCF files mainly in terms of: precision, recall and F-measure.

Additionally, the `config.json` contains the paths to store our data and `run_pipeline.sh` allows us to run our workflow. To do so, execute the command: `sh run_pipeline.sh`