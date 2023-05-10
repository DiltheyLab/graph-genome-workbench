configfile: "config.json"

include: "rules/download-references.smk"
include: "rules/download-reads-NA24385.smk"
include: "rules/download-pangenome-graphs.smk"
include: "rules/download-giab-variants.smk"

sample = config['sample']

rule all:
	input:
		# references
		'downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa',
		'downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai',
		
		# reads
		expand("downloaded/reads/{sample}/raw_reads.fastq.gz", sample=sample),
		
		# Input Pangenome Graph (PanGenie Github)
		'downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz',
		'downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze4_64haplotypes.vcf.gz',
		'downloaded/vcf/HPRC-GRCh38/Pangenome_graph_88haplotypes.vcf.gz',
		'downloaded/vcf/HPRC-CHM13/Pangenome_graph_88haplotypes.vcf.gz', 

		# Benchmark set HG002_NA24385_son
		'downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz', 
		'downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi', 
		'downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.bed' 
