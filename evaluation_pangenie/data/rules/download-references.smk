# get reference sequence (hg38)
rule download_reference_hg38:
	output:
		fasta="downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa"
	shell:
		"wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

rule remove_alt:
	input:
		fa="downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa",
		fai="downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
	output:
		"downloaded/fasta/no-alt/no-alt.fa"
	shell:
		"samtools faidx {input.fa} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > {output}"


# index hg38 fasta
rule index:
	input:
		"{filename}.fa"
	output:
		"{filename}.fa.fai"
	shell:
		"samtools faidx {input}"



