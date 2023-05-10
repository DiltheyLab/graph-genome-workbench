configfile: "config.json"

pangenie = config['pangenie']

# Sample to be genotyped
sample = config['sample'] 

rule all:
	input:
		expand("genotyping/{sample}_genotyping.vcf.gz",  sample=sample)



################################################
# genotype a sample using PanGenie
################################################
rule genotyping:
	input:
		reads = config['reads'],
		reference = config['reference'],
		input_pangenome_graph = config['input_pangenome_graph']
	output:
		reads = temp("genotyping/{sample}-reads.fastq"),
		uncompressed_input_graph = temp("genotyping/Uncompressed_input_graph_for_{sample}.vcf"),
		genotypes = "genotyping/{sample}_genotyping.vcf"
	threads: 24
	log:
		"genotyping/{sample}.log"
	params:
		prefix = "genotyping/{sample}"
	resources:
		mem_total_mb=100000,
		runtime_hrs=4,
		runtime_min=59
	run:
		# jellyfish requires that the files are uncompressed
		shell("gunzip -c {input.reads} > {output.reads}")
		shell("gunzip -c {input.input_pangenome_graph} > {output.uncompressed_input_graph}")
		shell("(/usr/bin/time -v {pangenie} -i {output.reads} -v {output.uncompressed_input_graph} -r {input.reference} -o {params.prefix} -s {wildcards.sample} -j {threads} -t {threads} -g) &> {log}")


rule compress:
	input:
		"genotyping/{sample}_genotyping.vcf"
	output:
		"genotyping/{sample}_genotyping.vcf.gz"
	shell:
		"""
		bgzip -c {input} > {output}
		tabix -p vcf {output}
		"""
