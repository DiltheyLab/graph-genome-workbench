configfile: "config.json"

input_reference=config['data']['reference']
reference_version=config['data']['reference_version']
outname=config['parameters']['outname_reads']

prefix = "" if "37" in reference_version or "19" in reference_version else "chr"


##############################################################
#################   prepare read data     ##################
##############################################################

###compute coverage of full data ####
#rule compute_bam_coverage:
#	input:
#		full_cov_bam='{results}/{sample}/aligned/{sample}-full.bam',
#		full_cov_bai='{results}/{sample}/aligned/{sample}-full.bam.bai'
#	output:
#		"{results}/{sample}/raw/{sample}-coverage.cov"
#	conda:
#		'../env/genotyping.yml'
#	resources:
#		mem_total_mb=10000,
#		runtime_hrs=5,
#		runtime_min=1
#	log:
#		"{results}/{sample}/raw/{sample}-coverage.log"
#	shell:
#		"bash ../scripts/compute-coverage.sh {input.full_cov_bam} {output} &> {log}"
		

### downsample reads ####
# rule downsample_reads:
#	input:
#		reads1=lambda wildcards: config['data'][wildcards.sample]['reads'][0],
#		reads2=lambda wildcards: config['data'][wildcards.sample]['reads'][1],
#		coverage="{results}/{sample}/raw/{sample}-coverage.cov"
#	output:
#		sampled1=temp("{results}/{sample}/raw/{sample}-{fraction, [0-9.]+}_1.fastq"),
#		sampled2=temp("{results}/{sample}/raw/{sample}-{fraction, [0-9.]+}_2.fastq")
#	conda:
#		'../env/genotyping.yml'
#	resources:
#		mem_total_mb=20000,
#		runtime_hrs=5,
#		runtime_min=1
#	log:
#		"{results}/{sample}/raw/{sample}-{fraction, [0-9.]+}.downsample.log"
#	shell:
#		"bash ../scripts/downsample-fasta.sh {input.coverage} {wildcards.fraction} {input.reads1} {input.reads2} {output.sampled1} {output.sampled2} &> {log}"


### data for mapping free approaches ####

#def combine_reads_input(wildcards):
#	if wildcards.fraction == "full":
#		return config['data'][wildcards.sample]['reads']
#	else:
#		return expand("{results}/{sample}/raw/{sample}-{fraction}_{r}.fastq", results=wildcards.results, sample=wildcards.sample, fraction=wildcards.fraction, r=[1,2])

## generate combined fastq file
#rule combine_reads:
#	input:
#		combine_reads_input
#	output:
#		"{results}/{sample}/raw/{sample}-{fraction, (full|[0-9.]+)}.fastq"
#	shell:
#		"cat {input} > {output}"


### data for mapping based approaches ####

## index fasta
rule bwa_index:
	input:
		input_reference
	output:
		input_reference + ".ann"
	log:
		"{results}/reference-indexing.log".format(results=outname)
	conda:
		'../env/genotyping.yml'
	resources:
		mem_total_mb=5000
	shell:
		"(/usr/bin/time -v bwa index {input}) &> {log}"

## create fasta.fai file
rule samtools_faidx:
	input:
		input_reference
	output:
		input_reference + '.fai'
	conda:
		'../env/genotyping.yml'
	shell:
		"samtools faidx {input}"

## align illumina reads
#def bwa_mem_input(wildcards):
#	if wildcards.fraction == 'full':
#		return [config['data'][wildcards.sample]['reads'][0], config['data'][wildcards.sample]['reads'][1]]
#	else:
#		return expand("{results}/{sample}/raw/{sample}-{fraction}_{r}.fastq", results=wildcards.results, sample=wildcards.sample, fraction=wildcards.fraction, r=[1,2])

rule bwa_mem:
	input:
		reads=lambda wildcards: config['data'][wildcards.sample]['reads'],
		fasta=input_reference,
		index=input_reference + '.ann',
		fai=input_reference + '.fai'
		#reads=bwa_mem_input,
	output:
		'{results}/{sample}/aligned/{sample}-{fraction, (full|[0-9.]+)}.bam'
	log:
		'{results}/{sample}/aligned/{sample}-{fraction, (full|[0-9.]+)}.log'
	threads: 24
	resources:
		mem_total_mb=60000,
		runtime_hrs=25,
		runtime_min=1
	conda:
		'../env/genotyping.yml'
	shell:
		'(/usr/bin/time -v bwa mem -t {threads} -M {input.fasta} -R "@RG\\tID:{wildcards.sample}\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:{wildcards.sample}" {input.reads} | samtools view -bS | samtools sort -o {output} - ) &> {log}'

## index BAM file
rule samtools_index:
	input:
		"{filename}.bam"
	output:
		"{filename}.bam.bai"
	log:
		"{filename}-index.log"
	conda:
		'../env/genotyping.yml'
	shell:
		"(/usr/bin/time -v samtools index {input}) &> {log}"

## split BAM by chromosome
rule split_bam_by_chromosome:
	input:
		bam='{results}/{sample}/aligned/{sample}-{fraction, [0-9.]+}.bam',
		bai='{results}/{sample}/aligned/{sample}-{fraction, [0-9.]+}.bam.bai'
	output:
		bam='{results}/{sample}/aligned/{sample}-{fraction}.chr{chrom, X|Y|[0-9]+}.bam'
	conda:
		'../env/genotyping.yml'
	shell:
		"""
        samtools view -h {input.bam} {prefix}{wildcards.chrom} | samtools view -Sb -> {output.bam}
        """

#############################################################
##### old execution with copying BAM files from graph-genome-workbench
#############################################################

#### My rule to copy data from graph-genome-workbench
#rule copying_files:
#    input: 
#        bam = "/vol/whopper/graph-genome-workbench/evaluation_pangenie/downloaded/BAM/aligned_reads.bam",
#        bai = "/vol/whopper/graph-genome-workbench/evaluation_pangenie/downloaded/BAM/aligned_reads.bam.bai",
#        fastq_gz = "/vol/whopper/graph-genome-workbench/evaluation_pangenie/downloaded/reads/NA24385/raw_reads.fastq.gz"
#    output:
#        bam = "reads/NA24385/aligned/NA24385-full.bam",
#        bai = "reads/NA24385/aligned/NA24385-full.bam.bai",
#        fastq = "reads/NA24385/raw/NA24385-full.fastq", 
#        fastq_gz = "reads/NA24385/raw/NA24385-full.fastq.gz" 
#    shell:
#        """
#        cp {input.bam} {output.bam}
#        cp {input.bai} {output.bai}
#        cp {input.fastq_gz} {output.fastq_gz}
#        gunzip -c {output.fastq_gz} > {output.fastq}
#        """


## split BAM by chromosome
#rule split_bam_by_chromosome:
#	input:
#		bam='{results}/{sample}/aligned/{sample}-{fraction, [0-9.]+}.bam',
#		bai='{results}/{sample}/aligned/{sample}-{fraction, [0-9.]+}.bam.bai'
#	output:
#		bam='{results}/{sample}/aligned/{sample}-{fraction}.chr{chrom, X|Y|[0-9]+}.bam',
#		bai='{results}/{sample}/aligned/{sample}-{fraction}.chr{chrom, X|Y|[0-9]+}.bam.bai' 
#	conda:
#		'../env/genotyping.yml'
#	shell:
#		"""
#        samtools view -h {input.bam} {prefix}{wildcards.chrom} | samtools view -Sb -> {output.bam}
#        samtools index {output.bam} > {output.bai}
#        """
