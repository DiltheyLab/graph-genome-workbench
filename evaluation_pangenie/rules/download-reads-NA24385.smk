numbers = ['00' + str(i) for i in range(1,10)] + ['0' + str(i) for i in range(10,18)]
input_reference = config['reference']
reads = config['reads']
sample = config['reads']

rule download_giab_fastq:
	output:
		"downloaded/reads/NA24385/raw/D1_S1_{l}_R{r}_{n}.fastq.gz"
	shell:
		"wget -O {output} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/D1_S1_{wildcards.l}_R{wildcards.r}_{wildcards.n}.fastq.gz"

rule combine_giab_fastq:
	input:
		expand("downloaded/reads/NA24385/raw/D1_S1_{l}_R{r}_{n}.fastq.gz", l = ['L001', 'L002'], r = ['1', '2'], n = numbers)
	output:
		"downloaded/reads/NA24385/raw_reads.fastq.gz"
	shell:
		"cat {input} > {output}"

rule uncompress_combined_reads:
    input:
        "downloaded/reads/NA24385/raw_reads.fastq.gz"
    output:
        "downloaded/reads/NA24385/raw_reads.fastq"
    shell:
        "gunzip {input} > {output}"

#### data for mapping based approaches ####

# index fasta
rule bwa_index:
	input:
		input_reference
	output:
		input_reference + ".ann"
	log:
		"downloaded/logs/reference-indexing.log"
	resources:
		mem_total_mb=5000
	shell:
		"(/usr/bin/time -v bwa index {input}) &> {log}"

# create fasta.fai file
rule samtools_faidx:
	input:
		input_reference
	output:
		input_reference + '.fai'
	shell:
		"samtools faidx {input}"

# align illumina reads
rule bwa_mem:
	input:
		reads=reads,
		fasta=input_reference,
		index=input_reference + '.ann',
		fai=input_reference + '.fai'
	output:
		'downloaded/BAM/aligned_reads.bam'
	log:
		'downloaded/logs/bwa_mem.log'
	threads: 24
	resources:
		mem_total_mb=60000,
		runtime_hrs=25,
		runtime_min=1
	shell:
		'(/usr/bin/time -v bwa mem -t {threads} -M {input.fasta} -R "@RG\\tID:{sample}\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:{sample}" {input.reads} | samtools view -bS | samtools sort -o {output} - ) &> {log}'

# index BAM file
rule samtools_index:
	input:
		"downloaded/BAM/aligned_reads.bam"
	output:
		"downloaded/BAM/aligned_reads.bam.bai"
	log:
		"downloaded/logs/indexing.log"
	shell:
		"(/usr/bin/time -v samtools index {input}) &> {log}"
