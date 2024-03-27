tmp_workdir = config['utils']['tmp_workdir_computations']
bwa=config['programs']['bwa']
samtools=config['programs']['samtools']

downsampling_reads = {}
for line in open(config['reads'], 'r'):
        if line.startswith('#'):
                continue
        fields = line.strip().split()
        downsampling_reads[fields[1]] = fields[7] 


# align the reads to the reference to estimate the coverage of the full data set
# need to do it per callset because the reference might be different
rule align_reads:
	input:
		reads = lambda wildcards: downsampling_reads[wildcards.sample],
		fasta = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		fasta_ann = lambda wildcards: config['callsets'][wildcards.callset]['reference'] + ".ann"
	output:
		bam = "preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam"
		#unsorted_bam = "preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full_unsorted.bam",
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb = 60000,
		runtime_hrs = 25,
		runtime_min = 1
	threads: 28
	params:
		tmp_workdir = tmp_workdir 
	log:
		sorted_mem="logs/preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full_mem.log"
		#bwa_mem="logs/preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full_mem_unsorted.log",
	shell:
		"""
		(/usr/bin/time -v {bwa} mem -t {threads} -M {input.fasta} -R "@RG\\tID:{wildcards.sample}\\tLB:lib1\\tPL:illumina\\tPU:unit1\\tSM:{wildcards.sample}" {input.reads} | {samtools} view -bS | {samtools} sort -T {params.tmp_workdir} -o {output.bam} - ) &> {log.sorted_mem}
		"""
		#(/usr/bin/time -v  {output.unsorted_bam} > {output.bam} ) &> {log.sorted_mem}

## split BAM by chromosome
rule split_bam_by_chromosome:
	input:
		bam='preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam',
		bai='preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam.bai'
	output:
		bam='preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.chr{chrom}.bam'
	conda:
		'../envs/genotyping.yml'
	shell:
		"""
        {samtools} view -h {input.bam} chr{wildcards.chrom} | {samtools} view -Sb - > {output.bam}
        """

## index BAM file
rule samtools_index:
	input:
		"preprocessing/downsampling/{callset}/{coverage}/aligned/{filename}.bam"
	output:
		"preprocessing/downsampling/{callset}/{coverage}/aligned/{filename}.bam.bai"
	log:
		"logs/preprocessing/downsampling/{callset}/{coverage}/aligned/{filename}-index.log"
	shell:
		"(/usr/bin/time -v {samtools} index {input}) &> {log}"


# estimate the coverage of the aligned data
rule compute_bam_coverage:
	input:
		full_cov_bam = "preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam",
		full_cov_bai = "preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.bam.bai"
	output:
		"preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.cov"
	conda:
		'../env/genotyping.yml'
	resources:
		mem_total_mb = 10000,
		runtime_hrs = 5,
		runtime_min = 1
	log:
		"logs/preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full_cov.log"
	shell:
		"bash workflow/scripts/compute-coverage.sh {input.full_cov_bam} {output} &> {log}"


# downsample fastq to desired coverage
rule downsample_reads:
	input:
		lambda wildcards: downsampling_reads[wildcards.sample]
	output:
		"preprocessing/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 5,
		runtime_min = 1
	log:
		"logs/preprocessing/downsampling/{callset}/{coverage}/{sample}_{coverage}.log"
	shell:
		"bash workflow/scripts/downsample-fasta.sh {input.coverage} {wildcards.fraction} {input} {output} &> {log}"
