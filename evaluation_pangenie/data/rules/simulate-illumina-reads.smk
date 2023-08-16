configfile:"config.json"

# parameters
reference=config['reference']
coverage=int(config['coverage'])
chromosomes=[i for i in range(1,23)] + ['X', 'Y']


############################# generate haplotype sequences ###################################


# genomesimulator script cannot handle multi-sample VCFs. By splitting the VCF by sample, we get biallelic VCFs,
# since the VCFs used to create this VCFs were biallelic as well (so no 1|2, 2|1 genotypes exist).
rule prepare_vcfs:
	input:
		lambda wildcards: config['dataset'][wildcards.dataset]['merged'] 
	output:
		'generated/reads/{dataset}/{sample}/haplotypes/input-variants.{sample}.vcf'
	shell:
		'bcftools view -s {wildcards.sample} -a {input} > {output}'

rule generate_haplotypes:
	input:
		vcf='generated/reads/{dataset}/{sample}/haplotypes/input-variants.{sample}.vcf',
		ref=reference
	output:
		hap1='generated/reads/{dataset}/{sample}/haplotypes/{sample}.chr{chromosome}.hap1.fasta',
		hap2='generated/reads/{dataset}/{sample}/haplotypes/{sample}.chr{chromosome}.hap2.fasta'
	log:
		'generated/reads/{dataset}/{sample}/haplotypes/{sample}.{chromosome}.log'
	params:
		scripts=config['scripts']
	run:
		shell('mkdir -p generated/reads/{wildcards.dataset}/{wildcards.sample}/haplotypes/tmp/{wildcards.sample}-{wildcards.chromosome}')
		shell('python3 {params.scripts}/genomesimulator.py -c {wildcards.chromosome} {input.vcf} {input.ref} generated/reads/{wildcards.dataset}/{wildcards.sample}/haplotypes/tmp/{wildcards.sample}-{wildcards.chromosome} > {log} 2>&1')
		shell('mv generated/reads/{wildcards.dataset}/{wildcards.sample}/haplotypes/tmp/{wildcards.sample}-{wildcards.chromosome}/{wildcards.sample}.chr{wildcards.chromosome}.1.fasta {output.hap1}')
		shell('mv generated/reads/{wildcards.dataset}/{wildcards.sample}/haplotypes/tmp/{wildcards.sample}-{wildcards.chromosome}/{wildcards.sample}.chr{wildcards.chromosome}.2.fasta {output.hap2}')

############################# simulate illumina reads ###################################

rule simulate_illumina_reads:
	input:
		'generated/reads/{dataset}/{sample}/haplotypes/{sample}.chr{chromosome, X|Y|[0-9]+}.hap{hap, [0-9]+}.fasta'
	output:
		fastq1='generated/reads/{dataset}/{sample}/simulation/{sample}.chr{chromosome, X|Y|[0-9]+}.hap{hap, [0-9]+}.1.tagged.fasta',
		fastq2='generated/reads/{dataset}/{sample}/simulation/{sample}.chr{chromosome, X|Y|[0-9]+}.hap{hap, [0-9]+}.2.tagged.fasta',
		aln1='generated/reads/{dataset}/{sample}/simulation/{sample}.chr{chromosome, X|Y|[0-9]+}.hap{hap, [0-9]+}.1.aln',
		aln2='generated/reads/{dataset}/{sample}/simulation/{sample}.chr{chromosome, X|Y|[0-9]+}.hap{hap, [0-9]+}.2.aln'
	log:
		'generated/reads/{dataset}/{sample}/simulation/{sample}.chr{chromosome}.hap{hap, [0-9]+}.log'
	params:
		cov=coverage/2,
		mean_is=750,
		mean_sd=50
	run:
		print('Simulating reads with parameters: cov={cov}, mean_is={mean_i}, mean_sd={mean_s}'.format(cov=params.cov, mean_i=params.mean_is, mean_s=params.mean_sd))
		shell('art_illumina -ss MSv3 -i {input} -p -l 250 -f {params.cov} -m {params.mean_is} -s {params.mean_sd} -o generated/reads/{wildcards.dataset}/{wildcards.sample}/simulation/{wildcards.sample}.chr{wildcards.chromosome}.hap{wildcards.hap}. &> {log}')
		shell('awk \'{{ gsub("chr{wildcards.chromosome}-", "chr{wildcards.chromosome}_HAP{wildcards.hap}-", $0); print $0 }}\' generated/reads/{wildcards.dataset}/{wildcards.sample}/simulation/{wildcards.sample}.chr{wildcards.chromosome}.hap{wildcards.hap}.1.fq > {output.fastq1}')
		shell('awk \'{{ gsub("chr{wildcards.chromosome}-", "chr{wildcards.chromosome}_HAP{wildcards.hap}-", $0); print $0 }}\' generated/reads/{wildcards.dataset}/{wildcards.sample}/simulation/{wildcards.sample}.chr{wildcards.chromosome}.hap{wildcards.hap}.2.fq > {output.fastq2}')
		shell('rm generated/reads/{wildcards.dataset}/{wildcards.sample}/simulation/{wildcards.sample}.chr{wildcards.chromosome}.hap{wildcards.hap}.1.fq')
		shell('rm generated/reads/{wildcards.dataset}/{wildcards.sample}/simulation/{wildcards.sample}.chr{wildcards.chromosome}.hap{wildcards.hap}.2.fq')

rule merge_chromosome:
	input:
		expand('generated/reads/{{dataset}}/{{sample}}/simulation/{{sample}}.chr{{chromosome}}.hap{hap}.{{number}}.tagged.fasta', hap=[1,2])
	output:
		'generated/reads/{dataset}/{sample}/reads-per-chromosome/{sample}.chr{chromosome, X|Y|[0-9]+}.{number, [0-9]+}.tagged.fasta'
	shell:
		'cat {input} > {output}'

rule merge_genome:
	input:
		expand('generated/reads/{{dataset}}/{{sample}}/reads-per-chromosome/{{sample}}.chr{chromosome}.{{number}}.tagged.fasta', chromosome=chromosomes)
	output:
		'generated/reads/{dataset}/{sample}/reads/{sample}.{number, [0-9]+}.tagged.fasta'
	shell:
		'cat {input} > {output}'
			

