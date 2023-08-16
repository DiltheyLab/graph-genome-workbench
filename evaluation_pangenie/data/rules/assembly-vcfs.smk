################################ Call variants from haplotype-resolved assemblies ################################
#
# steps:
#
#  1.) align contigs to reference using minimap2
#  2.) determine regions (uniquely) covered contigs haplotypes
#  3.) call variants from assemblies using paftools
#  4.) generate bi-allelic vcf file
#  5.) check mendelian consistency in trios and construct graph (multi-allelic vcf)
#  6.) compute some statistics
#
# output:
#
# multisample vcf file containing phased variants
#
##################################################################################################################


#configfile: "config-no-alt.json"

# all samples
#samples = config['dataset']['assemblies']['contigs'].keys()
#samples_parents = [s for s in samples if not s in config['dataset']['assemblies']['trios']]
samples_parents = config['samples_HGSVC']
samples = config['samples_HGSVC']

subset_to_samples = {
	'all' : samples_parents
}

for sample in config['dataset']['assemblies']['leave_out_sample']:
	skipped = [s for s in samples_parents if s != sample]
	subset_to_samples['only' + sample] = [sample]
	subset_to_samples['no' + sample] = skipped
 
scripts = config['scripts']
outdir = config['outdir']
chromosomes = ['chr' + str(i) for i in range(1,23)] + ['chrX']
frac_missing = 0.2 # skip positions with more than this fraction of missing alleles

# paftools skips contig-alignments shorter than this threshold
min_alignment_len = 50000


###########################################
#		1) Alignment
###########################################

# align assemblies to reference genome
rule align_assemblies_paf:
	input:
		contigs = lambda wildcards: config['dataset']['assemblies']['contigs'][wildcards.sample][int(wildcards.haplotype)],
		reference = config['reference']
	output:
		 outdir + "paf/{sample}-hap{haplotype}.paf"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=50000,
		runtime_hrs=3,
		runtime_min=1
	threads: 24
	shell:
		"minimap2 -cx asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5  --cs -t {threads} {input.reference} {input.contigs} | sort -k6,6 -k8,8n > {output}"


###########################################
#		2) Callable regions
###########################################

# align assemblies to reference and produce BAM output
rule align_assemblies_bam:
	input:
		contigs = lambda wildcards: config['dataset']['assemblies']['contigs'][wildcards.sample][int(wildcards.haplotype)],
		reference = config['reference']
	output:
		outdir + "bam/{sample}-hap{haplotype}.bam"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=50000,
		runtime_hrs=3,
		runtime_min=1
	log:
		outdir + "bam/{sample}-hap{haplotype}.log"
	threads: 24
	shell:
		"""
		minimap2 -ax asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 -t {threads} {input.reference} {input.contigs} | samtools view -bS | samtools sort -o {output} - &> {log}
		samtools index {output}
		"""


# compute regions covered by at least one contig
rule compute_covered_regions:
	input:
		outdir + "bam/{sample}-hap{haplotype}.bam"
	output:
		outdir + "bed/{sample}-hap{haplotype}_covered.bed"
	conda:
		"../../env/calling.yml"
	shell:
		"bedtools bamtobed -i {input} | awk '($3-$2) >= {min_alignment_len}' | bedtools merge > {output}"


# compute regions with per-base coverage < 2 
# NOTE: this will NOT remove cases in which there are more than one contig, but all except one contain a deletion)
# 	CCCCCCCCCCCC
# 	C----------C
rule compute_coverage:
	input:
		outdir + "bam/{sample}-hap{haplotype}.bam"
	output:
		outdir + "bed/{sample}-hap{haplotype}_unique.bed"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=10000,
		runtime_hrs=3,
		runtime_min=1
	shell:
		"bedtools genomecov -bga -ibam {input} | awk '$4 < 2' | bedtools merge > {output}"


rule intersect_beds:
	input:
		covered= outdir + "bed/{sample}-hap{haplotype}_covered.bed",
		unique= outdir + "bed/{sample}-hap{haplotype}_unique.bed"
	output:
		outdir + "bed/{sample}-hap{haplotype}_callable.bed"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=10000,
		runtime_hrs=3,
		runtime_min=1
	shell:
		"bedtools intersect -a {input.covered} -b {input.unique} > {output}"

rule merge_beds:
	input:
		expand(outdir + "bed/{sample}-hap{haplotype}_callable.bed", sample=samples_parents, haplotype=[0,1])
	output:
		outdir + "bed/callable-regions.bed"
	conda:
		"../../env/calling.yml"
	shell:
		"cat {input} | bedtools sort | bedtools merge > {output}"

rule sort_bed:
	input:
		"{filename}.bed"
	output:
		"{filename}-sorted.bed"
	conda:
		"../../env/calling.yml"
	shell:
		"bedtools sort -i {input} > {output}"

rule callable_genotyping:
	input:
		expand(outdir + "bed/{sample}-hap{haplotype}_callable-sorted.bed", sample=samples, haplotype=[0,1])
	output:
		outdir + "bed/callable-regions-missing.bed"
	conda:
		"../../env/calling.yml"
	params:
		covered = int((1-frac_missing) * len(samples)*2)
	shell:
		"bedtools multiinter -i {input} | awk '$4 > {params.covered}' | bedtools merge > {output}"


rule sample_callable:
	input:
		hap0=outdir + "bed/{sample}-hap0_callable.bed",
		hap1=outdir + "bed/{sample}-hap1_callable.bed"
	output:
		outdir + "bed/{sample}_callable.bed"
	wildcard_constraints:
		sample = "|".join(config['dataset']['assemblies']['leave_out_sample'])
	conda:
		"../../env/calling.yml"
	shell:
		"bedtools intersect -a {input.hap0} -b {input.hap1} > {output}"


##########################################
#		3) Variant Calling 
##########################################

# call variants from alignments
rule paftools:
	input:
		paf= outdir + "paf/{sample}-hap{haplotype}.paf",
		reference = config['reference']
	output:
		outdir + "calls/{sample}-hap{haplotype, [0,1]}.vcf"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=50000,
		runtime_hrs=3,
		runtime_min=1
	shell:
		"paftools.js call -L {min_alignment_len} -s {wildcards.sample}_{wildcards.haplotype} -f {input.reference} {input.paf} | sed 's|1/1|1|g' > {output}"


############################################################
#  4) Determine callable regions and create bi-allelic VCF
############################################################


rule compress_vcf:
	input:
		"{filename}.vcf"
	output:
		gz="{filename}.vcf.gz",
		tbi="{filename}.vcf.gz.tbi"
	conda:
		"../../env/calling.yml"
	shell:
		"""
		bgzip -c {input} > {output.gz}
		tabix -p vcf {output.gz}
		"""

# create a multisample VCF containing all haplotypes
rule collect_all_haplotypes:
	input:
		vcfs=expand("{outdir}calls/{sample}-hap{haplotype}.vcf.gz", outdir=outdir, sample=samples, haplotype=[0,1]),
		tbi=expand("{outdir}calls/{sample}-hap{haplotype}.vcf.gz.tbi", outdir=outdir, sample=samples, haplotype=[0,1])
	output:
		outdir + "calls/all-haplotypes.vcf"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=50000,
		runtime_hrs=3,
		runtime_min=1
	shell:
		"bcftools merge -m none --missing-to-ref {input.vcfs} | python3 {scripts}/assign-variant-ids.py > {output}"


# extract variant ids of callable regions
rule extract_covered_ids:
	input:
		bed = outdir + "bed/{sample}-hap{haplotype}_callable.bed",
		vcf = outdir + "calls/all-haplotypes.vcf"
	output:
		outdir + "bed/{sample}_{haplotype}.txt"
	conda:
		"../../env/calling.yml"
	shell:
		"bedtools intersect -a {input.vcf} -b {input.bed} -wa -f 1.0 | cut -f3 > {output}"


# set alleles outside of callable regions to missing
rule set_to_missing:
	input:
		vcf = outdir + "calls/all-haplotypes.vcf",
		bed = expand("{outdir}bed/{sample}_{haplotype}.txt", outdir=outdir, sample=samples, haplotype=[0,1])
	output:
		outdir + "calls/all-haplotypes-callable.vcf"
	log:
		outdir + "calls/all-haplotypes-callable.log"
	resources:
		mem_mb=150000,
		runtime_hrs=10,
		runtime_min=1
	shell:
		"python3 {scripts}/set-to-missing.py -v {input.vcf} -m {frac_missing} -f {input.bed} 2> {log} 1> {output}"


# convert haploid VCF into a diploid one by combining haplotypes of each sample
rule write_input:
	output:
		outdir + "samples.txt"
	run:
		with open(output[0], 'w') as txt_output:
			for sample in samples:
				txt_output.write('\t'.join([sample, sample + '_0', sample + '_1']) + '\n')


rule combine_haplotypes:
	input:
		haps = outdir + "calls/all-haplotypes-callable.vcf",
		samples = outdir + "samples.txt"
	output:
		outdir + "multisample-vcfs/assemblies-all-samples-biallelic.vcf"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=100000,
		runtime_hrs=10,
		runtime_min=1
	shell:
		'python3 {scripts}/merge_vcfs.py combine_columns -samples {input.samples} -vcf {input.haps} > {output}'



#################################################################
#  5) Check mendelian consistency for trios and construct graph
#################################################################

# generate a file specifying the trio relationships
rule generate_ped_file:
	output:
		"{outdir}trios.ped"
	run:
		with open(output[0], "w") as ped_output:
			for trio in config['dataset']['assemblies']['trios']:
				father=config['dataset']['assemblies']['trios'][trio][0]
				mother=config['dataset']['assemblies']['trios'][trio][1]
				ped_output.write('\t'.join([trio, trio, father, mother]) + '\n')
				
rule generate_samples_file:
	output:
		"{outdir}trio-samples.txt"
	run:
		with open(output[0], "w") as sample_output:
			for trio in config['dataset']['assemblies']['trios']:
				sample_output.write(trio + '\n')
				for sample in config['dataset']['assemblies']['trios'][trio]:
					sample_output.write(sample + '\n')

rule generate_sample_list:
	output:
		"{outdir}sex-samples.tsv"
	run:
		with open(output[0], 'w') as sample_output:
			for sample in config['dataset']['assemblies']['sex']:
				sample_output.write(sample + '\t' + config['dataset']['assemblies']['sex'][sample] + '\n')


# remove all variants where there is a mendelian conflict in at least one of the trios
# if no trios are given in config, the vcf does not change.
rule check_mendelian_consistency:
	input:
		vcf="{outdir}multisample-vcfs/assemblies-all-samples-biallelic.vcf",
		ped="{outdir}trios.ped",
		samples="{outdir}trio-samples.txt",
		sexes = "{outdir}sex-samples.tsv"
	output:
		vcf="{outdir}multisample-vcfs/assemblies-all-samples-biallelic-filtered.vcf",
		tmp="{outdir}multisample-vcfs/assemblies-all-samples-biallelic-mendel-filtered.vcf",
		tsv="{outdir}multisample-vcfs/assemblies-all-samples-biallelic-filtered.tsv"
	log:
		mendel="{outdir}multisample-vcfs/assemblies-all-samples-filtered.log",
		sex_chrom="{outdir}multisample-vcfs/assemblies-all-samples-filtered-sex.log"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=100000,
		runtime_hrs=10,
		runtime_min=1
	shell:
		"""
		python3 {scripts}/mendelian-consistency.py filter -vcf {input.vcf} -samples {input.samples} -ped {input.ped} -o {output.tsv} --remove-children 2> {log.mendel}  1> {output.tmp}
		cat {output.tmp} | python3 {scripts}/filter-sex-chromosomes.py {input.sexes} 2> {log.sex_chrom} 1> {output.vcf} 
		"""


rule merge_haplotypes:
	input:
		vcf = outdir + "multisample-vcfs/assemblies-all-samples-biallelic-filtered.vcf.gz",
		reference = config['reference']
	output:
		biallelic=outdir + "multisample-vcfs/assemblies-{subset}-biallelic-filtered.vcf",
		vcf=outdir + "multisample-vcfs/assemblies-{subset}-filtered.vcf"
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	params:
		chrom = ','.join([c for c in chromosomes]),
		samples = lambda wildcards : ','.join([s for s in subset_to_samples[wildcards.subset]])
	log:
		outdir + "multisample-vcfs/assemblies-{subset}-filtered.log"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=900000,
		runtime_hrs=30,
		runtime_min=1
	shell:
		"""
		bcftools view -s {params.samples} {input.vcf} | bcftools view --min-ac 1 > {output.biallelic}
		python3 {scripts}/merge_vcfs.py merge -vcf {output.biallelic} -r {input.reference} -ploidy 2 -chromosomes {params.chrom} 2> {log} 1> {output.vcf}
		"""



#############################################
#  6) Generate statistics and create plots 
#############################################

### Here there is an "error". I added to the command -c x to exclude bad sites, otherwise it was not working (see error in screenshot error_normalizing_input_graph). However, the resultant VCF file contains unphased genotypes. We need to either exclude those rows or just transform them to phased genotypes. ---> See genotyping.smk. 
### There are some differences in the ordering of phased haplotypes when using the -c x option and without. See screenshot: differences_normalized_VCFs. I'll go with the approach WITHOUT -c x, however, to do so the input file must contain in the header all the chromosomes! So add it manually (maybe for future also in previous rule, i.e. merge haplotypes, but so far it would take too much time). So erase all .normalized.vcf files, inlucde the necessary header into filtered.vcf.gz and recreate them. 

### SOLUTION: avoid the bcftools from the conda file. It's probably an old version... It works without the -c x option!!

rule normalize_vcf:
    input:
        vcf = outdir + "multisample-vcfs/assemblies-{subset}-filtered.vcf.gz",
        tbi = outdir + "multisample-vcfs/assemblies-{subset}-filtered.vcf.gz.tbi",
        reference = config['reference']
    output:
        vcf = outdir + "multisample-vcfs/assemblies-{subset}-filtered.normalized.vcf"
    wildcard_constraints:
        subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
    conda:
        "../../env/calling.yml"
    log: 
        outdir + "multisample-vcfs/log-{subset}-filtered.normalized.log"
    shell:
        """
        (/usr/bin/time -v bcftools norm -m +any -d all -f {input.reference} {input.vcf} | bcftools sort > {output.vcf}) &> {log}
        """



rule untypable_ids:
	input:
		config["untypable_input"]
		#outdir + "multisample-vcfs/assemblies-all-samples-biallelic-filtered.vcf.gz"
	output:
		lists=expand(outdir + "statistics/untypable-ids/{sample}-untypable.tsv", sample=samples_parents),
		summary= outdir + "statistics/untypable-ids.tsv"
	params:
		out= outdir + "statistics/untypable-ids"
	conda:
		"../../env/calling.yml"
	shell:
		"zcat {input} | python3 {scripts}/untypable-ids.py {params.out} > {output.summary}"


rule untypable_bubbles:
	input:
		config["input_graph"]
		#outdir + "multisample-vcfs/assemblies-all-filtered.normalized.vcf"
	output:
		table=outdir + "statistics/untypable-bubbles.tsv",
		lists=expand(outdir + "statistics/untypable-bubbles/{sample}-untypable-bubble.bed", sample=samples_parents)
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	log:
		outdir + "statistics/untypable-bubbles.log"
	conda:
		"../../env/calling.yml"
	resources:
		runtime_hrs=4
	params:
		out= outdir + 'statistics'
	shell:
		"python3 {scripts}/untypable-bubbles.py {input} {params.out} &> {log}"

rule bcftools_statistics:
	input:
#		outdir + "multisample-vcfs/assemblies-{subset}-filtered.normalized.vcf"
		outdir + "multisample-vcfs/assemblies-{subset}-biallelic-filtered.vcf"
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	output:
		txt= outdir +"statistics/bcftools-{subset}.txt",
		plots= outdir +"statistics/bcftools-plots/summary-{subset}.pdf"
	conda:
		"../../env/calling.yml"
	shell:
		"""
		bcftools stats {input} > {output.txt}
		plot-vcfstats -p {outdir}statistics/bcftools-plots {output.txt}
		"""

rule vcfstats_statistics:
	input:
#		outdir + "multisample-vcfs/assemblies-{subset}-filtered.normalized.vcf"
		outdir + "multisample-vcfs/assemblies-{subset}-biallelic-filtered.vcf"
	output:
		txt= outdir + "statistics/vcfstats-{subset}.txt"
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	conda:
		"../../env/calling.yml"
	shell:
		"/home/ubuntu/rtg-tools-3.12.1/rtg vcfstats {input} > {output}"

rule vcfstats_plot_het_hom:
	input:
		outdir + "statistics/vcfstats-{subset}.txt"
	output:
		outdir + "statistics/vcfstats-plots/het-hom-{subset}.pdf"
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	conda:
		"../../env/calling.yml"
	shell:
		"python3 {scripts}/plot-callset-statistics.py ratios {input} {output}"

rule vcfstats_plot_numbers:
	input:
		outdir + "statistics/vcfstats-{subset}.txt"
	output:
		outdir + "statistics/vcfstats-plots/variant-numbers-{subset}.pdf"
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	conda:
		"../../env/calling.yml"
	shell:
		"python3 {scripts}/plot-callset-statistics.py numbers {input} {output}"

rule indel_histogram:
	input:
#		outdir + "multisample-vcfs/assemblies-{subset}-filtered.normalized.vcf"
		outdir + "multisample-vcfs/assemblies-{subset}-biallelic-filtered.vcf"
	output:
		histo= outdir + "statistics/vcftools-{subset}.indel.hist",
		plot= outdir + "statistics/vcftools-plots/indel-histogram-{subset}.pdf"
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	conda:
		"../../env/calling.yml"
	shell:
		"""
		vcftools --vcf {input} --out {outdir}statistics/vcftools-{wildcards.subset} --hist-indel-len
		cat {output.histo} | python3 {scripts}/plot-callset-statistics.py length {output.plot} 20000
		"""

###################### collect statistics #########################


# for all subsets total number of alleles detected across all haplotypes
rule paftools_statistics_all:
	input:
		vcf=lambda wildcards : expand(outdir + "calls/{sample}-hap{haplotype}.vcf.gz", sample=[s for s in subset_to_samples[wildcards.subset] if s not in config['dataset']['assemblies']['trios'].keys()], haplotype=[0,1])
	output:
		vcf=outdir + "statistics/raw-callset-{subset}.vcf",
		stats= outdir + "statistics/raw-{subset}.tsv"
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	conda:
		"../../env/calling.yml"
	shell:
		"""
		bcftools merge -m none {input.vcf} > {output.vcf}
		python3 {scripts}/variant-counter.py {output.vcf} {output.stats}
		"""

# total number of alleles in callable regions
rule paftools_statistics_callable:
	input:
		outdir + "multisample-vcfs/assemblies-all-samples-biallelic.vcf.gz"
	output:
		vcf=outdir + "statistics/raw-callset-callable-{subset}.vcf",
		stats=outdir + "statistics/raw-callable-{subset}.tsv"
	conda:
		"../../env/calling.yml"
	params:
		samples = lambda wildcards: ','.join([s for s in subset_to_samples[wildcards.subset] if s not in config['dataset']['assemblies']['trios'].keys()])
	shell:
		"""
		bcftools view -s {params.samples} {input} | bcftools view --min-ac 1 > {output.vcf}
		python3 {scripts}/variant-counter.py {output.vcf} {output.stats}
		"""


# total number of alleles after mendelian filtering
rule merged_statistics_all:
	input:
		 outdir + "multisample-vcfs/assemblies-{subset}-biallelic-filtered.vcf"
	output:
		stats=outdir + "statistics/filtered-{subset}.tsv"
	conda:
		"../../env/calling.yml"
	shell:
		"""
		python3 {scripts}/variant-counter.py {input} {output.stats}
		"""

# total number of bubbles after mendelian filtering
rule merged_statistics_all_filtered:
	input:
		outdir + "multisample-vcfs/assemblies-{subset}-filtered.normalized.vcf"
	output:
		stats= outdir + "statistics/merged-filtered-{subset}.tsv"
	wildcard_constraints:
		subset = "|".join([s for s in subset_to_samples.keys()] + ['all'])
	params:
		samples = lambda wildcards: ','.join([s for s in subset_to_samples[wildcards.subset] if s not in config['dataset']['assemblies']['trios'].keys()])
	conda:
		"../../env/calling.yml"
	shell:
		"""
		python3 {scripts}/variant-counter.py {input} {output.stats}
		"""

rule alleles_per_bubble:
	input:
		config['input_graph'],
		#outdir + "multisample-vcfs/assemblies-{subset}-filtered.normalized.vcf"
	output:
		plot=outdir + "statistics/alleles-per-bubble-{subset}.pdf",
		bed=outdir + "statistics/complex-bubbles-{subset}.bed"
	conda:
		"../../env/calling.yml"
	resources:
		mem_mb=5000,
		runtime_hrs=1
	shell:
		"cat {input} | python3 {scripts}/variant-statistics.py {output.plot} 1 > {output.bed}"
