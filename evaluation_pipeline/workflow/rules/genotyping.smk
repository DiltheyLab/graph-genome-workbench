chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
ref_chromosomes = {}
for c in config['callsets'].keys():
    ref_chromosomes[c] = []
    for line in open(config['callsets'][c]['reference_fai'], 'r'):
        fields = line.split()
        if fields[0] in ['chr' + i for i in chromosomes]:
            chr_number = fields[0].split('chr')[1].strip()
            ref_chromosomes[c].append(chr_number)

bayestyper_reference_canon=config['utils']['bayestyper_reference_canon']
bayestyper_reference_decoy=config['utils']['bayestyper_reference_decoy']

# stores paths to reads
reads_leave_one_out = {}
sex_per_sample = {}

for line in open(config['reads'], 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	sample_name = fields[1]
	sample_sex = int(fields[4])
	read_path = fields[7]
	sex_per_sample[sample_name] = 'M' if sample_sex == 1 else 'F'
	reads_leave_one_out[sample_name] = read_path


########################################################
##################    Util functions    ################
########################################################
rule compress_vcf:
	input:
		"preprocessing/{filename}.vcf"
	output:
		vcf = "preprocessing/{filename}.vcf.gz",
		tbi = "preprocessing/{filename}.vcf.gz.tbi"
	priority: 1
	shell:
		"""
		bgzip -c {input} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

########################################################
##################    run PanGenie    ##################
########################################################

rule pangenie:
	input:
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "preprocessing/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		vcf="genotyping/{callset}/{sample}/input-panel/panel_multi_norm.vcf"
	output:
		genotyping = temp("genotyping/{callset}/{sample}/{version}/{coverage}/temp/{version}-{sample}_genotyping.vcf")
	log:
		"logs/genotyping/{callset}/{sample}/{version}/{coverage}/{version}-{sample}.log"
	threads: 24
	resources:
		mem_total_mb=80000,
		runtime_hrs=5,
		runtime_min=1
	priority: 1
	params:
		out_prefix="genotyping/{callset}/{sample}/{version}/{coverage}/temp/{version}-{sample}",
		pangenie = lambda wildcards: config['pangenie'][wildcards.version]
	wildcard_constraints:
		version = "|".join([k for k in config['pangenie'].keys()] + ['^' + k for k in config['pangenie-modules']])
	shell:
		"""
		(/usr/bin/time -v {params.pangenie} -i <(zcat {input.reads}) -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g ) &> {log}
		"""


# run pangenie in the modularized way (> v2.1.1)
rule pangenie_modules:
	input:
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "preprocessing/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		vcf = "preprocessing/{callset}/{sample}/input-panel/panel_multi_norm.vcf"
	output:
		genotyping = temp("genotyping/{callset}/{sample}/{version}/{coverage}/genotyping.vcf")
	log:
		index = "logs/genotyping/{callset}/{sample}/{version}/{coverage}/index.log",
		genotype = "logs/genotyping/{callset}/{sample}/{version}/{coverage}/genotyping.log"
	threads: 24
	resources:
		mem_total_mb=50000,
		runtime_hrs=5,
		runtime_min=1
	priority: 1
	params:
		tmpdir="genotyping/{callset}/{sample}/{version}/{coverage}/temp",
        out_prefix="genotyping/{callset}/{sample}/{version}/{coverage}/temp/{version}-{sample}",
		pangenie = lambda wildcards: config['pangenie-modules'][wildcards.version]
	wildcard_constraints:
		version = "|".join([k for k in config['pangenie-modules'].keys()] + ['^' + k for k in config['pangenie']])
	shell:
		"""
		mkdir -p {params.tmpdir}
        (/usr/bin/time -v {params.pangenie}-index -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -t {threads} ) &> {log.index}
		(/usr/bin/time -v {params.pangenie} -f {params.out_prefix} -i <(gunzip -c {input.reads}) -s {wildcards.sample} -o {params.out_prefix} -j {threads} -t {threads}  ) &> {log.genotype}
        mv {params.out_prefix}_genotyping.vcf {output.genotyping}
        rm -r {params.tmpdir}
		"""


########################################################
##################    run BayesTyper    ##################
########################################################

# run kmc to count kmers
rule run_kmc:
    input:
        reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "preprocessing/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa"
    output:
        suf=temp("genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.kmc_suf"),
        pre=temp("genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.kmc_pre")
    log:
        "logs/genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}_kmc.log"
    params:
        out_prefix="genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/tempout",
        temp_suf="genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/tempout.kmc_suf",
        temp_pre="genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/tempout.kmc_pre",
        tmp="genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/tmp/"
    threads: 24
    resources:
        mem_mb=15000,
        runtime_hrs=4,
        runtime_min=59
    shell:
        """
        mkdir -p {params.tmp}
        /usr/bin/time -v {kmc} -k55 -t{threads} -ci1 {input} {params.out_prefix} {params.tmp} > {log} 2>&1
        mv {params.temp_suf} {output.suf}
        mv {params.temp_pre} {output.pre}
        """

# create bloomfilter
rule create_bloomfilter:
    input:
        suf = "genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.kmc_suf",
        pre = "genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.kmc_pre"
    output:
        temp("genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.bloomData"),
        temp("genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.bloomMeta")
    log:
        "logs/genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}_bloom.log"
    params:
        out_prefix="genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}"
    threads: 24
    resources:
        mem_mb=20000,
        runtime_hrs=3,
        runtime_min=59
    shell:
        "/usr/bin/time -v {bayestyper_tools} makeBloom -p {threads} -k {params.out_prefix} > {log} 2>&1"

# create samples file
rule bayestyper_make_samples_file:
	output:
		temp("genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.tsv")
	params:
		prefix="genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}",
		sex=lambda wildcards: sex_per_sample[wildcards.sample]
	run:
		with open(output[0], "w") as bayestyper_samples_file:
			bayestyper_samples_file.write("{sample}\t{sex}\t{prefix}\n".format(sample=wildcards.sample, sex=params.sex, prefix=params.prefix))

# bayestyper cluster
checkpoint bayestyper_cluster:
	input:
		"data/downloaded/bayestyper_utils/bayestyper_GRCh38_bundle.tar.gz",
		"genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.bloomData",
		"genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.bloomMeta",
		vcf = "preprocessing/{callset}/{sample}/input-panel/panel_multi_norm.vcf",
		samples = "genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.tsv"
	output:
		dir=directory("genotyping/{callset}/{sample}/bayestyper/{coverage}/clusters")
	params:
		out_prefix="genotyping/{callset}/{sample}/bayestyper/{coverage}/clusters/bayestyper"
	log:
		"logs/genotyping/{callset}/{sample}/bayestyper/{coverage}/clusters/clusters-log.log"
	resources:
		mem_mb=60000,
		runtime_hrs=4,
		runtime_min=59
	threads: 24
	shell:
		"""
		mkdir -p {output.dir}
		/usr/bin/time -v {bayestyper} cluster -v {input.vcf} -s {input.samples} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} \
					-p {threads} -o {params.out_prefix} > {log} 2>&1
		"""

# bayestyper genotype
rule run_bayestyper_genotype:
	input:
		suf = "genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.kmc_suf",
        pre = "genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.kmc_pre",
		samples = "genotyping/{callset}/{sample}/bayestyper/{coverage}/temp/kmers/{sample}.tsv",
		unit = "genotyping/{callset}/{sample}/bayestyper/{coverage}/clusters/bayestyper_unit_{unit_id}/variant_clusters.bin"
	output:
		genotypes="genotyping/{callset}/{sample}/bayestyper/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
	log:
		"logs/genotyping/{callset}/{sample}/bayestyper/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.log"
	params:
		cluster_data_dir="genotyping/{callset}/{sample}/bayestyper/{coverage}/clusters/bayestyper_cluster_data",
		out_prefix="genotyping/{callset}/{sample}/bayestyper/{coverage}/genotype/bayestyper_unit_{unit_id}/tmp",
		out_zipped_vcf="genotyping/{callset}/{sample}/bayestyper/{coverage}/genotype/bayestyper_unit_{unit_id}/tmp.vcf.gz"
	threads: 24
	resources:
		mem_mb=50000,
		runtime_hrs=1,
		runtime_min=59
	shell:
		"""
        /usr/bin/time -v {bayestyper} genotype -v {input.unit} -s {input.samples} -c {params.cluster_data_dir} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} -p {threads} -z -o {params.out_prefix} > {log} 2>&1
		gunzip -c {params.out_zipped_vcf} | bgzip -c > {output.genotypes}
        tabix -p vcf {output.genotypes}
        rm {params.out_zipped_vcf}
        """

# combine vcfs
def aggregate_input_vcfs(wildcards):
	checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output[0]
	result = expand("genotyping/{callset}/{sample}/bayestyper/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
								callset=wildcards.callset,
								sample=wildcards.sample,
								coverage=wildcards.coverage,
								unit_id=glob_wildcards(os.path.join(checkpoint_output, "bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
	return sorted(result)

def aggregate_input_tbis(wildcards):
	checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output[0]
	result = expand("genotyping/{callset}/{sample}/bayestyper/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz.tbi",
								callset=wildcards.callset,
								sample=wildcards.sample,
								coverage=wildcards.coverage,
								unit_id=glob_wildcards(os.path.join(checkpoint_output, "bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
	return sorted(result)

rule bcftools_concat_units:
	input:
		vcfs=aggregate_input_vcfs,
		#tbis=aggregate_input_tbis
	output:
		"genotyping/{callset}/{sample}/bayestyper/{coverage}/genotyping.vcf"
	log:
		"logs/genotyping/{callset}/{sample}/bayestyper/{coverage}/genotyping.log"
	shell:
		"/usr/bin/time -v bcftools concat -a -o {output} {input.vcfs} &> {log}"

#########################################################
#################    run GraphTyper    ##################
#########################################################

# Split Input panel up per chromosome 
rule split_vcf_by_chromosome:
    input:
        vcf = lambda wildcards: "preprocessing/{callset}/{sample}/input-panel/panel_bi.vcf.gz",
        tbi = lambda wildcards: "preprocessing/{callset}/{sample}/input-panel/panel_bi.vcf.gz.tbi"
    output:
        "preprocessing/{callset}/{sample}/input-panel/splitted/panel_bi_chr{chrom}.vcf"
    wildcard_constraints:
        chrom="|".join(chromosomes)
    conda:
        "../envs/genotyping.yml"
    shell:
        """
        bcftools view {input.vcf} -r chr{wildcards.chrom}  > {output}
        """

# Input VCF for graphtyper must be biallelic
rule graphtyper_preprocess:
	input:
		vcf="preprocessing/{callset}/{sample}/input-panel/splitted/panel_bi_chr{chrom}.vcf.gz",
		tbi="preprocessing/{callset}/{sample}/input-panel/splitted/panel_bi_chr{chrom}.vcf.gz.tbi"
	output:
		"preprocessing/{callset}/{sample}/input-panel/panel_bi_chr{chrom}_{variant}.vcf.gz"
	wildcard_constraints:
		chrom="|".join(chromosomes),
		variant="snp-indel|sv"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view {input.vcf} | python workflow/scripts/extract-varianttype.py {wildcards.variant} | bgzip -c > {output}
        tabix -p vcf {output}
		"""

rule graphtyper_genotype:
    input:
        vcf = "preprocessing/{callset}/{sample}/input-panel/panel_bi_chr{chrom}_{variant}.vcf.gz",
        bam = "preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.chr{chrom}.bam",
        bai = "preprocessing/downsampling/{callset}/{coverage}/aligned/{sample}_full.chr{chrom}.bam.bai",
        fasta= lambda wildcards: config['callsets'][wildcards.callset]['reference']
    output:
        vcf="genotyping/{callset}/{sample}/graphtyper/{coverage}/temp/{variant}/genotyping.chr{chrom}.vcf.gz"
    params:
        dir="genotyping/{callset}/{sample}/graphtyper/{coverage}/temp/{variant}"
    wildcard_constraints:
        chrom="|".join(chromosomes),
        variant="snp-indel|sv"
    log:
        "logs/genotyping/{callset}/{sample}/graphtyper/{coverage}/temp/{variant}/chr{chrom}.log"
    conda:
        "../envs/genotyping.yml"
    threads:
        1
    resources:
        mem_total_mb=30000,
        runtime_hrs=4,
        runtime_min=59
    shell:
        """
        ulimit -n 5000
		if [ "{wildcards.variant}" == "sv" ]; then 
            (/usr/bin/time -v {graphtyper} genotype_sv {input.fasta} {input.vcf} --sam={input.bam} --region=chr{wildcards.chrom}  --output={params.dir} --threads={threads} --force_no_filter_zero_qual) &> {log}
            bcftools concat -a {params.dir}/chr{wildcards.chrom}/*.vcf.gz | python workflow/scripts/graphtyper-postprocess.py {input.vcf} | bgzip -c > {output.vcf}
            tabix -p vcf {output.vcf}
        else
            if [ "{wildcards.variant}" == "snp-indel" ]; then
                (/usr/bin/time -v {graphtyper} genotype {input.fasta} --vcf={input.vcf} --sam={input.bam} --region=chr{wildcards.chrom} --no_decompose --verbose --output={params.dir} --threads={threads} --force_no_filter_zero_qual) &> {log}
                bcftools concat -a {params.dir}/chr{wildcards.chrom}/*.vcf.gz | bgzip -c > {output.vcf}
                tabix -p vcf {output.vcf}
            fi
        fi
        """
        
		
rule merge_vcfs_all_chromosomes_and_normalize:
    input:
        #vcfs=lambda wildcards: expand("genotyping/{{callset}}/{{sample}}/graphtyper/{{coverage}}/temp/{variant}/genotyping.chr{chrom}.vcf.gz", chrom=ref_chromosomes[wildcards.callset], variant=["snp-indel", "sv"]),
        vcfs=lambda wildcards: expand("genotyping/{{callset}}/{{sample}}/graphtyper/{{coverage}}/temp/{variant}/genotyping.chr{chrom}.vcf.gz", chrom=chromosomes, variant=["snp-indel", "sv"]),
        reference=lambda wildcards: config['callsets'][wildcards.callset]['reference']
    output: 
        "genotyping/{callset}/{sample}/graphtyper/{coverage}/genotyping.vcf"
    shell: 
        """
        bcftools concat -a {input.vcfs} | bcftools norm -f {input.reference} -m -any | bcftools sort > {output}
        """



