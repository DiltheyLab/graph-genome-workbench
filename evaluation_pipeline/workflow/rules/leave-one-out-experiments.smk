kmc=config['programs']['kmc']
bayestyper=config['programs']['bayestyper']
bayestyper_tools=config['programs']['bayestyper_tools']
graphtyper=config['programs']['graphtyper']
rtg=config['programs']['rtg']
truvari = config['programs']['truvari']

# parameters
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
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

allowed_variants = ['snp', 'indels', 'large-deletion', 'large-insertion']
callsets_leave_one_out = [s for s in config['callsets'].keys()]
coverages_leave_one_out = ['full'] + [c for c in config['downsampling']]
versions_leave_one_out = [v for v  in config['pangenie'].keys()] + [v for v in config['pangenie-modules'].keys()] + ['bayestyper', 'graphtyper'] # , 'graphtyper'


################################################################
######   prepare input panel and ground truth genotypes  #######
################################################################

# remove positions that are ".|." in the left out sample. These cannot be used for evaluation, as the true genotype is unknown
# for now, also remove CHM13, because PanGenie cannot handle haploid samples
# remove sample from the panel
rule remove_missing:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['multi'] if wildcards.representation == 'multi' else config['callsets'][wildcards.callset]['bi']
	output:
		temp("results/leave-one-out/{callset}/preprocessed-vcfs/{sample}-{callset}_{representation}_no-missing.vcf")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000
	priority: 1
	wildcard_constraints:
		representation = "bi|multi"
	shell:
		"zcat {input} | python3 workflow/scripts/remove-missing.py {wildcards.sample} > {output}"

### We need to annotated the biallelic input panel to not forget ID variants in graphtyper postprocess
rule prepare_panel:
	input:
		"results/leave-one-out/{callset}/preprocessed-vcfs/{sample}-{callset}_{representation}_no-missing.vcf.gz"
	output:
		"results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_{representation}.vcf"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	wildcard_constraints:
		representation = "bi|multi"
	params: 
		lambda wildcards: "| python3 workflow/scripts/annotate-ids.py" if wildcards.representation == 'bi' else ""
	log:
		"results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_{representation}.log"
	resources:
		mem_total_mb=20000
	shell:
		"""
		(/usr/bin/time -v bcftools view --samples ^{wildcards.sample} {input} | bcftools view --min-ac 1 {params} > {output} ) &> {log}
		"""


rule normalize_input_panel:
	input:
		vcf="results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_multi.vcf.gz",
		reference=lambda wildcards: config['callsets'][wildcards.callset]['reference']
	output:
		"results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_multi_norm.vcf"
	shell:
		"bcftools norm -m +any -d all -f {input.reference} {input.vcf} | bcftools sort > {output}"

# extract ground truth genotypes for sample 
rule prepare_truth:
	input:
		"results/leave-one-out/{callset}/preprocessed-vcfs/{sample}-{callset}_bi_no-missing.vcf.gz"
	output:
		"results/leave-one-out/{callset}/truth/truth-{sample}-{callset}.vcf"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	resources:
		mem_total_mb=20000
	log:
		"results/leave-one-out/{callset}/truth/truth-{sample}-{callset}.log"
	shell:
		"bcftools view --samples {wildcards.sample} {input} 2> {log} 1> {output}"


rule compress_vcf:
	input:
		temp("results/leave-one-out/{filename}.vcf")
	output:
		vcf = "results/leave-one-out/{filename}.vcf.gz",
		tbi = "results/leave-one-out/{filename}.vcf.gz.tbi"
	priority: 1
	shell:
		"""
		bgzip -c {input} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


########################################################
##################    run PanGenie    ##################
########################################################


# run pangenie
## If in input reads = lambda wildcards: "results/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
## the reads will be aligned before, but that's what we need, since graphtyper is a mapped-based approach (i.e. we need the aligned reads anyways if we want to apply Graphtyper). Just for Pangenie, it's not necessary.
rule pangenie:
	input:
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "results/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		vcf="results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_multi_norm.vcf"
	output:
		genotyping = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/{version}-{sample}_genotyping.vcf")
	log:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}.log"
	threads: 24
	resources:
		mem_total_mb=80000,
		runtime_hrs=5,
		runtime_min=1
	priority: 1
	params:
		out_prefix="results/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/{version}-{sample}",
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
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "results/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		vcf="results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_multi_norm.vcf"
	output:
		genotyping = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/{version}-{sample}_genotyping.vcf")
	log:
		index = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_index.log",
		genotype = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping.log"
	threads: 24
	resources:
		mem_total_mb=50000,
		runtime_hrs=5,
		runtime_min=1
	priority: 1
	params:
		out_prefix="results/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/{version}-{sample}",
		pangenie = lambda wildcards: config['pangenie-modules'][wildcards.version]
	wildcard_constraints:
		version = "|".join([k for k in config['pangenie-modules'].keys()] + ['^' + k for k in config['pangenie']])
	shell:
		"""
		(/usr/bin/time -v {params.pangenie}-index -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -t {threads} ) &> {log.index}
		(/usr/bin/time -v {params.pangenie} -f {params.out_prefix} -i <(gunzip -c {input.reads}) -s {wildcards.sample} -o {params.out_prefix} -j {threads} -t {threads}  ) &> {log.genotype}
		"""


########################################################
##################    run BayesTyper    ##################
########################################################

# run kmc to count kmers
rule run_kmc:
    input:
        reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "results/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa"
    output:
        suf=temp("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.kmc_suf"),
        pre=temp("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.kmc_pre")
    log:
        "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}_kmc.log"
    params:
        out_prefix="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/tempout",
        temp_suf="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/tempout.kmc_suf",
        temp_pre="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/tempout.kmc_pre",
        tmp="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/tmp/"
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
        suf = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.kmc_suf",
        pre = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.kmc_pre"
    output:
        temp("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.bloomData"),
        temp("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.bloomMeta")
    log:
        "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}_bloom.log"
    params:
        out_prefix="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}"
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
		temp("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.tsv")
	params:
		prefix="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}",
		sex=lambda wildcards: sex_per_sample[wildcards.sample]
	run:
		with open(output[0], "w") as bayestyper_samples_file:
			bayestyper_samples_file.write("{sample}\t{sex}\t{prefix}\n".format(sample=wildcards.sample, sex=params.sex, prefix=params.prefix))

# bayestyper cluster
checkpoint bayestyper_cluster:
	input:
		"data/downloaded/bayestyper_utils/bayestyper_GRCh38_bundle.tar.gz",
		"results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.bloomData",
		"results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.bloomMeta",
		vcf="results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_multi_norm.vcf",
		samples="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.tsv"
	output:
		dir=directory("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/clusters")
	params:
		out_prefix="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/clusters/bayestyper"
	log:
		"results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/clusters/clusters-log.log"
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
		suf = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.kmc_suf",
        pre = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.kmc_pre",
		samples = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}.tsv",
		unit = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/clusters/bayestyper_unit_{unit_id}/variant_clusters.bin"
	output:
		genotypes=temp("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf"),
	log:
		"results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.log"
	params:
		cluster_data_dir="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/clusters/bayestyper_cluster_data",
		out_prefix="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper_zipped",
		out_zipped_vcf="results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper_zipped.vcf.gz"
	threads: 24
	resources:
		mem_mb=50000,
		runtime_hrs=1,
		runtime_min=59
	run:
		shell("/usr/bin/time -v {bayestyper} genotype -v {input.unit} -s {input.samples} -c {params.cluster_data_dir} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} \
		-p {threads} -z -o {params.out_prefix} > {log} 2>&1")
		# fix the vcf ...
		shell("gunzip -c {params.out_zipped_vcf} > {output.genotypes}")

# combine vcfs
def aggregate_input_vcfs(wildcards):
	checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output[0]
	result = expand("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
								callset=wildcards.callset,
								sample=wildcards.sample,
								coverage=wildcards.coverage,
								unit_id=glob_wildcards(os.path.join(checkpoint_output, "bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
	return sorted(result)

def aggregate_input_tbis(wildcards):
	checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output[0]
	result = expand("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz.tbi",
								callset=wildcards.callset,
								sample=wildcards.sample,
								coverage=wildcards.coverage,
								unit_id=glob_wildcards(os.path.join(checkpoint_output, "bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
	return sorted(result)

rule bcftools_concat_units:
	input:
		vcfs=aggregate_input_vcfs,
		tbis=aggregate_input_tbis
	output:
		"results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/bayestyper-{sample}_genotyping.vcf"
	log:
		"results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/bayestyper-{sample}_concat-units.log"
	shell:
		"/usr/bin/time -v bcftools concat -a -o {output} {input.vcfs} &> {log}"

#########################################################
#################    run GraphTyper    ##################
#########################################################

# Split Input panel up per chromosome 
rule split_vcf_by_chromosome:
    input:
        vcf=lambda wildcards: "results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_{representation}.vcf.gz",
        tbi=lambda wildcards: "results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_{representation}.vcf.gz.tbi"
    output:
        temp("results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_{representation}_chr{chrom}.vcf")
    wildcard_constraints:
        chrom="X|Y|[0-9]+",
        representation="multi|bi"
    conda:
        "../envs/genotyping.yml"
    shell:
        """
        bcftools view {input.vcf} -r chr{wildcards.chrom} > {output}
        """

# Input VCF for graphtyper must be biallelic
rule graphtyper_preprocess:
	input:
		vcf="results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_bi_chr{chrom}.vcf.gz",
		tbi="results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_bi_chr{chrom}.vcf.gz.tbi"
	output:
		vcf=temp("results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_bi_chr{chrom}_{variant}.vcf")
	wildcard_constraints:
		chrom="X|Y|[0-9]+",
		variant="indel|sv"
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		zcat {input.vcf} | python3 workflow/scripts/extract-varianttype.py {wildcards.variant} > {output}
		"""

## Remark: indel contains all variants that are not SVs, including SNPs, as well.
rule graphtyper_genotype:
    input:
        vcf = "results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_bi_chr{chrom}_{variant}.vcf.gz",
        tbi = "results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_bi_chr{chrom}_{variant}.vcf.gz.tbi",
        bam = "results/downsampling/{callset}/{coverage}/aligned/{sample}_full.chr{chrom}.bam",
        bai = "results/downsampling/{callset}/{coverage}/aligned/{sample}_full.chr{chrom}.bam.bai",
        fasta= lambda wildcards: config['callsets'][wildcards.callset]['reference']
    output:
        vcf=temp("results/leave-one-out/{callset}/graphtyper/{sample}/{coverage}/temp/{variant}/graphtyper-{sample}_genotyping.chr{chrom}.vcf")
    params:
        dir="results/leave-one-out/{callset}/graphtyper/{sample}/{coverage}/temp/{variant}"
    wildcard_constraints:
        chrom="|".join(chromosomes),
        variant="indel|sv"
        #chrom="(X|Y|[0-9]+)",
    log:
        "results/leave-one-out/{callset}/graphtyper/{sample}/{coverage}/temp/{variant}/chr{chrom}.log"
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
            bcftools concat -a {params.dir}/chr{wildcards.chrom}/*.vcf.gz | python3 workflow/scripts/graphtyper-postprocess.py {input.vcf} > {output.vcf}
        else
            if [ "{wildcards.variant}" == "indel" ]; then
                (/usr/bin/time -v {graphtyper} genotype {input.fasta} --vcf={input.vcf} --sam={input.bam} --region=chr{wildcards.chrom} --no_decompose --verbose --output={params.dir} --threads={threads} --force_no_filter_zero_qual) &> {log}
                bcftools concat -a {params.dir}/chr{wildcards.chrom}/*.vcf.gz > {output.vcf}
            fi
        fi
        """
        
		
rule merge_vcfs_all_chromosomes_and_normalize:
    input:
        vcfs=lambda wildcards: expand("results/leave-one-out/{{callset}}/graphtyper/{{sample}}/{{coverage}}/temp/{variant}/graphtyper-{{sample}}_genotyping.chr{chrom}.vcf.gz", chrom=chromosomes, variant=["indel", "sv"]),
        reference=lambda wildcards: config['callsets'][wildcards.callset]['reference']
    output: 
        "results/leave-one-out/{callset}/graphtyper/{sample}/{coverage}/temp/graphtyper-{sample}_genotyping.vcf"
    shell: 
        """
        bcftools concat -a {input.vcfs} | bcftools norm -f {input.reference} -m -any | bcftools sort > {output}
        """


########################################################
##################    Evaluation      ##################
########################################################


rule alleles_per_bubble:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['multi']
	output:
		plot = "results/leave-one-out/{callset}/alleles-per-bubble.pdf",
		bed = "results/leave-one-out/{callset}/complex-bubbles.bed"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000,
		runtime_hrs=1
	shell:
		"zcat {input} | python3 workflow/scripts/variant-statistics.py {output.plot} 1 > {output.bed}"


# prepare beds for biallelic and complex graph regions
rule prepare_beds:
	input:
		bed = "results/leave-one-out/{callset}/complex-bubbles.bed",
		fai = lambda wildcards: config['callsets'][wildcards.callset]['reference'] + '.fai'
	output:
		bed = "results/leave-one-out/{callset}/biallelic-bubbles.bed",
		tmp = temp("results/leave-one-out/{callset}/biallelic-bubbles.fai")
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		sort -k1,1d -k 2,2n -k 3,3n {input.fai} > {output.tmp}
		bedtools complement -i {input.bed} -g {output.tmp} > {output.bed}
		"""

def get_panel_vcf(wildcards):
	if wildcards.version == 'graphtyper':
		return "results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_bi.vcf"
	else:
		return "results/leave-one-out/{callset}/input-panel/panel-{sample}-{callset}_multi_norm.vcf"

# convert genotyped VCF to biallelic representation
rule convert_genotypes_to_biallelic:
	input:
		genotyped_vcf = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/temp/{version}-{sample}_genotyping.vcf.gz",
        panel_vcf = get_panel_vcf,
		biallelic = lambda wildcards: config['callsets'][wildcards.callset]['bi']
	output:
		temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping-biallelic.vcf")
	log:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping-biallelic.log"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=30000
	priority: 1
	shell:
		"(zcat {input.genotyped_vcf} | python3 workflow/scripts/annotate.py {input.panel_vcf} | python3 workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}) &> {log}"


# determine untypable ids
rule untypable_ids:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['bi']
	output:
		lists = "results/leave-one-out/{callset}/untypable-ids/{sample}-untypable.tsv",
		summary = "results/leave-one-out/{callset}/untypable-ids-{sample}.tsv"
	params:
		out = "results/leave-one-out/{callset}/untypable-ids"
	conda:
		"../envs/genotyping.yml"
	shell:
		"zcat {input} | python3 workflow/scripts/untypable-ids-single.py {params.out} {wildcards.sample} > {output.summary}"


# determine untypable IDs
rule remove_untypable:
	input:
		vcf = "results/leave-one-out/{callset}/{path}{sample}{other}.vcf",
		ids = "results/leave-one-out/{callset}/untypable-ids/{sample}-untypable.tsv"
	output:
		vcf = "results/leave-one-out/{callset}/{path}{sample}{other}-typable-{vartype}.vcf.gz",
		tbi = "results/leave-one-out/{callset}/{path}{sample}{other}-typable-{vartype}.vcf.gz.tbi"
	wildcard_constraints:
		callset = "|".join([v for v in config['callsets'].keys()]),
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		vartype = "|".join(allowed_variants)
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 1
	priority: 1
	shell:
		"""
		cat {input.vcf} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

rule rtg_format:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['reference']
	output:
		directory("results/leave-one-out/{callset}/SDF")
	resources:
		mem_total_mb=20000
	priority: 1
	shell:
		'{rtg} format -o {output} {input}'


def region_to_bed(wildcards):
	if wildcards.regions == "biallelic":
		return "results/leave-one-out/{callset}/biallelic-bubbles.bed".format(callset=wildcards.callset)
	if wildcards.regions == "multiallelic":
		return "results/leave-one-out/{callset}/complex-bubbles.bed".format(callset=wildcards.callset)
	assert(False)


# precision-recall
rule vcfeval:
	input:
		callset = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline = "results/leave-one-out/{callset}/truth/truth-{sample}-{callset}-typable-{vartype}.vcf.gz",
		baseline_tbi = "results/leave-one-out/{callset}/truth/truth-{sample}-{callset}-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		sdf = "results/leave-one-out/{callset}/SDF"
	output:
		summary = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "biallelic|multiallelic",
		vartype = "snp|indels"
	params:
		tmp = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}_tmp",
		outname = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}",
		which = "--all-records"
	log:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/logs/{regions}_{vartype}.log"
	threads: 27
	resources:
		mem_total_mb = 120000,
		runtime_hrs = 1,
		runtime_min = 40
	shell:
		"""
		{rtg} vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 --threads {threads} &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""

# precision-recall
rule truvari:
	input:
		callset = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline = "results/leave-one-out/{callset}/truth/truth-{sample}-{callset}-typable-{vartype}.vcf.gz",
		baseline_tbi = "results/leave-one-out/{callset}/truth/truth-{sample}-{callset}-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		reference = lambda wildcards: config['callsets'][wildcards.callset]['reference']
	output:
		summary_json = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}/summary.json",
		summary_txt = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "biallelic|multiallelic",
		vartype = "large-deletion|large-insertion"
	params:
		tmp = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}_tmp",
		outname = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/{regions}_{vartype}",
	log:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/precision-recall-typable/logs/{regions}_{vartype}.log"
	threads: 1
	resources:
		mem_total_mb = 40000,
		runtime_hrs = 1,
		runtime_min = 40
	shell:
		"""
		{truvari} bench -b {input.baseline} -c {input.callset} -f {input.reference} -o {params.tmp} --pick multi --includebed {input.regions} -r 2000 --no-ref a -C 2000 --passonly &> {log} 
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		python3 workflow/scripts/parse_json_to_txt.py {output.summary_json} > {output.summary_txt}
		"""


# determine the variants that went into re-typing per category
rule collect_typed_variants:
	input:
		callset = "results/leave-one-out/{callset}/preprocessed-vcfs/{sample}-{callset}_bi_no-missing.vcf.gz",
		regions= region_to_bed,
		ids="results/leave-one-out/{callset}/untypable-ids/{sample}-untypable.tsv"
	output:
		"results/leave-one-out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "biallelic|multiallelic",
		vartype = "|".join(allowed_variants)
	resources:
		mem_total_mb=50000
	priority: 1
	shell:
		"zcat {input.callset} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 workflow/scripts/get_ids.py > {output}"


# compute concordances
rule genotype_concordances:
	input:
		callset = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline = "results/leave-one-out/{callset}/truth/truth-{sample}-{callset}-typable-{vartype}.vcf.gz",
		baseline_tbi = "results/leave-one-out/{callset}/truth/truth-{sample}-{callset}-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		typed_ids = "results/leave-one-out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	output:
		tmp_vcf1 = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}_base.vcf"),
		tmp_vcf2 = temp("results/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}_call.vcf"),
		summary = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "biallelic|multiallelic",
		vartype = "|".join(allowed_variants)
	log:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/concordance/{regions}_{vartype}/summary.log"
	resources:
		mem_total_mb = 40000,
		runtime_hrs = 0,
		runtime_min = 40
	priority: 1
	shell:
		"""
		bedtools intersect -header -a {input.baseline} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf1}
		bedtools intersect -header -a {input.callset} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf2}
		python3 workflow/scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual 0 2> {log} 1> {output.summary}
		"""


########################################################
##################     Plotting      ###################
########################################################


# collect results across all samples
rule collect_results_vcfeval:
	input:
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{{version}}/{sample}/{{coverage}}/{{metric}}/{{regions}}_{{vartype}}/summary.txt", sample = config['callsets'][wildcards.callset]['leave_one_out_samples'])
	output:
		"results/leave-one-out/{callset}/{version}/plots/{coverage}/{metric}-{callset}-{version}-{coverage}_{regions}_{vartype}.tsv"
	params:
		samples = lambda wildcards: ','.join([c for c in config['callsets'][wildcards.callset]['leave_one_out_samples']]),
		outfile = "results/leave-one-out/{callset}/{version}/plots/{coverage}/{metric}-{callset}-{version}-{coverage}_{regions}_{vartype}",
		folder = "results/leave-one-out/{callset}/{version}"
	priority: 1
	shell:
		"python3 workflow/scripts/collect-results.py {wildcards.metric} {wildcards.coverage} {params.samples} {wildcards.regions} -variants {wildcards.vartype} -folder {params.folder} -outfile {params.outfile}"



# plot results of different subsampling runs
rule plotting_versions:
	input:
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{version}/plots/{{coverage}}/{m}-{{callset}}-{version}-{{coverage}}_{{regions}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', version=versions_leave_one_out, vartype=config['callsets'][wildcards.callset]['variants'])
	output:
		"results/leave-one-out/{callset}/plots/comparison-versions/{metric}/{metric}_{coverage}_{regions}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped",
		regions="biallelic|multiallelic"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + v + '-' + wildcards.coverage + '_' + wildcards.regions for v in versions_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"



# plot results of different subsampling runs, comparing concordance and typed variants per sample
rule plotting_versions_conc_vs_untyped:
	input:
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{version}/plots/{{coverage}}/concordance-{{callset}}-{version}-{{coverage}}_{{regions}}_{vartype}.tsv", version=versions_leave_one_out, vartype=config['callsets'][wildcards.callset]['variants'])
	output:
		"results/leave-one-out/{callset}/plots/comparison-versions/concordance-vs-untyped/concordance-vs-untyped_{coverage}_{regions}.pdf"
	wildcard_constraints:
		regions="biallelic|multiallelic"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + v + '-' + wildcards.coverage + '_' + wildcards.regions for v in versions_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric concordance-vs-untyped"
		


# plot results of different coverages
rule plotting_coverages:
	input:
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{{version}}/plots/{coverage}/{m}-{{callset}}-{{version}}-{coverage}_{{regions}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', coverage=coverages_leave_one_out, vartype=config['callsets'][wildcards.callset]['variants'])
	output:
		"results/leave-one-out/{callset}/plots/comparison-coverages/{metric}/{metric}_{version}_{regions}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped",
		regions="biallelic|multiallelic"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + wildcards.version + '-' + c + '_' + wildcards.regions for c in coverages_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"

def collect_logs_all_unit_ids_bayestyper(wildcards):
	checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output[0]
	result = expand("results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.log",
					callset=wildcards.callset,
					sample=wildcards.sample,
					coverage=wildcards.coverage,
					unit_id=glob_wildcards(os.path.join(checkpoint_output, "bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
	return sorted(result)

rule collect_runtime_info_bayestyper:
	input:
		bayestyper_kmers = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}_kmc.log",
		bayestyper_bloomfilter = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}_bloom.log",
		bayestyper_clusters = "results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/clusters/clusters-log.log",
		bayestyper_genotype = collect_logs_all_unit_ids_bayestyper
	output:
		"results/leave-one-out/{callset}/bayestyper/{sample}/{coverage}/bayestyper-{sample}.log"
	shell:
		"cat {input.bayestyper_kmers} {input.bayestyper_bloomfilter} {input.bayestyper_clusters} {input.bayestyper_genotype} > {output}"

rule collect_runtime_info_graphtyper:
	input:
		graphtyper_aligning_reads = "results/downsampling/{callset}/{coverage}/aligned/{sample}_full_mem.log",
		graphtyper_indexing_reads = lambda wildcards: expand("results/downsampling/{{callset}}/{{coverage}}/aligned/{{sample}}_full.chr{chrom}-index.log", chrom=chromosomes),
		graphtyper_genotype = lambda wildcards: expand("results/leave-one-out/{{callset}}/graphtyper/{{sample}}/{{coverage}}/temp/{variant}/chr{chrom}.log", chrom=chromosomes, variant=["indel", "sv"])
	output:
		"results/leave-one-out/{callset}/graphtyper/{sample}/{coverage}/graphtyper-{sample}.log"
	conda:
		"../envs/genotyping.yml"
	shell:
		"cat {input.graphtyper_aligning_reads} {input.graphtyper_indexing_reads} {input.graphtyper_genotype} > {output}"

rule collect_runtime_info_pangenie:
	input:
		pangenie_index = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_index.log",
		pangenie_genotype = "results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}_genotyping.log"
	output:
		"results/leave-one-out/{callset}/{version}/{sample}/{coverage}/{version}-{sample}.log"
	wildcard_constraints:
		version = "|".join([k for k in config['pangenie-modules'].keys()] + ['^' + k for k in config['pangenie']])
	shell:
		"cat {input.pangenie_index} {input.pangenie_genotype} > {output}"
        
        

# plot resources (single core CPU time and max RSS) for different subsampling runs
rule plotting_resources:
	input:
		lambda wildcards: expand("results/leave-one-out/{{callset}}/{version}/{sample}/{{coverage}}/{version}-{sample}.log", version = versions_leave_one_out, sample = config['callsets'][wildcards.callset]['leave_one_out_samples'])
	output:
		"results/leave-one-out/{callset}/plots/resources/resources_{callset}-{coverage}.pdf"
	conda:
		"../envs/genotyping.yml"
	params:
		outname = "results/leave-one-out/{callset}/plots/resources/resources_{callset}-{coverage}",
		samples	= lambda wildcards: " ".join(config['callsets'][wildcards.callset]['leave_one_out_samples']),
		versions = " ".join(versions_leave_one_out)
	shell:
		"python3 workflow/scripts/plot-resources.py -files {input} -outname {params.outname} -samples {params.samples} -sizes {params.versions}"
