configfile: "config.json"

pangenie = config['pangenie']
graphtyper = config['graphtyper']
bayestyper=config['bayestyper']
bayestyper_tools=config['bayestyper_tools']
kmc=config['kmc']

# Sample to be genotyped
samples = config['sample'] 
reads = config['reads']
reference = config['reference']
input_pangenome_graph = config['input_pangenome_graph']
input_callset = config['input_callset']
graphtypervariants = config['graphtyper_variants']
variants = config['variants']
tools = config['tools']
chroms = config['chromosomes']
bayestyper_reference_canon=config['bayestyper_reference_canon']
bayestyper_reference_decoy=config['bayestyper_reference_decoy']


rule all:
    input:
        expand("genotyping/{sample}/{tool}/genotyping_{variant}.vcf.gz", sample=samples, tool=tools, variant=variants)
        #"genotyping/NA24385/graphtyper/indel/genotyping.chr1.vcf.gz",
        #expand("genotyping/{sample}/bayestyper/genotyping.vcf", sample=samples),
        #
        #expand("genotyping/{sample}/{tool}/{sample}_genotyping_{variant}.vcf.gz", sample=samples, tool=tools, variant=variants),
        #expand("genotyping/{sample}/{tool}/{sample}_genotyping.vcf", sample=samples, tool=tools)
        #expand("genotyping/{sample}/{tool}/{sample}/{graphtypervariant}/genotyping.chr{chrom}.vcf.gz", sample=samples, graphtypervariant=graphtypervariants, tool=tools, chrom=chroms)


### Add uncompress all vcfs and extract sample from panel reference (in case of one-leave-out experiments)
# uncompress vcf
rule uncompress_vcfs:
	input:
		"{filename}.vcf.gz"
	output:
		"{filename}.vcf"
	conda:
		"../env/genotyping.yml"
	shell:
		"gunzip -c {input} > {output}"

# rule tabix
rule tabix:
	input:
		"{filename}.vcf.gz"
	output:
		"{filename}.vcf.gz.tbi"
	conda:
		"../env/genotyping.yml"
	shell:
		"tabix -p vcf {input}"


################################################
# genotype a sample using PanGenie
################################################
rule genotyping_pangenie:
	input:
		reads = lambda wildcards: config["data"][wildcards.sample]["reads"],
		reference = reference,
		input_pangenome_graph = input_pangenome_graph
	output:
		#reads = temp("genotyping/{sample}/pangenie/{sample}-reads.fastq"),
		#uncompressed_input_graph = temp("genotyping/{sample}/pangenie/Uncompressed_input_graph_for_{sample}.vcf"),
		genotypes = "genotyping/{sample}/pangenie/{sample}_genotyping.vcf"
	threads: 24
	log:
		"genotyping/{sample}/pangenie/{sample}.log"
	params:
		prefix = "genotyping/{sample}/pangenie/{sample}"
	resources:
		mem_total_mb=100000,
		runtime_hrs=4,
		runtime_min=59
	run:
		# jellyfish requires that the files are uncompressed
		#shell("gunzip -c {input.reads} > {output.reads}")
		#shell("gunzip -c {input.input_pangenome_graph} > {output.uncompressed_input_graph}")
		shell("(/usr/bin/time -v {pangenie} -i {input.reads} -v {input.input_pangenome_graph} -r {input.reference} -o {params.prefix} -s {wildcards.sample} -j {threads} -t {threads} -g) &> {log}")
		#shell("(/usr/bin/time -v {pangenie} -i {output.reads} -v {output.uncompressed_input_graph} -r {input.reference} -o {params.prefix} -s {wildcards.sample} -j {threads} -t {threads} -g) &> {log}")

################################################
# genotype a sample using GraphTyper
################################################

#### Note: To genotype using GraphTyper we need to consider split between
####   indels and SVs, as they have two modules (genotype and genotype_sv), 
####   each optimized for the corresponding case

rule graphtyper_preprocess:
	input:
		vcf = input_callset,
	output:
		vcf="genotyping/{sample}/graphtyper/input_variants/input_callset_{graphtypervariant}_chr{chrom}.vcf.gz"
	shell:
		"""
		cat {input.vcf} | python3 scripts/extract-varianttype.py {wildcards.graphtypervariant} | bgzip > {output.vcf}
		tabix -f -p vcf {output.vcf} 
		"""
		#tabix -f -p vcf {output.vcf} > {output.tbi}

rule graphtyper_genotype:
	input:
		vcf="genotyping/{sample}/graphtyper/input_variants/input_callset_{graphtypervariant}_chr{chrom}.vcf.gz",
		tbi="genotyping/{sample}/graphtyper/input_variants/input_callset_{graphtypervariant}_chr{chrom}.vcf.gz.tbi",
		bam="downloaded/BAM/aligned_reads.bam",
		bai="downloaded/BAM/aligned_reads.bam.bai",
		fasta=reference
	output:
		vcf="genotyping/{sample}/graphtyper/{graphtypervariant}/genotyping.chr{chrom}.vcf.gz",
		tbi="genotyping/{sample}/graphtyper/{graphtypervariant}/genotyping.chr{chrom}.vcf.gz.tbi"
	params:
		prefix = "genotyping/{sample}/graphtyper/{graphtypervariant}"
	log:
		"genotyping/{sample}/graphtyper/genotyping_{graphtypervariant}_{chrom}.log"
	threads:
		24
	resources:
		mem_total_mb=30000,
		runtime_hrs=4,
		runtime_min=59
	shell:
		"""
		if [ "{wildcards.graphtypervariant}" == "sv" ]; then 
            (/usr/bin/time -v {graphtyper} genotype_sv {input.fasta} {input.vcf} --sam={input.bam} --region=chr{wildcards.chrom}  --output={params.prefix} --threads={threads}) &> {log}
            bcftools concat -a {params.prefix}/chr{wildcards.chrom}/*.vcf.gz | python3 scripts/graphtyper-postprocess.py {input.vcf} | bgzip -c > {output.vcf}
		else
            if [ "{wildcards.graphtypervariant}" == "indel" ]; then
                (/usr/bin/time -v {graphtyper} genotype {input.fasta} --vcf={input.vcf} --sam={input.bam} --region=chr{wildcards.chrom} --no_decompose --verbose --output={params.prefix} --threads={threads}) &> {log}
                bcftools concat -a {params.prefix}/chr{wildcards.chrom}/*.vcf.gz | bgzip > {output.vcf}
            fi
        fi
        tabix -f -p vcf {output.vcf} > {output.tbi}
        """


### include graphtyper postprocessing: merging of vcf files  
## final output: genotyping/{tool}/{sample}_genotyping_{graphtypervariant}.vcf.gz

rule merge_vcfs_all_chromosomes:
    input:
       vcfs=lambda wildcards: expand("genotyping/{{sample}}/graphtyper/{graphtypervariant}/genotyping.chr{chrom}.vcf.gz", chrom=chroms, graphtypervariant=graphtypervariants)
    output: 
        "genotyping/{sample}/graphtyper/{sample}_genotyping.vcf"
    shell: 
        """
        bcftools concat -a {input.vcfs} > {output}
        """

##############################################################################
##################    run BayesTyper genotyping pipeline    ##################
##############################################################################


# run kmc to count kmers
rule run_kmc:
    input:
        reads= lambda wildcards: config["data"][wildcards.sample]["reads"]
    output:
        #uncompressed_reads="genotyping/{sample}/bayestyper/reads.fastq",
        suf="genotyping/{sample}/bayestyper/kmers/sample.kmc_suf",
        pre="genotyping/{sample}/bayestyper/kmers/sample.kmc_pre"
    log:
        "genotyping/{sample}/bayestyper/kmers/sample_kmc.log"
    params:
        out_prefix="genotyping/{sample}/bayestyper/kmers/tempout",
        temp_suf="genotyping/{sample}/bayestyper/kmers/tempout.kmc_suf",
        temp_pre="genotyping/{sample}/bayestyper/kmers/tempout.kmc_pre",
        tmp="genotyping/{sample}/bayestyper/kmers/tmp/"
    threads: 24
    resources:
        mem_total_mb=15000,
        runtime_hrs=4,
        runtime_min=59
    shell:
        """
        mkdir -p {params.tmp}
        /usr/bin/time -v {kmc} -k55 -t{threads} -ci1 {input.reads} {params.out_prefix} {params.tmp} > {log} 2>&1
        mv {params.temp_suf} {output.suf}
        mv {params.temp_pre} {output.pre}
        """
        #/usr/bin/time -v {kmc} -k55 -t{threads} -ci1 {output.uncompressed_reads} {params.out_prefix} {params.tmp} > {log} 2>&1
        #gunzip -c {input.reads} > {output.uncompressed_reads}

# create bloomfilter
rule create_bloomfilter:
    input:
        "genotyping/{sample}/bayestyper/kmers/sample.kmc_pre",
        "genotyping/{sample}/bayestyper/kmers/sample.kmc_suf"
    output:
        "genotyping/{sample}/bayestyper/kmers/sample.bloomData",
        "genotyping/{sample}/bayestyper/kmers/sample.bloomMeta"
    log:
        "genotyping/{sample}/bayestyper/kmers/sample.log"
    params:
        out_prefix="genotyping/{sample}/bayestyper/kmers/sample"
    threads: 24
    resources:
        mem_total_mb=20000,
        runtime_hrs=3,
        runtime_min=59
    shell:
        "/usr/bin/time -v {bayestyper_tools} makeBloom -p {threads} -k {params.out_prefix} > {log} 2>&1"

# create samples file
rule bayestyper_make_samples_file:
    output:
        "genotyping/{sample}/bayestyper/sample.tsv"
    params:
        prefix="genotyping/{sample}/bayestyper/kmers/sample",
        sex=lambda wildcards: config['data'][wildcards.sample]['sex']
    run:
        with open(output[0], "w") as bayestyper_samples_file:
            bayestyper_samples_file.write("{sample}\t{sex}\t{prefix}\n".format(sample=wildcards.sample, sex=params.sex, prefix=params.prefix))

# bayestyper cluster
checkpoint bayestyper_cluster:
    input:
        "genotyping/{sample}/bayestyper/kmers/sample.bloomData",
        "genotyping/{sample}/bayestyper/kmers/sample.bloomMeta",
        variants=input_pangenome_graph,
        samples="genotyping/{sample}/bayestyper/sample.tsv"
    output:
        dir=directory("genotyping/{sample}/bayestyper/clusters/")
    params:
        out_prefix="genotyping/{sample}/bayestyper/clusters/bayestyper"
    log:
        "genotyping/{sample}/bayestyper/clusters/clusters-log.log"
    resources:
        mem_total_mb=60000,
        runtime_hrs=4,
        runtime_min=59
    threads: 24
    shell:
        """
        /usr/bin/time -v {bayestyper} cluster -v {input.variants} -s {input.samples} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} \
                    -p {threads} -o {params.out_prefix} > {log} 2>&1
        """

# bayestyper genotype
rule run_bayestyper_genotype:
	input:
		kmc_pre="genotyping/{sample}/bayestyper/kmers/sample.kmc_pre",
		kmc_suf="genotyping/{sample}/bayestyper/kmers/sample.kmc_suf",
		samples="genotyping/{sample}/bayestyper/sample.tsv",
		unit="genotyping/{sample}/bayestyper/clusters/bayestyper_unit_{unit_id}/variant_clusters.bin"
	output:
		genotypes="genotyping/{sample}/bayestyper/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
		kmer_coverage_file="genotyping/{sample}/bayestyper/genotype/bayestyper_unit_{unit_id}/bayestyper_genomic_parameters.txt"
	log:
		"genotyping/{sample}/bayestyper/genotype/bayestyper_unit_{unit_id}/bayestyper.log"
	params:
		cluster_data_dir="genotyping/{sample}/bayestyper/clusters/bayestyper_cluster_data",
		out_prefix="genotyping/{sample}/bayestyper/genotype/bayestyper_unit_{unit_id}/bayestyper"
	threads: 24
	resources:
		mem_total_mb=50000,
		runtime_hrs=1,
		runtime_min=59
	run:
		shell("/usr/bin/time -v {bayestyper} genotype -v {input.unit} -s {input.samples} -c {params.cluster_data_dir} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} \
		-p {threads} -z -o {params.out_prefix} > {log} 2>&1")
		# fix the vcf ...
		shell("gunzip -c {output.genotypes} | bgzip > {output.genotypes}-tmp")
		shell("mv {output.genotypes}-tmp  {output.genotypes}")
		shell("tabix -p vcf {output.genotypes}")


# combine vcfs
def aggregate_input(wildcards):
	checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output.dir
	result = expand("genotyping/{sample}/bayestyper/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
                                sample=wildcards.sample,
								unit_id=glob_wildcards(os.path.join(checkpoint_output, "bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
	return sorted(result)

rule bcftools_concat_units:
    input:
        aggregate_input
    output:
        "genotyping/{sample}/bayestyper/{sample}_genotyping.vcf"
    log:
        "genotyping/{sample}/bayestyper/genotyping.log"
    shell:
        "/usr/bin/time -v bcftools concat -a -o {output} {input} &> {log}"


################################################
# split variants and compress genotyped files
################################################

### NOTE: Be sure of convert genotyped files to biallelic before extracting variants
### GraphTyper's output is already biallelic, but not for BayesTyper and PanGenie 

rule split_variants_and_compress_genotyped:
	input:
		vcf="genotyping/{sample}/{tool}/{sample}_genotyping.vcf"
	output:
		splitted_vcf="genotyping/{sample}/{tool}/genotyping_{variant}.vcf.gz" 
	shell:
		"""
		cat {input.vcf} | python3 scripts/extract-varianttype.py {wildcards.variant} | bgzip > {output.splitted_vcf}
		tabix -p vcf {output}
		"""