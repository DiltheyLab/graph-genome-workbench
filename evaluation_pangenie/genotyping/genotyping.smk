configfile:'config.json'
include: "prepare-reads.smk"

### Parameters ###

# input data for genotyping
samples=config['data']['samples']
input_reference=config['data']['reference']
reference_version=config['data']['reference_version']
repeats_bed=config['data']['repeats']
complex_bed=config['data']['complex']

# programs
pangenie=config['programs']['pangenie']
kmc=config['programs']['kmc']
bayestyper=config['programs']['bayestyper']
bayestyper_tools=config['programs']['bayestyper_tools']
graphtyper=config['programs']['graphtyper']
rtg=config['programs']['rtg']

# parameters
chromosomes=config['parameters']['chromosomes']
outname=config['parameters']['outname']
outname_reads=config['parameters']['outname_reads']
downsampling=['full']
bayestyper_reference_canon=config['parameters']['bayestyper_reference_canon']
bayestyper_reference_decoy=config['parameters']['bayestyper_reference_decoy']
other_methods=['bayestyper', 'graphtyper']
variants=['snp', 'small-deletion', 'small-insertion', 'midsize-deletion', 'midsize-insertion', 'large-deletion', 'large-insertion', 'indel', 'sv', 'small', 'midsize', 'large']

metric_to_script = {
	'precision-recall-all' : '../scripts/plot-precision-recall.py',
	'precision-recall-typable' : '../scripts/plot-precision-recall.py',
	'concordance' : '../scripts/plot-concordances.py',
	'fscore' : '../scripts/plot-fscores.py'
}

##############################################################
##################    prepare input VCF    ##################
##############################################################

# uncompress vcf
rule extract_vcf:
	input:
		lambda wildcards: config['data'][wildcards.sample]['graph'] if wildcards.mode == "reference-panel" else config['data'][wildcards.sample]['biallelic'] + ".gz"
		#lambda wildcards: config['data'][wildcards.sample]['graph'] + '.gz' if wildcards.mode == "reference-panel" else config['data'][wildcards.sample]['biallelic'] + ".gz"
	output:
		"{results}/{sample}/{mode, biallelic|reference-panel}/{sample}-all.vcf"
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


# split VCF by chromosome
rule split_vcf_by_chromosome:
    input:
        vcfs=lambda wildcards: config['data'][wildcards.sample]['graph'] if wildcards.mode == "reference-panel" else (config['data'][wildcards.sample]['biallelic'] + '.gz'), 
        tbi=lambda wildcards: config['data'][wildcards.sample]['graph'] + '.tbi' if wildcards.mode == "reference-panel" else (config['data'][wildcards.sample]['biallelic'] + '.gz.tbi')
    output:
        vcf="{results}/{sample}/{mode}/{sample}-chr{chrom}.vcf",
        gz="{results}/{sample}/{mode}/{sample}-chr{chrom}.vcf.gz"
    wildcard_constraints:
        chrom="X|Y|[0-9]+",
        mode="reference-panel|biallelic"
    params:
        prefix='' if ('37' in reference_version) or ('19' in reference_version) else 'chr'
    conda:
        "../env/genotyping.yml"
    shell:
        """
        bcftools view {input.vcfs} -r {params.prefix}{wildcards.chrom} > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        #		tabix -p vcf {output.gz}
        """


##############################################################################
##################    run BayesTyper genotyping pipeline    ##################
##############################################################################


# run kmc to count kmers
rule run_kmc:
    input:
        lambda wildcards: config['data'][wildcards.sample]['reads'],
        #outname_reads + "/{sample}/raw/{sample}-{fraction, [0-9.]+}.fastq",
    output:
        suf="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_suf",
        pre="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_pre",
    log:
        "{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}_kmc.log"
    params:
        out_prefix="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/tempout",
        temp_suf="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/tempout.kmc_suf",
        temp_pre="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/tempout.kmc_pre",
        tmp="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/tmp/"
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
        "{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_pre",
        "{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_suf"
    output:
        "{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.bloomData",
        "{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.bloomMeta"
    log:
        "{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}_bloom.log"
    params:
        out_prefix="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}"
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
		"{results}/{sample}/bayestyper/{sample}_{fraction}/{sample}.tsv"
	params:
		prefix="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}",
		sex=lambda wildcards: config['data'][wildcards.sample]['sex']
	run:
		with open(output[0], "w") as bayestyper_samples_file:
			bayestyper_samples_file.write("{sample}\t{sex}\t{prefix}\n".format(sample=wildcards.sample, sex=params.sex, prefix=params.prefix))

# bayestyper cluster
checkpoint bayestyper_cluster:
	input:
		"{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.bloomData",
		"{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.bloomMeta",
		variants="{results}/{sample}/reference-panel/{sample}-all.vcf",
		#variants=lambda wildcards: config['data'][wildcards.sample]['graph'],
		samples="{results}/{sample}/bayestyper/{sample}_{fraction}/{sample}.tsv"
	output:
		dir=directory("{results}/{sample}/bayestyper/{sample}_{fraction}/clusters/")
	params:
		out_prefix="{results}/{sample}/bayestyper/{sample}_{fraction}/clusters/bayestyper"
	log:
		"{results}/{sample}/bayestyper/{sample}_{fraction}/clusters/clusters-log.log"
	resources:
		mem_mb=60000,
		runtime_hrs=4,
		runtime_min=59
	threads: 24
	shell:
		"""
		mkdir -p {output.dir}
		/usr/bin/time -v {bayestyper} cluster -v {input.variants} -s {input.samples} -g {bayestyper_reference_canon} -d {bayestyper_reference_decoy} \
					-p {threads} -o {params.out_prefix} > {log} 2>&1
		"""

# bayestyper genotype
rule run_bayestyper_genotype:
	input:
		kmc_pre="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_pre",
		kmc_suf="{results}/{sample}/bayestyper/{sample}_{fraction}/kmers/{sample}.kmc_suf",
		samples="{results}/{sample}/bayestyper/{sample}_{fraction}/{sample}.tsv",
		unit="{results}/{sample}/bayestyper/{sample}_{fraction}/clusters/bayestyper_unit_{unit_id}/variant_clusters.bin"
	output:
		genotypes="{results}/{sample}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
		kmer_coverage_file="{results}/{sample}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper_genomic_parameters.txt"
	log:
		"{results}/{sample}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper.log"
	params:
		cluster_data_dir="{results}/{sample}/bayestyper/{sample}_{fraction}/clusters/bayestyper_cluster_data",
		out_prefix="{results}/{sample}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper"
	threads: 24
	resources:
		mem_mb=50000,
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
	checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output[0]
	result = expand("{results}/{sample}/bayestyper/{sample}_{fraction}/genotype/bayestyper_unit_{unit_id}/bayestyper.vcf.gz",
								results=wildcards.results,
								sample=wildcards.sample,
								fraction=wildcards.fraction,
								unit_id=glob_wildcards(os.path.join(checkpoint_output, "bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
	return sorted(result)

rule bcftools_concat_units:
	input:
		aggregate_input
	output:
		"{results}/{sample}/bayestyper/{sample}_{fraction}_bayestyper_genotyping.vcf"
	log:
		"{results}/{sample}/bayestyper/{sample}_{fraction}_bayestyper_genotyping.log"
	conda:
		"../env/genotyping.yml"
	shell:
		"/usr/bin/time -v bcftools concat -a -o {output} {input} &> {log}"




########################################################
##################    run PanGenie    ##################
########################################################


# run pangenie
rule pangenie:
	input:
		reads=lambda wildcards: config['data'][wildcards.sample]['reads'],
        fasta=input_reference,
		vcf="{results}/{sample}/reference-panel/{sample}-all.vcf",
		#reads= outname_reads + "/{sample}/raw/{sample}-{fraction, [0-9.]+}.fastq",
		#vcf=lambda wildcards: config['data'][wildcards.sample]['graph'],
	output:
		"{results}/{sample}/pangenie/{sample}_{fraction}_pangenie_genotyping.vcf"
	log:
		"{results}/{sample}/pangenie/{sample}_{fraction}_pangenie_log.log"
	threads: 24
	params:
		out_prefix="{results}/{sample}/pangenie/{sample}_{fraction}_pangenie"
	resources:
		mem_mb=120000,
		runtime_hrs=9,
		runtime_min=59
	shell:
		"(/usr/bin/time -v {pangenie} -i {input.reads} -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g) &> {log}"



#########################################################
#################    run GraphTyper    ##################
#########################################################

# input VCF for graphtyper (biallelic, untypables removed since representation is changed)
rule graphtyper_preprocess:
	input:
		vcf="{results}/{sample}/biallelic/{sample}-chr{chrom}.vcf.gz",
		tbi="{results}/{sample}/biallelic/{sample}-chr{chrom}.vcf.gz.tbi"
	output:
		"{results}/{sample}/graphtyper/variants/variants_chr{chrom}_{variant}.vcf.gz"
	wildcard_constraints:
		chrom="X|Y|[0-9]+",
		variant="indel|sv"
	conda:
		"../env/genotyping.yml"
	shell:
		"""
		zcat {input.vcf} | python3 ../scripts/extract-varianttype.py {wildcards.variant} | bgzip > {output}
#		tabix -p vcf {output}
		"""

#### new graphtyper
## Remark: indel contains all variant that are not SVs, including SNPs, as well.
# genotype SNPs/indels/SVs (/usr/bin/time -v {graphtyper} genotype {input.fasta} --vcf={input.vcf_small} --sam={input.bam} --region=chr{wildcards.chrom} --no_decompose --verbose --output={params.dir_small} --threads={threads}) &> {log.small} 
rule graphtyper_genotype:
    input:
        vcf="{results}/{sample}/graphtyper/variants/variants_chr{chrom}_{variant}.vcf.gz",
        tbi="{results}/{sample}/graphtyper/variants/variants_chr{chrom}_{variant}.vcf.gz.tbi",
        bam='/vol/whopper/genotyping-experiments/genotyping/reads/{sample}/aligned/{sample}-{fraction}.chr{chrom}.bam',
        bai='/vol/whopper/genotyping-experiments/genotyping/reads/{sample}/aligned/{sample}-{fraction}.chr{chrom}.bam.bai',
        ##### Benefiting from alignments done for genotyping-experiments
        #bam=outname_reads + '/{sample}/aligned/{sample}-{fraction}.chr{chrom}.bam',
        #bai= outname_reads + '/{sample}/aligned/{sample}-{fraction}.chr{chrom}.bam.bai',
        fasta= input_reference
    output:
        vcf="{results}/{sample}/graphtyper/{variant}/{sample}_{fraction}_graphtyper_genotyping.chr{chrom, X|Y|[0-9]+}.vcf.gz",
        #tbi="{results}/{sample}/graphtyper/{variant}/{sample}_{fraction}_graphtyper_genotyping.chr{chrom, X|Y|[0-9]+}.vcf.gz.tbi"
    params:
        dir="{results}/{sample}/graphtyper/{variant}/{sample}_{fraction}_{chrom}"
    wildcard_constraints:
        variant="indel|sv"
    log:
        "{results}/{sample}/graphtyper/logs/{sample}_{fraction}_graphtyper_discovery.chr{chrom}.{variant}.log"
    conda:
        "../env/genotyping.yml"
    threads:
        24
    resources:
        mem_mb=30000,
        runtime_hrs=4,
        runtime_min=59
    shell:
        """
        if [ "{wildcards.variant}" == "sv" ]; then 
            (/usr/bin/time -v {graphtyper} genotype_sv {input.fasta} {input.vcf} --sam={input.bam} --region=chr{wildcards.chrom}  --output={params.dir} --threads={threads}) &> {log}
            bcftools concat -a {params.dir}/chr{wildcards.chrom}/*.vcf.gz | python3 ../scripts/graphtyper-postprocess.py {input.vcf} | bgzip -c > {output.vcf}
        else
            if [ "{wildcards.variant}" == "indel" ]; then
                (/usr/bin/time -v {graphtyper} genotype {input.fasta} --vcf={input.vcf} --sam={input.bam} --region=chr{wildcards.chrom} --no_decompose --verbose --output={params.dir} --threads={threads}) &> {log}
                bcftools concat -a {params.dir}/chr{wildcards.chrom}/*.vcf.gz | bgzip -c > {output.vcf}
            fi
        fi
        """
        #tabix -f -p vcf {output.vcf} > {output.tbi}

rule merge_vcfs_all_chromosomes:
    input:
       vcfs=lambda wildcards: expand("{{results}}/{{sample}}/graphtyper/{variant}/{{sample}}_{{fraction}}_graphtyper_genotyping.chr{chrom}.vcf.gz", chrom=chromosomes, variant=["indel", "sv"]),
       tbis=lambda wildcards: expand("{{results}}/{{sample}}/graphtyper/{variant}/{{sample}}_{{fraction}}_graphtyper_genotyping.chr{chrom}.vcf.gz.tbi", chrom=chromosomes, variant=["indel", "sv"])
    output: 
        "{results}/{sample}/graphtyper/{sample}_{fraction}_graphtyper_genotyping.vcf"
    shell: 
        """
        bcftools concat -a {input.vcfs} > {output}
        """


########################################################
##################      Evaluation    ##################
########################################################

def get_bed(wildcards):
	if wildcards.regions == "repeats-complex":
		return "{results}/bed/complex-rep.bed"
	elif wildcards.regions == "nonrep-complex":
		return "{results}/bed/complex-nonrep.bed"
	elif wildcards.regions == "repeats-simple":
		return "{results}/bed/simple-rep.bed"
	elif wildcards.regions == "nonrep-simple":
		return "{results}/bed/simple-nonrep.bed"
	elif wildcards.regions == "external":
		return config['data'][wildcards.sample]['external_bed']
	else:
		assert(False)

# prepare all regions
rule prepare_regions:
	input:
		fai=input_reference + '.fai',
		repeats=repeats_bed,
		complex=complex_bed
	output:
		tmp1=temp("{results}/bed/tmp1.txt"),
		tmp2=temp("{results}/bed/tmp2.txt"),
		tmp3=temp("{results}/bed/tmp4.txt"),
		complex_nonrep="{results}/bed/complex-nonrep.bed",
		complex_rep="{results}/bed/complex-rep.bed",
		simple_nonrep="{results}/bed/simple-nonrep.bed",
		simple_rep="{results}/bed/simple-rep.bed"
	conda:
		"../env/genotyping.yml"
	shell:
		"bash ../scripts/prepare-beds.sh {input.repeats} {input.complex} {input.fai} {output.tmp1} {output.tmp2} {output.tmp3} {output.complex_rep} {output.simple_rep} {output.complex_nonrep} {output.simple_nonrep}"


# prepare evaluation region (graph)
rule evaluation_region_graph:
	input:
		regions=get_bed,
		callable= lambda wildcards: config['data'][wildcards.sample]['bed']
	output:
		"{results}/{sample}/bed/graph/{sample}_{regions}.bed"
	wildcard_constraints:
		regions = "repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external"
	conda:
		"../env/genotyping.yml"
	shell:
		"bedtools intersect -a {input.regions} -b {input.callable} > {output}"


# prepare evaluation region (external)
rule evaluation_region_external:
	input:
		regions=get_bed,
		callable="{results}/{sample}/bed/graph/{sample}_external.bed"
	output:
		"{results}/{sample}/bed/external/{sample}_{regions}.bed"
	wildcard_constraints:
		regions = "repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external"
	conda:
		"../env/genotyping.yml"
	shell:
		"bedtools intersect -a {input.regions} -b {input.callable} > {output}"


# annotate genotyped VCFs and turn them to bi-allelic
rule convert_genotyping_to_biallelic:
	input:
		callset="{results}/{sample}/{method}/{sample}_{fraction}_{method}_genotyping.vcf",
		graph=lambda wildcards: "{results}/{sample}/biallelic/{sample}-all.vcf" if wildcards.method in ["graphtyper"] else "{results}/{sample}/reference-panel/{sample}-all.vcf",
		biallelic= lambda wildcards: config['data'][wildcards.sample]['biallelic']
	output:
		"{results}/{sample}/{method}/{sample}_{fraction}_{method}_genotyping-biallelic.vcf.gz"
	wildcard_constraints:
		method = "pangenie|bayestyper|graphtyper"
	conda:
		"../env/genotyping.yml"
	resources:
		mem_mb=20000,
		runtime_hrs=2
	shell:
		"""
		cat {input.callset} | python3 ../scripts/annotate.py {input.graph} | python3 ../scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' | bgzip > {output}
#		tabix -f -p vcf {output}
		"""


# normalize and annotate discovery sets and assign variant IDs
rule convert_discovery_to_biallelic:
	input:
		callset="{results}/{sample}/{method}/{sample}_{fraction}_{method}_discovery.vcf",
		biallelic= config['data']['full_callset'],
		reference = input_reference
	output:
		"{results}/{sample}/{method}/{sample}_{fraction}_{method}_discovery-biallelic.vcf.gz"
	conda:
		"../env/genotyping.yml"
	resources:
		mem_mb=20000,
		runtime_hrs=2
	shell:
		"""
		cat {input.callset} | python3 ../scripts/annotate.py {input.biallelic} | bgzip > {output}
#		tabix -p vcf {output}
		"""

# prepare external ground truth sets
rule prepare_ground_truth:
    input:
        callset= lambda wildcards: config['data'][wildcards.sample]['external'],
        truth = config['data']['full_callset'],
        reference=input_reference
    output:
        vcf="{results}/{sample}/external-truth/{sample}-truth.vcf.gz"
    params:
        working_dir = "/vol/whopper/genotyping-experiments/genotyping/"
    conda:
        "../env/genotyping.yml"
    resources:
        mem_mb=20000,
        runtime_hrs=2
    shell:
        """
        bcftools norm -f {input.reference} -m -any {input.callset} | bcftools sort -T {params.working_dir} | python3 ../scripts/annotate.py {input.truth} | bgzip > {output.vcf}
        """

# determine untypable IDs
rule remove_untypable:
	input:
		vcf="{path}{sample}{other}.vcf.gz",
		tbi="{path}{sample}{other}.vcf.gz.tbi",
		ids= lambda wildcards: config['data'][wildcards.sample]['untypable'] if wildcards.mode == "graph" else config['data'][wildcards.sample]['untypable_external']
	output:
		"{path}{sample}{other}-typable-{vartype}-{mode}.vcf.gz"
	wildcard_constraints:
		sample="|".join(samples),
		vartype="|".join(variants),
		mode = "external|graph"
	resources:
		mem_mb=20000,
		runtime_hrs=1
	shell:
		"""
		zcat {input.vcf} | python3 ../scripts/skip-untypable.py {input.ids} | python3 ../scripts/extract-varianttype.py {wildcards.vartype} | bgzip > {output}
#		tabix -p vcf {output}
		"""

rule prepare_all:
	input:
		vcf="{path}{sample}{other}.vcf.gz",
		tbi="{path}{sample}{other}.vcf.gz.tbi"
	output:
		"{path}{sample}{other}-all-{vartype}-{mode}.vcf.gz"
	wildcard_constraints:
		sample="|".join(samples),
		vartype="|".join(variants),
		mode= "external|graph"
	resources:
		mem_mb=20000,
		runtime_hrs=1
	shell:
		"""
		zcat {input.vcf} | python3 ../scripts/extract-varianttype.py {wildcards.vartype} | bgzip > {output}
#		tabix -p vcf {output}
		"""


########################## compute precision/recall ##############################

rule rtg_format:
	input:
		input_reference
	output:
		directory("{results}/evaluation/SDF")
	shell:
		'{rtg} format -o {output} {input}'


# precision-recall
rule vcfeval:
	input:
		callset="{results}/{sample}/{method}/{sample}_{fraction}_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz",
		callset_tbi="{results}/{sample}/{method}/{sample}_{fraction}_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz.tbi",
		baseline=lambda wildcards: "{results}/{sample}/external-truth/{sample}-truth-{variantset}-{vartype}-{mode}.vcf.gz" if wildcards.mode == "external" else config['data'][wildcards.sample]['truth'][:-4] + "-{variantset}-{vartype}-{mode}.vcf.gz",
		baseline_tbi=lambda wildcards: "{results}/{sample}/external-truth/{sample}-truth-{variantset}-{vartype}-{mode}.vcf.gz.tbi" if wildcards.mode == "external" else config['data'][wildcards.sample]['truth'][:-4] + "-{variantset}-{vartype}-{mode}.vcf.gz.tbi",
#		regions=lambda wildcards: "{results}/{sample}/bed/{sample}_{regions}.bed" if wildcards.mode == "graph" else config['data'][wildcards.sample]['external_bed'],
		regions="{results}/{sample}/bed/{mode}/{sample}_{regions}.bed",
		sdf="{results}/evaluation/SDF"
	output:
		summary="{results}/{sample}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_0/summary.txt"
	conda:
		"../env/genotyping.yml"
	wildcard_constraints:
		sample = "|".join(samples),
		mode = "external|graph",
		run_mode = "genotyping-biallelic",
		fraction = "|".join([str(f) for f in downsampling]),
		regions = "repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external",
		vartype = "|".join(variants),
		variantset = "typable|all"
	params:
		tmp = "{results}/{sample}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_0_tmp",
		outname = "{results}/{sample}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_0",
		which = lambda wildcards: "--all-records" if wildcards.run_mode == 'genotyping-biallelic' else ""
	resources:
		mem_mb=120000,
		runtime_hrs=0,
		runtime_min=20
	shell:
		"""
		if [ -d {params.tmp} ]; then
            rm -r {params.tmp}
        fi
        {rtg} vcfeval -b {input.baseline} -c {input.callset} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""


# precision-recall quality cutoff
rule vcfeval_cutoff:
	input:
		callset="{results}/{sample}/{method}/{sample}_{fraction}_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz",
		callset_tbi="{results}/{sample}/{method}/{sample}_{fraction}_{method}_{run_mode}-{variantset}-{vartype}-{mode}.vcf.gz.tbi",
		baseline=lambda wildcards: "{results}/{sample}/external-truth/{sample}-truth-{variantset}-{vartype}-{mode}.vcf.gz" if wildcards.mode == "external" else config['data'][wildcards.sample]['truth'][:-4] + "-{variantset}-{vartype}-{mode}.vcf.gz",
		baseline_tbi=lambda wildcards: "{results}/{sample}/external-truth/{sample}-truth-{variantset}-{vartype}-{mode}.vcf.gz.tbi" if wildcards.mode == "external" else config['data'][wildcards.sample]['truth'][:-4] + "-{variantset}-{vartype}-{mode}.vcf.gz.tbi",
#		regions=lambda wildcards: "{results}/{sample}/bed/{sample}_{regions}.bed" if wildcards.mode == "graph" else config['data'][wildcards.sample]['external_bed'],
		regions="{results}/{sample}/bed/{mode}/{sample}_{regions}.bed",
		sdf="{results}/evaluation/SDF"
	output:
		tmp_vcf=temp("{results}/{sample}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}_coverage-{fraction}_{regions}_{vartype}_qual_{qual}.vcf.gz"),
		tmp_tbi=temp("{results}/{sample}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}_coverage-{fraction}_{regions}_{vartype}_qual_{qual}.vcf.gz.tbi"),
		summary="{results}/{sample}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt"
	conda:
		"../env/genotyping.yml"
	wildcard_constraints:
		sample = "|".join(samples),
		mode = "external|graph",
		run_mode = "genotyping-biallelic",
		fraction = "|".join([str(f) for f in downsampling]),
		regions = "repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external",
		vartype = "|".join(variants),
		qual = "200",
		variantset = "typable|all"
	params:
		tmp = "{results}/{sample}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}_tmp",
		outname = "{results}/{sample}/evaluation/precision-recall-{variantset}/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}",
		which = lambda wildcards: "--all-records" if wildcards.run_mode == 'genotyping-biallelic' else ""
	resources:
		mem_mb=120000,
		runtime_hrs=0,
		runtime_min=20
	shell:
		"""
		bcftools view -i 'FMT/GQ>={wildcards.qual}' -O z {input.callset} > {output.tmp_vcf}
		tabix -p vcf {output.tmp_vcf}
		{rtg} vcfeval -b {input.baseline} -c {output.tmp_vcf} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {input.regions} {params.which} --Xmax-length 30000 > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""


################################## compute genotype concordance ########################################


# determine the variants that went into re-typing per category
rule collected_typed_variants:
	input:
		graph= lambda wildcards: config['data'][wildcards.sample]['biallelic'] + ".gz",
#		regions = lambda wildcards: "{results}/{sample}/bed/{sample}_{regions}.bed" if wildcards.mode == "graph" else config['data'][wildcards.sample]['external_bed']
		regions="{results}/{sample}/bed/{mode}/{sample}_{regions}.bed"
	output:
		"{results}/{sample}/genotyped-ids/{mode}_{regions}_{vartype}.tsv"
	conda:
		"../env/genotyping.yml"
	wildcard_constraints:
		sample = "|".join(samples),
		mode = "external|graph",
		regions = "repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external",
		vartype = "|".join(variants)
	shell:
		"zcat {input.graph} | python3 ../scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 ../scripts/get_ids.py > {output}"


# compute concordances
rule genotype_concordances:
	input:
		callset="{results}/{sample}/{method}/{sample}_{fraction}_{method}_{run_mode}-typable-{vartype}-{mode}.vcf.gz",
		callset_tbi="{results}/{sample}/{method}/{sample}_{fraction}_{method}_{run_mode}-typable-{vartype}-{mode}.vcf.gz.tbi",
		baseline=lambda wildcards: "{results}/{sample}/external-truth/{sample}-truth-typable-{vartype}-{mode}.vcf.gz" if wildcards.mode == "external" else config['data'][wildcards.sample]['truth'][:-4] + "-typable-{vartype}-{mode}.vcf.gz",
		baseline_tbi=lambda wildcards: "{results}/{sample}/external-truth/{sample}-truth-typable-{vartype}-{mode}.vcf.gz.tbi" if wildcards.mode == "external" else config['data'][wildcards.sample]['truth'][:-4] + "-typable-{vartype}-{mode}.vcf.gz.tbi",
#		regions=lambda wildcards: "{results}/{sample}/bed/{sample}_{regions}.bed" if wildcards.mode == "graph" else config['data'][wildcards.sample]['external_bed'],
		regions="{results}/{sample}/bed/{mode}/{sample}_{regions}.bed",
		typed_ids = "{results}/{sample}/genotyped-ids/{mode}_{regions}_{vartype}.tsv"
	output:
		tmp_vcf1=temp("{results}/{sample}/evaluation/concordance/{mode}/{method}-{run_mode}_coverage-{fraction}_{regions}_{vartype}_qual_{qual}_base.vcf"),
		tmp_vcf2=temp("{results}/{sample}/evaluation/concordance/{mode}/{method}-{run_mode}_coverage-{fraction}_{regions}_{vartype}_qual_{qual}_call.vcf"),
		summary="{results}/{sample}/evaluation/concordance/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt"
	conda:
		"../env/genotyping.yml"
	wildcard_constraints:
		sample = "|".join(samples),
		mode = "external|graph",
		run_mode = "genotyping-biallelic",
		fraction = "|".join([str(f) for f in downsampling]),
		regions = "repeats-complex|repeats-simple|nonrep-complex|nonrep-simple|external",
		vartype = "|".join(variants),
		qual= "0|200"
	params:
		which = lambda wildcards: "" if wildcards.run_mode == 'genotyping-biallelic' else "--only_pass"
	log:
		"{results}/{sample}/evaluation/concordance/{mode}/{method}-{run_mode}/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.log"
	resources:
		mem_mb=20000,
		runtime_hrs=0,
		runtime_min=20
	shell:
		"""
		bedtools intersect -header -a {input.baseline} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf1}
		bedtools intersect -header -a {input.callset} -b {input.regions} -u -f 0.5 | bgzip > {output.tmp_vcf2}
		python3 ../scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual {wildcards.qual} {params.which} 2> {log} 1> {output.summary}
		"""



########################################################
##################       Plotting     ##################
########################################################

def plot_input(wildcards):
	output = []
	variant_list = ['indel', 'sv', 'large-insertion', 'large-deletion', 'large'] if (wildcards.run_mode == "discovery-biallelic") or (wildcards.mode != "graph") else variants
	all_regions = [wildcards.regions] if 'external' in wildcards.regions else [wildcards.regions + '-simple', wildcards.regions + '-complex']
	metric = 'precision-recall-typable' if wildcards.metric == 'fscore' else wildcards.metric
	for reg in all_regions:
		if wildcards.run_mode == "genotyping-biallelic":
			for fraction in downsampling:
				for var in variant_list:
					for method in ['pangenie'] + other_methods:	
						output.append("{results}/{sample}/evaluation/{metric}/{mode}/{method}-genotyping-biallelic/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt".format(
							results=wildcards.results,
							vartype=var,
							sample=wildcards.sample,
							mode=wildcards.mode,
							method=method,
							fraction=fraction,
							regions=reg,
							qual="0",
							metric=metric
							)
						)
					output.append("{results}/{sample}/evaluation/{metric}/{mode}/pangenie-genotyping-biallelic/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt".format(
							results=wildcards.results,
							vartype=var,
							sample=wildcards.sample,
							mode=wildcards.mode,
							fraction=fraction,
							regions=reg,
							qual="200",
							metric=metric
							)
	
						)
		elif wildcards.run_mode == "discovery-biallelic":
			for fraction in downsampling:
					for var in variant_list:
						output.append("{results}/{sample}/evaluation/{metric}/{mode}/discovery-biallelic/coverage-{fraction}_{regions}_{vartype}/qual_{qual}/summary.txt".format(
							results=wildcards.results,
							sample=wildcards.sample,
							vartype=var,
							mode=wildcards.mode,
							fraction=fraction,
							regions=reg,
							qual="0",
							metric=metric
							)
						)
		else:
			assert(False)
	return output


# plot precision/recall or genotype concordance over all coverages
rule plot_results:
	input:
		plot_input
	output:
		"{results}/{sample}/evaluation/{metric}/{mode}/{run_mode}_{regions}_{metric}.pdf"
	conda:
		"../env/genotyping.yml"
	wildcard_constraints:
		run_mode = "genotyping-biallelic",
		regions = "repeats|nonrep|external",
		metric = "precision-recall-typable|precision-recall-all|concordance|fscore"
	params:
		folder = "{results}/",
		script = lambda wildcards: metric_to_script[wildcards.metric],
		coverages = lambda wildcards: [max(downsampling)] if wildcards.metric == "fscore" else downsampling,
		run_mode = lambda wildcards: wildcards.run_mode if 'genotyping-biallelic' == wildcards.run_mode else [wildcards.run_mode, 'genotyping-biallelic'],
		varset = lambda wildcards: "-variantset all" if wildcards.metric == "precision-recall-all" else ""
	shell:
		"python3 {params.script} {wildcards.mode} {wildcards.regions} -run_mode {params.run_mode} -folder {params.folder} -coverages {params.coverages} -outfile {output} -sample {wildcards.sample} {params.varset}"

