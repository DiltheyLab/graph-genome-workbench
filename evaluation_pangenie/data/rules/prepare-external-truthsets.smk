#samples = config['dataset']['assemblies']['contigs'].keys()
#samples_parents = [s for s in samples if not s in config['dataset']['assemblies']['trios']]
samples = config["samples_HGSVC"]
datadir = config['outdir'] 

scripts=config['scripts']

truthset_to_panel = {
		'giab-small': datadir  + "multisample-vcfs/assemblies-noNA12878-biallelic-filtered.vcf",
		'syndip': datadir  + "multisample-vcfs/assemblies-all-biallelic-filtered.vcf",
		'giab-sv': datadir + "multisample-vcfs/assemblies-noNA24385-biallelic-filtered.vcf", 
		'GIAB': "/vol/whopper/graph-genome-workbench/evaluation_pangenie/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf"
	}

truthset_to_truth = {
		'giab-small': datadir  + "multisample-vcfs/assemblies-all-biallelic-filtered.vcf",
		'syndip': datadir  + "multisample-vcfs/assemblies-all-biallelic-filtered.vcf",
		'giab-sv': datadir + "multisample-vcfs/assemblies-all-biallelic-filtered.vcf",
		'GIAB':  "/vol/whopper/graph-genome-workbench/evaluation_pangenie/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf"
}

truthset_to_samples = {
		'giab-small': [s for s in samples_parents if not s=="NA12878"],
		'syndip': samples_parents,
		'giab-sv': [s for s in samples_parents if not s=="NA24385"],
		'GIAB': [s for s in samples if not s=="NA24385"],
	}



# remove variants containing alleles different from: ATCG
rule filter_syndip:
	input:
		lambda wildcards: config['dataset'][wildcards.truthset]
	output:
		config['outdir'] + "{truthset}/{truthset}-filtered.vcf.gz"
	conda:
		"../../env/calling.yml"
	shell:
		"""
		zcat {input} | python3 {scripts}/filter-calls.py | bgzip -c > {output}
		tabix -f -p vcf {output}
		"""

# normalize the ground truth variants and annotate them with IDs from the panel VCF
# (matching variants will be assigned the same IDs)
rule prepare_ground_truth:
    input:
        callset= config['outdir'] + "{truthset}/{truthset}-filtered.vcf.gz",
        truth = lambda wildcards: truthset_to_truth[wildcards.truthset],
        reference= config['reference']
    output:
        vcf= config['outdir'] + "{truthset}/{truthset}-annotated.vcf.gz"
    params:
        working_dir = "/vol/whopper/genotyping-experiments/genotyping/"
    conda:
        "../../env/genotyping.yml"
    wildcard_constraints:
        truthset = "GIAB",
        #truthset = "giab-small|giab-sv"
    threads: 24
    resources:
        mem_mb=20000,
        runtime_hrs=2
    shell:
        """
        bcftools norm -f {input.reference} -m -any --threads {threads} {input.callset} | bcftools sort -T {params.working_dir} | python3 {scripts}/annotate.py {input.truth} | bgzip -c > {output.vcf}
        tabix -f -p vcf {output.vcf}
        """


# compare truth set sample to any panel sample and determine false negatives
# i.e. variants only in panel sample (to be later intersected across all samples)
rule extract_sample:
	input:
		truthset = config['outdir'] + "{truthset}/{truthset}-annotated.vcf.gz",
		panel = lambda wildcards: truthset_to_panel[wildcards.truthset],
		sdf=config['outdir'] + 'SDF'
	output:
		sample_vcf=config['outdir'] + '{truthset}/samples/{sample}.vcf.gz',
		sample_vcf_tbi=config['outdir'] + '{truthset}/samples/{sample}.vcf.gz.tbi',
		fn = config['outdir'] + '{truthset}/samples/{sample}/vcfeval/fn.vcf.gz'
	wildcard_constraints:
		truthset="GIAB",
		#truthset="giab-small|giab-sv"
	conda:
		"../../env/genotyping.yml"
	params:
		tmp = config['outdir'] + '{truthset}/samples/{sample}/vcfeval_temp',
		outname = config['outdir'] + '{truthset}/samples/{sample}/vcfeval',
		sdf= config['outdir'] + 'SDF'
	log:
		config['outdir'] + '{truthset}/samples/{sample}/vcfeval.log'
	resources:
		mem_mb=120000,
		runtime_hrs=0,
		runtime_min=59
	shell:
		"""
		bcftools view --samples {wildcards.sample} {input.panel} | bgzip -c > {output.sample_vcf}
		tabix -f -p vcf {output.sample_vcf}
		/home/ubuntu/rtg-tools-3.12.1/rtg vcfeval -b {input.truthset} -c {output.sample_vcf} -t {params.sdf} -o {params.tmp} --squash-ploidy --ref-overlap --all-records --Xmax-length 30000 &> {log}
		mkdir -p {params.outname}
        mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""

# prepare reference for vcfeval
rule rtg_format:
	input:
		config['reference']
	output:
		directory(config['outdir'] + 'SDF')
	conda:
		"../../env/genotyping.yml"
	shell:
		'/home/ubuntu/rtg-tools-3.12.1/rtg format -o {output} {input}'



# determine variants unique to a truthset
rule determine_unique:
	input:
		truthset = config['outdir'] + "{truthset}/{truthset}-annotated.vcf.gz",
		samples = lambda wildcards: expand(config['outdir'] + '{{truthset}}/samples/{sample}/vcfeval/fn.vcf.gz', sample = truthset_to_samples[wildcards.truthset] )
	output:
		unique_tsv=config['outdir'] + "{truthset}/{truthset}-unique.tsv",
		unique_vcf=config['outdir'] + "{truthset}/{truthset}-unique.vcf"
	conda:
		"../../env/genotyping.yml"
	params:
		n_files = lambda wildcards: len(truthset_to_samples[wildcards.truthset])
	resources:
		mem_mb=30000,
		runtime_hrs=0,
		runtime_min=30
	shell:
		"""
		bcftools isec -n={params.n_files} -w1 {input.samples}  > {output.unique_vcf}
		grep -v '#' {output.unique_vcf} | cut -f 3  > {output.unique_tsv}
		"""


################################################
#  Compute some statistics 
################################################

rule count_unique_type:
	input:
		total=config['outdir'] + "{truthset}/{truthset}-annotated.vcf.gz",
		unique=config['outdir'] + "{truthset}/{truthset}-unique.vcf"
	output:
		total=config['outdir'] + "{truthset}/{truthset}-annotated-{vartype}.vcf",
		unique=config['outdir'] + "{truthset}/{truthset}-unique-{vartype}.vcf"
	conda:
		"../../env/genotyping.yml"
	wildcard_constraints:
		vartype="indel|sv"
	shell:
		"""
		bcftools view {input.total} | python3 {scripts}/extract-varianttype.py {wildcards.vartype} > {output.total}
		bcftools view {input.unique} | python3 {scripts}/extract-varianttype.py {wildcards.vartype} > {output.unique}
		"""

rule callset_precision_recall:
	input:
		callset = "/vol/whopper/graph-genome-workbench/evaluation_pangenie/downloaded/vcf/HGSVC-GRCh38/graph_annotated_onlyNA24385_biallelic.vcf.gz", 
		truth = config['outdir'] + "{truthset}/{truthset}-annotated.vcf.gz",
		sdf = config['outdir'] + 'SDF',
		truth_bed =  config['truth_bed'],
		callset_bed = datadir + "bed/NA24385_callable.bed"
		### three arguments remove. Check it out in original if needed
	output:
		summary=config['outdir'] + "{truthset}/comparison/{truthset}-calls-{vartype}/summary.txt",
		tmp1 = temp(config['outdir'] + "{truthset}/comparison/{truthset}-calls-callset-{vartype}.vcf.gz"),
		tmp2 = temp(config['outdir'] + "{truthset}/comparison/{truthset}-calls-truthset-{vartype}.vcf.gz"),
		tmp_bed = temp(config['outdir'] + "{truthset}/comparison/{truthset}-calls-{vartype}.bed")
	params:
		tmp= config['outdir'] + "{truthset}/comparison/{truthset}-calls-{vartype}_tmp",
		outname = config['outdir'] + "{truthset}/comparison/{truthset}-calls-{vartype}"
	conda:
		"../../env/genotyping.yml"
	wildcard_constraints:
		truthset = "GIAB"
	resources:
		mem_mb=120000
	shell:
		"""
		bedtools intersect -a {input.truth_bed} -b {input.callset_bed} > {output.tmp_bed}
		bcftools view {input.callset} | python3 {scripts}/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.tmp1}
		bcftools view {input.truth} | python3 {scripts}/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.tmp2}
		tabix -f -p vcf {output.tmp1}
		tabix -f -p vcf {output.tmp2}
		/home/ubuntu/rtg-tools-3.12.1/rtg vcfeval -b {output.tmp2} -c {output.tmp1} -t {input.sdf} -o {params.tmp} --ref-overlap --evaluation-regions {output.tmp_bed} --all-records --Xmax-length 30000 > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""
