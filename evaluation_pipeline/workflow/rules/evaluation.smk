kmc=config['programs']['kmc']
bayestyper=config['programs']['bayestyper']
bayestyper_tools=config['programs']['bayestyper_tools']
graphtyper=config['programs']['graphtyper']
rtg=config['programs']['rtg']
truvari = config['programs']['truvari']

###########################################################
#################      Evaluation        ##################
###########################################################

def get_panel_vcf(wildcards):
	if wildcards.version == 'graphtyper':
		return "preprocessing/{callset}/{sample}/input-panel/panel_bi.vcf.gz"
	else:
		return "preprocessing/{callset}/{sample}/input-panel/panel_multi_norm.vcf.gz"

# convert genotyped VCF to biallelic representation
rule convert_genotypes_to_biallelic:
	input:
		genotyped_vcf = "genotyping/{callset}/{sample}/{version}/{coverage}/genotyping.vcf",
        panel_vcf = get_panel_vcf,
		biallelic = lambda wildcards: config['callsets'][wildcards.callset]['bi']
	output:
		"genotyping/{callset}/{sample}/{version}/{coverage}/genotyping-biallelic.vcf.gz"
	log:
		"logs/genotyping/{callset}/{sample}/{version}/{coverage}/genotyping-biallelic.log"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=30000
	priority: 1
	shell:
		"(cat {input.genotyped_vcf} | python3 workflow/scripts/annotate.py {input.panel_vcf} | python3 workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' | bgzip -c > {output}) &> {log}"

## remove untypables and split variants
rule postprocessing_callset:
	input:
		vcf = "genotyping/{callset}/{sample}/{version}/{coverage}/genotyping-biallelic.vcf.gz",
		untypable_ids = "preprocessing/{callset}/{sample}/{pipeline}/untypables/untypables.tsv"
	output:
		vcf = "genotyping/{callset}/{sample}/{pipeline}/{version}/{coverage}/callset-{filter}_{vartype}.vcf.gz"
	wildcard_constraints:
		vartype = "snp-indel|sv|snp|indels|large-deletion|large-insertion",
		filter="typable|all"
	resources:
		mem_total_mb=30000
	conda:
		"../envs/genotyping.yml"
	params:
		untypables = lambda wildcards: f" python workflow/scripts/skip-untypable.py preprocessing/{wildcards.callset}/{wildcards.sample}/{wildcards.pipeline}/untypables/untypables.tsv |" if wildcards.filter == 'typable' else ""
	shell:
		"""
		bcftools view {input.vcf} | {params.untypables} python workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

## remove untypables and split variants
rule postprocessing_truthset:
	input:
		vcf = "preprocessing/{callset}/{sample}/{pipeline}/truthset/truthset-biallelic.vcf.gz",
		untypable_ids = "preprocessing/{callset}/{sample}/{pipeline}/untypables/untypables.tsv"
	output:
		vcf = "genotyping/{callset}/{sample}/{pipeline}/truthset/truthset-{filter}_{vartype}.vcf.gz"
	wildcard_constraints:
		vartype = "snp-indel|sv|snp|indels|large-deletion|large-insertion",
		filter="typable|all"
	resources:
		mem_total_mb=30000
	conda:
		"../envs/genotyping.yml"
	params:
		untypables = lambda wildcards: f" python workflow/scripts/skip-untypable.py preprocessing/{wildcards.callset}/{wildcards.sample}/{wildcards.pipeline}/untypables/untypables.tsv |" if wildcards.filter == 'typable' else ""
	shell:
		"""
		bcftools view {input.vcf} | {params.untypables} python workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

def region_to_bed(wildcards):
	if wildcards.region == "bi":
		return "preprocessing/{callset}/{sample}/regions/bi-bubbles.bed".format(callset=wildcards.callset, sample=wildcards.sample)
	elif wildcards.region == "multi":
		return "preprocessing/{callset}/{sample}/regions/multi-bubbles.bed".format(callset=wildcards.callset, sample=wildcards.sample)
	elif wildcards.region == "all":
		return "resources/ucsc-simple-repeats.bed"
	else:
		assert(False)

# compute precision/recall for small variants
rule evaluate_with_vcfeval:
	input:
		callset = "genotyping/{callset}/{sample}/{pipeline}/{version}/{coverage}/callset-{filter}_{vartype}.vcf.gz",
		baseline = "genotyping/{callset}/{sample}/{pipeline}/truthset/truthset-{filter}_{vartype}.vcf.gz",
		regions = region_to_bed,
		sdf = "preprocessing/{callset}/SDF",
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai']
	output:
		fixed_vcf="evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/fixed.vcf.gz",
		summary="evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "snp-indel|snp|indels",
        region = "all|multi|bi",
		filter = "all|typable"
	log:
		"logs/evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}.log"
	params:
		tmp = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}_tmp",
		outname = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}",
		bed = lambda wildcards: f"--evaluation-regions preprocessing/{wildcards.callset}/{wildcards.sample}/regions/{wildcards.region}-bubbles.bed" if wildcards.region != 'all' else ""
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
	threads: 24
	shell:
		"""
		bcftools view {input.callset} --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.fixed_vcf}
		tabix -p vcf {output.fixed_vcf}
		(/usr/bin/time -v {rtg} vcfeval -b {input.baseline} -c {output.fixed_vcf} -t {input.sdf} -o {params.tmp} {params.bed} --ref-overlap --all-records --Xmax-length 30000 --threads {threads} > {output.summary}.tmp) &> {log}
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""

# compute precision/recall for SVs
rule evaluate_with_truvari:
	input:
		callset = "genotyping/{callset}/{sample}/{pipeline}/{version}/{coverage}/callset-{filter}_{vartype}.vcf.gz",
		baseline = "genotyping/{callset}/{sample}/{pipeline}/truthset/truthset-{filter}_{vartype}.vcf.gz",
		regions = region_to_bed,
		reference = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai']
	output:
		fixed_vcf = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/fixed.vcf.gz",
		summary = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/summary.json"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "sv|large-deletion|large-insertion",
        region = "all|multi|bi",
		filter = "all|typable"
	log:
		"logs/evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}.log"
	params:
		tmp = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}_tmp",
		outname = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}",
		bed = lambda wildcards: f"--includebed preprocessing/{wildcards.callset}/{wildcards.sample}/regions/{wildcards.region}-bubbles.bed" if wildcards.region != 'all' else ""
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 0,
		runtime_min = 40
	threads: 1
	shell:
		"""
		bcftools view {input.callset} --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.fixed_vcf}
		tabix -p vcf {output.fixed_vcf}
		(/usr/bin/time -v {truvari} bench -b {input.baseline} -c {output.fixed_vcf} -f {input.reference} -o {params.tmp} --pick multi {params.bed} -r 2000 --no-ref a -C 2000 --passonly > {output.summary}.tmp) &> {log}
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""