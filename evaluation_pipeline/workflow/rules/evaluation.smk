kmc=config['programs']['kmc']
bcftools=config['programs']['bcftools']
bayestyper=config['programs']['bayestyper']
bayestyper_tools=config['programs']['bayestyper_tools']
graphtyper=config['programs']['graphtyper']
rtg=config['programs']['rtg']
truvari = config['programs']['truvari']
vartypes_leave1out = [var for c in config["callsets"].keys() for t in config["callsets"][c]["leave1out"].keys() for var in config["callsets"][c]["leave1out"][t]["vartype"]]
vartypes_external = [var for c in config["callsets"].keys() for t in config["callsets"][c]["external"].keys() for var in config["callsets"][c]["external"][t]["vartype"]]

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
		"(cat {input.genotyped_vcf} | python workflow/scripts/annotate.py {input.panel_vcf} | python workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' | bgzip -c > {output}) &> {log}"

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
		{bcftools} view {input.vcf} | {params.untypables} python workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
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
		{bcftools} view {input.vcf} | {params.untypables} python workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output.vcf}
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
		fixed_vcf="evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall/fixed.vcf.gz",
		summary="evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "snp-indel|snp|indels",
        region = "all|multi|bi",
		filter = "all|typable"
	log:
		"logs/evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall.log"
	params:
		tmp = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall_tmp",
		outname = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall",
		bed = lambda wildcards: f"--evaluation-regions preprocessing/{wildcards.callset}/{wildcards.sample}/regions/{wildcards.region}-bubbles.bed" if wildcards.region != 'all' else ""
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
	threads: 24
	shell:
		"""
		{bcftools} view {input.callset} --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.fixed_vcf}
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
		fixed_vcf = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall/fixed.vcf.gz",
		summary_json = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall/summary.json",
		summary_txt = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "sv|large-deletion|large-insertion",
        region = "all|multi|bi",
		filter = "all|typable"
	log:
		"logs/evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall.log"
	params:
		tmp = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall_tmp",
		outname = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/precision-recall",
		bed = lambda wildcards: f"--includebed preprocessing/{wildcards.callset}/{wildcards.sample}/regions/{wildcards.region}-bubbles.bed" if wildcards.region != 'all' else ""
	resources:
		mem_total_mb = 20000,
		runtime_hrs = 0,
		runtime_min = 40
	threads: 1
	shell:
		"""
		{bcftools} view {input.callset} --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.fixed_vcf}
		tabix -p vcf {output.fixed_vcf}
		(/usr/bin/time -v {truvari} bench -b {input.baseline} -c {output.fixed_vcf} -f {input.reference} -o {params.tmp} --pick multi {params.bed} -r 2000 --no-ref a -C 2000 --passonly) &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		python workflow/scripts/parse_json_to_txt.py {output.summary_json} > {output.summary_txt}
		"""
		


################ Plots ###################

# determine the variants that went into re-typing per category
rule collect_typed_variants:
	input:
		callset = "preprocessing/{callset}/{sample}/input-panel/panel_bi.vcf",
		regions= region_to_bed,
		ids="preprocessing/{callset}/{sample}/{pipeline}/untypables/untypables.tsv"
	output:
		"preprocessing/{callset}/{sample}/{pipeline}/genotyped-ids/typed_{region}_{vartype}.tsv"
		# "results/leave-one-out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		region = "bi|multi|all",
		vartype = "snp-indel|sv|snp|indels|large-deletion|large-insertion"
	resources:
		mem_total_mb=50000
	priority: 1
	shell:
		"cat {input.callset} | python workflow/scripts/skip-untypable.py {input.ids} | python workflow/scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python workflow/scripts/get_ids.py > {output}"


# compute concordances
rule genotype_concordances:
	input:
		callset = "genotyping/{callset}/{sample}/{pipeline}/{version}/{coverage}/callset-{filter}_{vartype}.vcf.gz",
		# callset_tbi = "genotyping/{callset}/{sample}/{pipeline}/{version}/{coverage}/callset-{filter}_{vartype}.vcf.gz.tbi",
		baseline = "genotyping/{callset}/{sample}/{pipeline}/truthset/truthset-{filter}_{vartype}.vcf.gz",
		# baseline_tbi = "genotyping/{callset}/{sample}/{pipeline}/truthset/truthset-{filter}_{vartype}.vcf.gz.tbi",
		region = region_to_bed,
		typed_ids = "preprocessing/{callset}/{sample}/{pipeline}/genotyped-ids/typed_{region}_{vartype}.tsv"
	output:
		tmp_vcf1 = temp("genotyping/{callset}/{sample}/{pipeline}/{version}/{coverage}/concordance/{region}-{filter}-{vartype}_base.vcf"),
		tmp_vcf2 = temp("genotyping/{callset}/{sample}/{pipeline}/{version}/{coverage}/concordance/{region}-{filter}-{vartype}_call.vcf"),
		summary = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/concordance/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		region = "bi|multi|all",
		vartype = "snp-indel|sv|snp|indels|large-deletion|large-insertion"
	log:
		"logs/evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/concordance.log"
	resources:
		mem_total_mb = 40000,
		runtime_hrs = 0,
		runtime_min = 40
	shell:
		"""
		bedtools intersect -header -a {input.baseline} -b {input.region} -u -f 0.5 | bgzip > {output.tmp_vcf1}
		bedtools intersect -header -a {input.callset} -b {input.region} -u -f 0.5 | bgzip > {output.tmp_vcf2}
		python workflow/scripts/genotype-evaluation.py {output.tmp_vcf1} {output.tmp_vcf2} {input.typed_ids} --qual 0 2> {log} 1> {output.summary}
		"""

# collect results across all metrics
rule collect_results_evaluation_leave1out:
	input:
		"evaluation/{callset}/{sample}/leave1out/{version}-{coverage}-{filter}-{region}-{vartype}/{metric}/summary.txt"
	output:
		"evaluation/{callset}/{sample}/leave1out/plots/{metric}-{version}-{coverage}-{filter}-{region}_{vartype}.tsv"
	params:
		outfile = "evaluation/{callset}/{sample}/leave1out/plots/{metric}-{version}-{coverage}-{filter}-{region}_{vartype}",
		folder = "evaluation/{callset}/{sample}/leave1out", 
	wildcard_constraints:
		region = "bi|multi|all",
		filter = "typable|all",
		metric = "concordance|precision-recall",
		vartype = "snp|indels|large-deletion|large-insertion"
	log:
		"logs/evaluation/{callset}/{sample}/leave1out/plots/{metric}-{version}-{coverage}-{filter}-{region}_{vartype}.log"
	shell:
		"(python workflow/scripts/collect-results.py {wildcards.metric} {wildcards.version} {wildcards.sample} {wildcards.coverage} {wildcards.region} {wildcards.filter} -variants {wildcards.vartype} -folder {params.folder} -outfile {params.outfile}) &> {log}"

rule collect_results_evaluation_external:
	input:
		"evaluation/{callset}/{sample}/external/{version}-{coverage}-{filter}-{region}-{vartype}/{metric}/summary.txt"
	output:
		"evaluation/{callset}/{sample}/external/plots/{metric}-{version}-{coverage}-{filter}-{region}_{vartype}.tsv"
	params:
		outfile = "evaluation/{callset}/{sample}/external/plots/{metric}-{version}-{coverage}-{filter}-{region}_{vartype}",
		folder = "evaluation/{callset}/{sample}/external", 
	wildcard_constraints:
		region = "bi|multi|all",
		filter = "typable|all",
		metric = "concordance|precision-recall",
		vartype = "snp-indel|sv"
	log:
		"logs/evaluation/{callset}/{sample}/external/plots/{metric}-{version}-{coverage}-{filter}-{region}_{vartype}.log"
	shell:
		"(python workflow/scripts/collect-results.py {wildcards.metric} {wildcards.version} {wildcards.sample} {wildcards.coverage} {wildcards.region} {wildcards.filter} -variants {wildcards.vartype} -folder {params.folder} -outfile {params.outfile}) &> {log}"


# plot results of different subsampling runs
rule plotting_versions_leave1out:
	input:
		lambda wildcards: expand("evaluation/{{callset}}/{{sample}}/leave1out/plots/{m}-{version}-{{coverage}}-{{filter}}-{{region}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', version=versions_to_run, vartype=vartypes_leave1out)
	output:
		"evalplots/{callset}/{sample}/leave1out/comparison-versions/{metric}-{coverage}-{filter}-{region}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall|untyped",
		regions="bi|multi|all",
		filter="typable|all"
	log:
		"logs/evalplots/{callset}/{sample}/leave1out/comparison-versions/{metric}-{coverage}-{filter}-{region}.log"
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([v + '-' + wildcards.coverage + '-' + wildcards.filter + '-' + wildcards.region for v in versions_to_run]),
		vartypes = vartypes_leave1out
	shell:
		"(python workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric} -vartypes {params.vartypes}) &> {log}"



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
		"python workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric concordance-vs-untyped"
		


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
		"python workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"
