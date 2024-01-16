kmc=config['programs']['kmc']
bayestyper=config['programs']['bayestyper']
bayestyper_tools=config['programs']['bayestyper_tools']
graphtyper=config['programs']['graphtyper']
rtg=config['programs']['rtg']
truvari = config['programs']['truvari']

# parameters
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
ref_chromosomes = {}
for c in config['callsets'].keys():
    ref_chromosomes[c] = []
    for line in open(config['callsets'][c]['reference_fai'], 'r'):
        fields = line.split()
        chr_number = fields[0].split('chr')[1].strip()
        if chr_number in chromosomes:
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

callsets_leave_one_out = [s for s in config['callsets'].keys()]
allowed_variants = ['snp', 'indels', 'large-insertion', 'large-deletion']
coverages_leave_one_out = ['full'] + [c for c in config['downsampling']]
versions_leave_one_out = [v for v  in config['pangenie'].keys()] + [v for v in config['pangenie-modules'].keys()] + ['bayestyper', 'graphtyper'] # , 'graphtyper'
callsets = [s for s in config['callsets'].keys()]


########################################################
##################    Util functions    ################
########################################################
rule compress_vcf:
	input:
		"genotyping/leave1out/{filename}.vcf"
	output:
		vcf = "genotyping/leave1out/{filename}.vcf.gz",
		tbi = "genotyping/leave1out/{filename}.vcf.gz.tbi"
	priority: 1
	shell:
		"""
		bgzip -c {input} > {output.vcf}
		tabix -p vcf {output.vcf}
		"""


###########################################################
##############    Leave1out Evaluation      ###############
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
		#vartype = "snp-indel|sv|snp|indels|large-deletion|large-insertion",
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
		fixed_vcf="evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/fixed.vcf.gz",
		summary="evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}/summary.json"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype = "sv|large-deletion|large-insertion",
		filter = "all|typable"
	log:
		"logs/evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}.log"
	params:
		tmp = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}_tmp",
		outname = "evaluation/{callset}/{sample}/{pipeline}/{version}-{coverage}-{filter}-{region}-{vartype}",
		bed = lambda wildcards: f"--includebed preprocessing/{wildcards.callset}/{wildcards.sample}/regions/{wildcards.region}-bubbles.bed" if wildcards.region != 'all' else ""
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
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










# determine the variants that went into re-typing per category
rule collect_typed_variants:
	input:
		callset = "genotyping/leave1out/{callset}/preprocessed-vcfs/{sample}-{callset}_bi_no-missing.vcf.gz",
		regions= region_to_bed,
		ids="genotyping/leave1out/{callset}/untypable-ids/{sample}-untypable.tsv"
	output:
		"genotyping/leave1out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "bi|multi",
		vartype = "|".join(allowed_variants)
	resources:
		mem_total_mb=50000
	priority: 1
	shell:
		"zcat {input.callset} | python3 workflow/scripts/skip-untypable.py {input.ids} | python3 workflow/scripts/extract-varianttype.py {wildcards.vartype} | bedtools intersect -header -a - -b {input.regions} -u -f 0.5 | python3 workflow/scripts/get_ids.py > {output}"


# compute concordances
rule genotype_concordances:
	input:
		callset = "genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/{version}-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz",
		callset_tbi = "genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/{version}-{sample}_genotyping-biallelic-typable-{vartype}.vcf.gz.tbi",
		baseline = "genotyping/leave1out/{callset}/truth/truth-{sample}-{callset}-typable-{vartype}.vcf.gz",
		baseline_tbi = "genotyping/leave1out/{callset}/truth/truth-{sample}-{callset}-typable-{vartype}.vcf.gz.tbi",
		regions = region_to_bed,
		typed_ids = "genotyping/leave1out/{callset}/genotyped-ids/{sample}_{regions}_{vartype}.tsv"
	output:
		tmp_vcf1 = temp("genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/concordance/{regions}_{vartype}_base.vcf"),
		tmp_vcf2 = temp("genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/concordance/{regions}_{vartype}_call.vcf"),
		summary = "genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/concordance/{regions}_{vartype}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		sample = "|".join([s for s in reads_leave_one_out.keys()]),
		regions = "bi|multi",
		vartype = "|".join(allowed_variants)
	log:
		"genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/concordance/{regions}_{vartype}/summary.log"
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
		lambda wildcards: expand("genotyping/leave1out/{{callset}}/{{version}}/{sample}/{{coverage}}/{{metric}}/{{regions}}_{{vartype}}/summary.txt", sample = config['callsets'][wildcards.callset]['leave_one_out_samples'])
	output:
		"genotyping/leave1out/{callset}/{version}/plots/{coverage}/{metric}-{callset}-{version}-{coverage}_{regions}_{vartype}.tsv"
	params:
		samples = lambda wildcards: ','.join([c for c in config['callsets'][wildcards.callset]['leave_one_out_samples']]),
		outfile = "genotyping/leave1out/{callset}/{version}/plots/{coverage}/{metric}-{callset}-{version}-{coverage}_{regions}_{vartype}",
		folder = "genotyping/leave1out/{callset}/{version}"
	priority: 1
	shell:
		"python3 workflow/scripts/collect-results.py {wildcards.metric} {wildcards.coverage} {params.samples} {wildcards.regions} -variants {wildcards.vartype} -folder {params.folder} -outfile {params.outfile}"



# plot results of different subsampling runs
rule plotting_versions:
	input:
		lambda wildcards: expand("genotyping/leave1out/{{callset}}/{version}/plots/{{coverage}}/{m}-{{callset}}-{version}-{{coverage}}_{{regions}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', version=versions_leave_one_out, vartype=config['callsets'][wildcards.callset]['variants'])
	output:
		"genotyping/leave1out/{callset}/plots/comparison-versions/{metric}/{metric}_{coverage}_{regions}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped",
		regions="bi|multi"
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
		lambda wildcards: expand("genotyping/leave1out/{{callset}}/{version}/plots/{{coverage}}/concordance-{{callset}}-{version}-{{coverage}}_{{regions}}_{vartype}.tsv", version=versions_leave_one_out, vartype=config['callsets'][wildcards.callset]['variants'])
	output:
		"genotyping/leave1out/{callset}/plots/comparison-versions/concordance-vs-untyped/concordance-vs-untyped_{coverage}_{regions}.pdf"
	wildcard_constraints:
		regions="bi|multi"
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
		lambda wildcards: expand("genotyping/leave1out/{{callset}}/{{version}}/plots/{coverage}/{m}-{{callset}}-{{version}}-{coverage}_{{regions}}_{vartype}.tsv", m = wildcards.metric if wildcards.metric != 'untyped' else 'concordance', coverage=coverages_leave_one_out, vartype=config['callsets'][wildcards.callset]['variants'])
	output:
		"genotyping/leave1out/{callset}/plots/comparison-coverages/{metric}/{metric}_{version}_{regions}.pdf"
	wildcard_constraints:
		metric="concordance|precision-recall-typable|untyped",
		regions="bi|multi"
	priority: 1
	conda:
		"../envs/genotyping.yml"
	params:
		sources = lambda wildcards: ' '.join([wildcards.callset + '-' + wildcards.version + '-' + c + '_' + wildcards.regions for c in coverages_leave_one_out])
	shell:
		"python3 workflow/scripts/plot-results.py -files {input} -outname {output} -sources {params.sources} -metric {wildcards.metric}"

def collect_logs_all_unit_ids_bayestyper(wildcards):
	checkpoint_output = checkpoints.bayestyper_cluster.get(**wildcards).output[0]
	result = expand("genotyping/leave1out/{callset}/bayestyper/{sample}/{coverage}/genotype/bayestyper_unit_{unit_id}/bayestyper.log",
					callset=wildcards.callset,
					sample=wildcards.sample,
					coverage=wildcards.coverage,
					unit_id=glob_wildcards(os.path.join(checkpoint_output, "bayestyper_unit_{unit_id}/variant_clusters.bin")).unit_id)
	return sorted(result)

rule collect_runtime_info_bayestyper:
	input:
		bayestyper_kmers = "genotyping/leave1out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}_kmc.log",
		bayestyper_bloomfilter = "genotyping/leave1out/{callset}/bayestyper/{sample}/{coverage}/temp/kmers/{sample}_bloom.log",
		bayestyper_clusters = "genotyping/leave1out/{callset}/bayestyper/{sample}/{coverage}/clusters/clusters-log.log",
		bayestyper_genotype = collect_logs_all_unit_ids_bayestyper
	output:
		"genotyping/leave1out/{callset}/bayestyper/{sample}/{coverage}/bayestyper-{sample}.log"
	shell:
		"cat {input.bayestyper_kmers} {input.bayestyper_bloomfilter} {input.bayestyper_clusters} {input.bayestyper_genotype} > {output}"

rule collect_runtime_info_graphtyper:
	input:
		graphtyper_aligning_reads = "results/downsampling/{callset}/{coverage}/aligned/{sample}_full_mem.log",
		graphtyper_indexing_reads = lambda wildcards: expand("results/downsampling/{{callset}}/{{coverage}}/aligned/{{sample}}_full.chr{chrom}-index.log", chrom=ref_chromosomes[wildcards.callset]),
		graphtyper_genotype = lambda wildcards: expand("genotyping/leave1out/{{callset}}/graphtyper/{{sample}}/{{coverage}}/temp/{variant}/chr{chrom}.log", chrom=ref_chromosomes[wildcards.callset], variant=["snp-indel", "sv"])
	output:
		"genotyping/leave1out/{callset}/graphtyper/{sample}/{coverage}/graphtyper-{sample}.log"
	conda:
		"../envs/genotyping.yml"
	shell:
		"cat {input.graphtyper_aligning_reads} {input.graphtyper_indexing_reads} {input.graphtyper_genotype} > {output}"

rule collect_runtime_info_pangenie:
	input:
		pangenie_index = "genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/{version}-{sample}_index.log",
		pangenie_genotype = "genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/{version}-{sample}_genotyping.log"
	output:
		"genotyping/leave1out/{callset}/{sample}/{version}/{coverage}/{version}-{sample}.log"
	wildcard_constraints:
		version = "|".join([k for k in config['pangenie-modules'].keys()] + ['^' + k for k in config['pangenie']])
	shell:
		"cat {input.pangenie_index} {input.pangenie_genotype} > {output}"
        
        

# plot resources (single core CPU time and max RSS) for different subsampling runs
rule plotting_resources:
	input:
		lambda wildcards: expand("genotyping/leave1out/{{callset}}/{version}/{sample}/{{coverage}}/{version}-{sample}.log", version = versions_leave_one_out, sample = config['callsets'][wildcards.callset]['leave_one_out_samples'])
	output:
		"genotyping/leave1out/{callset}/plots/resources/resources_{callset}-{coverage}.pdf"
	conda:
		"../envs/genotyping.yml"
	params:
		outname = "genotyping/leave1out/{callset}/plots/resources/resources_{callset}-{coverage}",
		samples	= lambda wildcards: " ".join(config['callsets'][wildcards.callset]['leave_one_out_samples']),
		versions = " ".join(versions_leave_one_out)
	shell:
		"python3 workflow/scripts/plot-resources.py -files {input} -outname {params.outname} -samples {params.samples} -sizes {params.versions}"
