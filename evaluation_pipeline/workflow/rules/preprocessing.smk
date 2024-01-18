import gzip
truthsets_external = [t for c in config["callsets"].keys() for t in config["callsets"][c]["external"].keys()]
truthsets_leave1out = [t for c in config["callsets"].keys() for t in config["callsets"][c]["leave1out"].keys()]
callsets = [s for s in config['callsets'].keys()]

def query_left_one_out_samples(c, t):
	for line in gzip.open(config["callsets"][c]['bi'], 'rt'):
		if line.startswith("##"):
			continue
		if line.startswith("#"):
			all_samples = line.strip().split()[9:]
			left_one_out_samples = [k for k in all_samples if k != t]
		else:
			break
	return left_one_out_samples

# determine list of samples present in panel vcf, except the sample we are genotyping
assembly_samples = {}
for c in callsets:
	assembly_samples[c] = {}
	for p in ['leave1out', 'external']:
		assembly_samples[c][p] = {}
		truthsets = truthsets_leave1out if p == 'leave1out' else truthsets_external 
		for t in truthsets:
			assembly_samples[c][p][t] = query_left_one_out_samples(c, t)

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
		temp("preprocessing/{callset}/{sample}/input-panel/panel_{representation}_no-missing.vcf")
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000
	priority: 1
	wildcard_constraints:
		representation = "bi|multi"
	shell:
		"bcftools view {input} | python3 workflow/scripts/remove-missing.py {wildcards.sample} > {output}"

### We need to annotate the biallelic input panel to not forget ID variants in graphtyper postprocess
rule prepare_panel:
	input:
		"preprocessing/{callset}/{sample}/input-panel/panel_{representation}_no-missing.vcf"
	output:
		"preprocessing/{callset}/{sample}/input-panel/panel_{representation}.vcf"
	conda:
		"../envs/genotyping.yml"
	priority: 1
	wildcard_constraints:
		representation = "bi|multi"
	params: 
		lambda wildcards: "| python3 workflow/scripts/annotate-ids.py" if wildcards.representation == 'bi' else ""
	log:
		"logs/preprocessing/{callset}/{sample}/input-panel/panel_{representation}.log"
	resources:
		mem_total_mb=20000
	shell:
		"""
		(/usr/bin/time -v bcftools view --samples ^{wildcards.sample} {input} | bcftools view --min-ac 1 {params} > {output} ) &> {log}
		"""

rule normalize_input_panel:
	input:
		vcf = "preprocessing/{callset}/{sample}/input-panel/panel_multi.vcf",
		reference = lambda wildcards: config['callsets'][wildcards.callset]['reference']
	output:
		vcf = "preprocessing/{callset}/{sample}/input-panel/panel_multi_norm.vcf",
	shell:
		"""
		bcftools norm -m +any -d all -f {input.reference} {input.vcf} | bcftools sort > {output.vcf}
		"""

rule prepare_truth_leave1out_to_biallelic:
	input:
		truthset = "preprocessing/{callset}/{sample}/input-panel/panel_bi_no-missing.vcf",
		panel = lambda wildcards: config['callsets'][wildcards.callset]['bi'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai']
	output:
		"preprocessing/{callset}/{sample}/leave1out/truthset/truthset-biallelic.vcf.gz",
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000
	shell:
		"""
        bcftools view --samples {wildcards.sample} {input.truthset} | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | python workflow/scripts/annotate.py {input.panel} | bgzip -c > {output}
        tabix -p vcf {output}
        """

rule prepare_truth_external_to_biallelic:
	input:
		truthset = lambda wildcards: config['callsets'][wildcards.callset]['external'][wildcards.sample]["path"],
		panel = lambda wildcards: config['callsets'][wildcards.callset]['bi'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai']
	output:
		"preprocessing/{callset}/{sample}/external/truthset/truthset-biallelic.vcf.gz"
	shell:
		"""
		bcftools norm -m -any {input.truthset} | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | python workflow/scripts/annotate.py {input.panel} | bgzip -c > {output}
		tabix -p vcf {output}
		"""

rule split_truthset_by_variants:
	input:
		"preprocessing/{callset}/{sample}/{pipeline}/truthset/truthset-biallelic.vcf.gz",
	output:
		"preprocessing/{callset}/{sample}/{pipeline}/truthset/truthset_{vartype}-biallelic.vcf.gz",
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype="snp-indel|sv"
	resources:
		mem_total_mb=20000
	shell:
		"""
        bcftools view {input} | python workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output}
        tabix -p vcf {output}
        """


########################################################
###########      Utils for VCFEval         #############
########################################################

rule rtg_format:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['reference']
	output:
		directory("preprocessing/{callset}/SDF")
	resources:
		mem_total_mb=20000
	priority: 1
	shell:
		'{rtg} format -o {output} {input}'


########################################################
###########       Determine untypables     #############
########################################################

######## Compute Untypables ########
####  Find out which variants in the truth set are not contained in the pangenome graph (=untypables, because our methods are re-genotypers, i.e. no new variants detected)
#### Assign each variant a unique ID (if a variant matches with panel, use the same ID as panel).
rule determine_false_negatives:
	input:
		truthset="preprocessing/{callset}/{sample}/{pipeline}/truthset/truthset_{vartype}-biallelic.vcf.gz",
		panel = lambda wildcards: config['callsets'][wildcards.callset]['bi'],
		reference = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai'],
		sdf="preprocessing/{callset}/SDF"
	output:
		sample_vcf="preprocessing/{callset}/{sample}/{pipeline}/untypables/samples/{panelsample}-{vartype}.vcf.gz",
		fn="preprocessing/{callset}/{sample}/{pipeline}/untypables/samples/{panelsample}/{vartype}/fn.vcf.gz"
	wildcard_constraints:
		vartype="snp-indel|sv"
	conda:
		"../envs/genotyping.yml"
	params:
		tmp="preprocessing/{callset}/{sample}/{pipeline}/untypables/samples/{panelsample}/{vartype}_temp",
		outname="preprocessing/{callset}/{sample}/{pipeline}/untypables/samples/{panelsample}/{vartype}"
	threads: 1
	resources:
		mem_total_mb=50000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"logs/preprocessing/{callset}/{sample}/{pipeline}/untypables/samples/{panelsample}/{vartype}.log"
	shell:
		"""
		bcftools view --samples {wildcards.panelsample} {input.panel} | bcftools view --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		if [ "{wildcards.vartype}" == "snp-indel" ]; then 
			(/usr/bin/time -v {rtg} vcfeval -b {input.truthset} -c {output.sample_vcf} -t {input.sdf} -o {params.tmp} --squash-ploidy ) &> {log}
		elif [ "{wildcards.vartype}" == "sv" ]; then
			( /usr/bin/time -v {truvari} bench -b {input.truthset} -c {output.sample_vcf} -f {input.reference} -o {params.tmp} --pick multi -r 2000 --no-ref a -C 2000 --passonly) &> {log}
		fi
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""

# intersect the sets of FNs computed for each panel sample. The intersection then defines the set of unique/untypable variants
rule determine_unique:
	input:
		truthset = "preprocessing/{callset}/{sample}/{pipeline}/truthset/truthset_{vartype}-biallelic.vcf.gz",
		samples = lambda wildcards: expand("preprocessing/{{callset}}/{{sample}}/{{pipeline}}/untypables/samples/{sample}/{{vartype}}/fn.vcf.gz", sample=[s for s in assembly_samples[wildcards.callset][wildcards.pipeline][wildcards.sample] if not s == "CHM13"])
	output:
		unique_tsv="preprocessing/{callset}/{sample}/{pipeline}/untypables/untypables_{vartype}.tsv",
		unique_vcf="preprocessing/{callset}/{sample}/{pipeline}/untypables/untypables_{vartype}.vcf.gz"
	conda:
		"../envs/genotyping.yml"
	params:
		n_files = lambda wildcards: len([s for s in assembly_samples[wildcards.callset][wildcards.pipeline][wildcards.sample]  if not s == "CHM13"])
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=30
	shell:
		"""
		bcftools isec -n={params.n_files} -w1 {input.samples} | bgzip -c  > {output.unique_vcf}
		tabix -p vcf {output.unique_vcf}
		bcftools view -H {output.unique_vcf} | cut -f 3  > {output.unique_tsv}
		"""

rule merge_unique:
	input:
		tsv=expand("preprocessing/{{callset}}/{{sample}}/{{pipeline}}/untypables/untypables_{vartype}.tsv", vartype=["snp-indel", "sv"]),
		vcf=expand("preprocessing/{{callset}}/{{sample}}/{{pipeline}}/untypables/untypables_{vartype}.vcf.gz",  vartype=["snp-indel", "sv"])
	output:
		tsv="preprocessing/{callset}/{sample}/{pipeline}/untypables/untypables.tsv",
		vcf="preprocessing/{callset}/{sample}/{pipeline}/untypables/untypables.vcf"		
	conda:
		"../envs/genotyping.yml"
	params:
		n_files = lambda wildcards: len([s for s in assembly_samples[wildcards.callset][wildcards.pipeline][wildcards.sample]  if not s == "CHM13"])
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=30
	shell:
		"""
		cat {input.tsv} > {output.tsv}
		bcftools concat -a {input.vcf} > {output.vcf}
		"""

########################################################
###########    Prepare evaluation regions  #############
########################################################

rule alleles_per_bubble:
	input:
		"preprocessing/{callset}/{sample}/input-panel/panel_multi_norm.vcf.gz"
	output:
		plot = "preprocessing/{callset}/{sample}/regions/alleles-per-bubble.pdf",
		bed = "preprocessing/{callset}/{sample}/regions/multi-bubbles.bed"
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
		bed = "preprocessing/{callset}/{sample}/regions/multi-bubbles.bed",
		fai = lambda wildcards: config['callsets'][wildcards.callset]['reference'] + '.fai'
	output:
		bed = "preprocessing/{callset}/{sample}/regions/bi-bubbles.bed",
		tmp = temp("preprocessing/{callset}/{sample}/regions/bi-bubbles.fai")
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		sort -k1,1d -k 2,2n -k 3,3n {input.fai} > {output.tmp}
		bedtools complement -i {input.bed} -g {output.tmp} > {output.bed}
		"""