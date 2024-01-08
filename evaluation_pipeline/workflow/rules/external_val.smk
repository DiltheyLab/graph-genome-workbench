import gzip
rtg=config['programs']['rtg']
truvari = config['programs']['truvari']

# parameters
chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]

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

allowed_variants = ['snp-indel', 'sv']
coverages = ['full'] + [c for c in config['downsampling']]
versions_external = [v for v  in config['pangenie'].keys()] + [v for v in config['pangenie-modules'].keys()]

truthsets_sv = [t for c in config["callsets"].keys() for t in config["callsets"][c]["truthsets"].keys() if 'sv' in config["callsets"][c]["truthsets"][t]["vartype"]]
truthsets_small = [t for c in config["callsets"].keys() for t in config["callsets"][c]["truthsets"].keys() if "snp-indel" in config["callsets"][c]["truthsets"][t]["vartype"]]

#truthsets = truthsets_sv + truthsets_small

truthsets = [t for c in config["callsets"].keys() for t in config["callsets"][c]["truthsets"].keys()]
callsets = [s for s in config['callsets'].keys()]
## callsets & regions to be determined during execution

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
	for t in truthsets:
		assembly_samples[c][t] = query_left_one_out_samples(c, t)


##########################################################################
##  Preprocessing (Callsets: Biallelic and Multiallelic) and Truthset  ###
##########################################################################

### Preprocess Multiallelic (Pangenome)
rule prepare_panel:
	input:
		vcf = lambda wildcards: config['callsets'][wildcards.callset][wildcards.representation]
	output:
		vcf = "external-eval/{callset}/{sample}/input-panel/panel_{representation}.vcf"
	wildcard_constraints:
		representation = "multi|bi"
	shell:
		"""
		bcftools view --samples ^{wildcards.sample} {input.vcf} | bcftools view --min-ac 1 > {output.vcf}
		"""


########################################################
##################    run PanGenie    ##################
########################################################

# run pangenie
## If in input reads = lambda wildcards: "results/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
## the reads will be aligned before, but that's what we need, since graphtyper is a mapped-based approach (i.e. we need the aligned reads anyways if we want to apply Graphtyper). Just for Pangenie, it's not necessary.
rule pangenie:
	input:
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "data/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		vcf = "external-eval/{callset}/{sample}/input-panel/panel_multi.vcf"
	output:
		genotyping = "external-eval/{callset}/{sample}/{version}/{coverage}/genotyping.vcf"
	log:
		"external-eval/{callset}/{sample}/{version}/{coverage}/{version}-{sample}.log"
	threads: 24
	resources:
		mem_total_mb=80000,
		runtime_hrs=5,
		runtime_min=1
	priority: 1
	params:
		tmpdir="external-eval/{callset}/{sample}/{version}/{coverage}/temp",
		out_prefix="external-eval/{callset}/{sample}/{version}/{coverage}/temp/{sample}-{version}",
		pangenie = lambda wildcards: config['pangenie'][wildcards.version]
	wildcard_constraints:
		version = "|".join([k for k in config['pangenie'].keys()] + ['^' + k for k in config['pangenie-modules']])
	shell:
		"""
        mkdir -p {params.tmpdir}
		(/usr/bin/time -v {params.pangenie} -i <(zcat {input.reads}) -v {input.vcf} -r {input.fasta} -o {params.out_prefix} -s {wildcards.sample} -j {threads} -t {threads} -g ) &> {log}
		mv {params.out_prefix}_genotyping.vcf {output.genotyping}
		rm -r {params.tmpdir}
        """

# run pangenie in the modularized way (> v2.1.1)
rule pangenie_modules:
	input:
		reads = lambda wildcards: reads_leave_one_out[wildcards.sample] if wildcards.coverage == 'full' else "results/downsampling/{callset}/{coverage}/{sample}_{coverage}.fa.gz",
		fasta = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		vcf = "external-eval/{callset}/{sample}/input-panel/panel_multi.vcf"
	output:
		genotyping = "external-eval/{callset}/{sample}/{version}/{coverage}/genotyping.vcf"
	log:
		index = "external-eval/{callset}/{sample}/{version}/{coverage}/index.log",
		genotype = "external-eval/{callset}/{sample}/{version}/{coverage}/genotyping.log"
	threads: 24
	resources:
		mem_total_mb=50000,
		runtime_hrs=5,
		runtime_min=1
	priority: 1
	params:
		tmpdir="external-eval/{callset}/{sample}/{version}/{coverage}/temp",
		out_prefix="external-eval/{callset}/{sample}/{version}/{coverage}/temp/{sample}-{version}",
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

# convert genotyped VCF to biallelic representation
rule convert_genotypes_to_biallelic:
	input:
		genotyped_vcf = "external-eval/{callset}/{sample}/{version}/{coverage}/genotyping.vcf",
		biallelic = lambda wildcards: config['callsets'][wildcards.callset]['bi']
	output:
		"external-eval/{callset}/{sample}/{version}/{coverage}/genotyping-biallelic.vcf"
	log:
		"external-eval/{callset}/{sample}/{version}/{coverage}/genotyping-biallelic.log"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=30000
	priority: 1
	shell:
		"(cat {input.genotyped_vcf} | python workflow/scripts/convert-to-biallelic.py {input.biallelic} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1 -k2,2n \"}}' > {output}) &> {log}"


####################################################################################################
######   Split Callset and Truthset benchmark depending on Variant Type (snp-indels vs. SVs)  ######  
####################################################################################################

# split between snp-indels and SVs
rule split_genotyped_by_variantsize_callset:
    input:
        genotyped_callset = "external-eval/{callset}/{sample}/{version}/{coverage}/genotyping-biallelic.vcf"
    output:
        "external-eval/{callset}/{sample}/{version}/{coverage}/genotyping-biallelic_{vartype}.vcf.gz"
    wildcard_constraints:
        vartype="snp-indel|sv"
    shell:
        """
        cat {input.genotyped_callset} | python workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output}
        tabix -p vcf {output}
        """

# split between snp-indels and SVs
rule convert_to_biallelic_and_split_by_variantsize_truthset:
    input:
        groundtruth = lambda wildcards: config['callsets'][wildcards.callset]['truthsets'][wildcards.sample]["path"]
    output:
        "external-eval/{callset}/vcf/truthset-{sample}_{vartype}.vcf.gz"
    wildcard_constraints:
        vartype="snp-indel|sv"
    shell:
        """
        bcftools norm -m -any {input.groundtruth} | python workflow/scripts/extract-varianttype.py {wildcards.vartype} | bgzip -c > {output}
        tabix -p vcf {output}
        """


####################################################################################################
#  find out which variants in the truth set are not contained in the pangenome graph (=untypables)
####################################################################################################

# assign each variant a unique ID (if a variant matches with panel, use the same ID as panel).
rule annotate_variants_from_callset_in_truthset:
	input:
		truthset = "external-eval/{callset}/vcf/truthset-{sample}_{vartype}.vcf.gz",
		panel = lambda wildcards: config['callsets'][wildcards.callset]['bi'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai']
	output:
		vcf="external-eval/{callset}/vcf/truthset-{sample}_{vartype}-annotated.vcf.gz",
		tbi="external-eval/{callset}/vcf/truthset-{sample}_{vartype}-annotated.vcf.gz.tbi"
	wildcard_constraints:
		callsets = "|".join(callsets),
		sample = "|".join(truthsets), 
		vartype="snp-indel|sv"
	resources:
		mem_total_mb=20000,
		runtime_hrs=1,
		runtime_min=59
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view {input.truthset} | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | python workflow/scripts/annotate.py {input.panel} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf} 
		"""


# format reference (needed for running vcfeval)
rule rtg_format_callsets:
	input:
		reference = lambda wildcards: config['callsets'][wildcards.callset]['reference']
	output:
		directory("external-eval/{callset}/SDF")
	resources:
		mem_total_mb=20000
	conda:
		"../envs/genotyping.yml"
	shell:
		"{rtg} format -o {output} {input}"


# for each panel sample determine its false negatives comparing to the truth set.
# these are variants only in the truth, but not detected in the sample itself.
# later, intersect all false negatives across the panel samples, to find out
# which variants are not present in the pangenome graph, but are in the truth.
# these variants are not accessible by re-genotyping methods, that are unable
# to detect variants themselves ("untypables").
rule determine_false_negatives_vcfeval:
	input:
		truthset="external-eval/{callset}/vcf/truthset-{sample}_snp-indel-annotated.vcf.gz",
		panel = lambda wildcards: config['callsets'][wildcards.callset]['bi'],
		reference = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai'],
		sdf="external-eval/{callset}/SDF"
	output:
		sample_vcf="external-eval/{callset}/vcf/{sample}/samples/{panelsample}-vcfeval.vcf.gz",
		fn="external-eval/{callset}/vcf/{sample}/samples/{panelsample}/vcfeval/fn.vcf.gz"
	wildcard_constraints:
		truthset="|".join(truthsets_small), 
        callsets = "|".join(callsets),
		vartype="snp-indel"
	conda:
		"../envs/genotyping.yml"
	params:
		tmp="external-eval/{callset}/vcf/{sample}/samples/{panelsample}/vcfeval_temp",
		outname="external-eval/{callset}/vcf/{sample}/samples/{panelsample}/vcfeval"
	threads: 28
	resources:
		mem_total_mb=50000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"external-eval/{callset}/vcf/{sample}/samples/{panelsample}/vcfeval.log"
	shell:
		"""
		bcftools view --samples {wildcards.panelsample} {input.panel} | bcftools view --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		(/usr/bin/time -v {rtg} vcfeval -b {input.truthset} -c {output.sample_vcf} -t {input.sdf} -o {params.tmp} --squash-ploidy ) &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""


# same as previous rule, but for SVs using truvari instead of vcfeval
rule determine_false_negatives_truvari:
	input:
		truthset = "external-eval/{callset}/vcf/truthset-{sample}_sv-annotated.vcf.gz",
		panel = lambda wildcards: config['callsets'][wildcards.callset]['bi'],
		reference = lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai'],
		sdf = "external-eval/{callset}/SDF"
	output:
		sample_vcf = "external-eval/{callset}/vcf/{sample}/samples/{panelsample}-truvari.vcf.gz",
		fn = "external-eval/{callset}/vcf/{sample}/samples/{panelsample}/truvari/fn.vcf.gz"
	wildcard_constraints:
		truthset="|".join(truthsets_sv), 
        callsets = "|".join(callsets),
		vartype="sv"
	conda:
		"../envs/genotyping.yml"
	params:
		tmp="external-eval/{callset}/vcf/{sample}/samples/{panelsample}/truvari_temp",
		outname="external-eval/{callset}/vcf/{sample}/samples/{panelsample}/truvari"
	threads: 1
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	log:
		"external-eval/{callset}/vcf/{sample}/samples/{panelsample}/truvari.log"
	shell:
		"""
		bcftools view --samples {wildcards.panelsample} {input.panel} | bcftools view --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.sample_vcf}
		tabix -p vcf {output.sample_vcf}
		( /usr/bin/time -v {truvari} bench -b {input.truthset} -c {output.sample_vcf} -f {input.reference} -o {params.tmp} --pick multi -r 2000 --no-ref a -C 2000 --passonly) &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""

def get_truthset_given_method(wildcards):
    if wildcards.method == 'vcfeval':
        return f"external-eval/{wildcards.callset}/vcf/truthset-{wildcards.sample}_snp-indel-annotated.vcf.gz"
    if wildcards.method == 'truvari':
        return f"external-eval/{wildcards.callset}/vcf/truthset-{wildcards.sample}_sv-annotated.vcf.gz"

# intersect the sets of FNs computed for each panel sample. The intersection then defines the set of unique/untypable variants
# Take into account that vartype=snp-indel relates to method=vcfeval and vartype=sv relates to method=truvari
rule determine_unique:
	input:
		truthset = get_truthset_given_method,
		samples = lambda wildcards: expand("external-eval/{{callset}}/vcf/{{sample}}/samples/{sample}/{{method}}/fn.vcf.gz", sample=[s for s in assembly_samples[wildcards.callset][wildcards.sample] if not s == "CHM13"])
	output:
		unique_tsv="external-eval/{callset}/vcf/{sample}-unique_{method}.tsv",
		unique_vcf="external-eval/{callset}/vcf/{sample}-unique_{method}.vcf"
	conda:
		"../envs/genotyping.yml"
	params:
		n_files = lambda wildcards: len([s for s in assembly_samples[wildcards.callset][wildcards.sample]  if not s == "CHM13"])
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=30
	shell:
		"""
		bcftools isec -n={params.n_files} -w1 {input.samples}  > {output.unique_vcf}
		grep -v '#' {output.unique_vcf} | cut -f 3  > {output.unique_tsv}
		"""


rule count_unique_type:
	input:
		total="external-eval/{callset}/vcf/truthset-{sample}.vcf.gz",
		unique="external-eval/{callset}/vcf/{sample}-unique_{method}.vcf"
	output:
		total="external-eval/{callset}/vcf/{sample}-total_{method}-{vartype}.vcf",
		unique="external-eval/{callset}/vcf/{sample}-unique_{method}-{vartype}.vcf"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		vartype="snp-indel|sv"
	shell:
		"""
		bcftools view {input.total} | python workflow/scripts/extract-varianttype.py {wildcards.vartype} > {output.total}
		bcftools view {input.unique} | python workflow/scripts/extract-varianttype.py {wildcards.vartype} > {output.unique}
		"""

####################################################################################################
#  compare genotyping results to the truthset, excluding untypable variants which cannot be 
#  genotyped correctly by a re-genotyper
####################################################################################################

# remove untypables from groundtruth and extract only variants of a specific type.
rule remove_untypables_if_apply_and_extract_variant_type_groundtruth:
	input:
		vcf = get_truthset_given_method,
		untypable = "external-eval/{callset}/vcf/{sample}-unique_{method}.tsv"
	output:
		vcf="external-eval/{callset}/{sample}/groundtruth-{filter}_{method}.vcf.gz",
		tbi="external-eval/{callset}/{sample}/groundtruth-{filter}_{method}.vcf.gz.tbi"
	wildcard_constraints:
		sample = "|".join(truthsets),
		callset= "|".join(callsets),
		filter="typable|all"
	resources:
		mem_total_mb=30000
	conda:
		"../envs/genotyping.yml"
	params:
		untypables = lambda wildcards: f" python workflow/scripts/skip-untypable.py external-eval/{wildcards.callset}/vcf/{wildcards.sample}-unique_{wildcards.method}.tsv |" if wildcards.filter == 'typable' else "",
		vartype = lambda wildcards: "snp-indel" if wildcards.method == 'vcfeval' else "sv"
	shell:
		"""
		bcftools norm -m -any {input.vcf} | {params.untypables} python workflow/scripts/extract-varianttype.py {params.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""
        
def get_callset_given_method(wildcards):
    if wildcards.method == 'vcfeval':
        return f"external-eval/{wildcards.callset}/{wildcards.sample}/{wildcards.version}/{wildcards.coverage}/genotyping-biallelic_snp-indel.vcf.gz"
    if wildcards.method == 'truvari':
        return f"external-eval/{wildcards.callset}/{wildcards.sample}/{wildcards.version}/{wildcards.coverage}/genotyping-biallelic_sv.vcf.gz"


# remove untypables from genotyped callset and extract only variants of a specific type.
rule remove_untypables_if_apply_and_extract_variant_type_genotyped_callset:
	input:
		vcf= get_callset_given_method,
		untypable = "external-eval/{callset}/vcf/{sample}-unique_{method}.tsv"
	output:
		vcf="external-eval/{callset}/{sample}/{version}/{coverage}/callset-{filter}_{method}.vcf.gz",
		tbi="external-eval/{callset}/{sample}/{version}/{coverage}/callset-{filter}_{method}.vcf.gz.tbi"
	wildcard_constraints:
		callset= "|".join(callsets),
		filter="typable|all",
		version="|".join(versions_external),
		coverage="|".join(coverages)
	resources:
		mem_total_mb=30000
	conda:
		"../envs/genotyping.yml"
	params:
		untypables = lambda wildcards: f"| python workflow/scripts/skip-untypable.py external-eval/{wildcards.callset}/vcf/{wildcards.sample}-unique_{wildcards.method}.tsv |" if wildcards.filter == "typable" else "",
		vartype = lambda wildcards: "snp-indel" if wildcards.method == 'vcfeval' else "sv"
	shell:
		"""
		bcftools norm -m -any {input.vcf} {params.untypables} python workflow/scripts/extract-varianttype.py {params.vartype} | bgzip -c > {output.vcf}
		tabix -p vcf {output.vcf}
		"""

####################################################################################################
# prepare beds for biallelic and complex graph regions
# AND prepare regions (stratification BED files intersected with callable regions as defined by the BED 
# files that come with each ground truth)
####################################################################################################

rule alleles_per_bubble:
	input:
		"external-eval/{callset}/{sample}/input-panel/panel_multi.vcf"
	output:
		plot = "external-eval/{callset}/{sample}/alleles-per-bubble.pdf",
		bed = "external-eval/{callset}/{sample}/complex-bubbles.bed"
	conda:
		"../envs/genotyping.yml"
	resources:
		mem_total_mb=20000,
		runtime_hrs=1
	shell:
		"cat {input} | python workflow/scripts/variant-statistics.py {output.plot} 1 > {output.bed}"


rule prepare_beds:
	input:
		bed = "external-eval/{callset}/{sample}/complex-bubbles.bed",
		fai = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai']
	output:
		bed = "external-eval/{callset}/{sample}/biallelic-bubbles.bed",
		tmp = temp("external-eval/{callset}/{sample}/biallelic-bubbles.fai")
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		sort -k1,1d -k 2,2n -k 3,3n {input.fai} > {output.tmp}
		bedtools complement -i {input.bed} -g {output.tmp} > {output.bed}
		"""


def define_bed_given_region(wildcards):
	if wildcards.region == 'all':
		return config["callsets"][wildcards.callset]["truthsets"][wildcards.sample]["callable_regions"]
	elif wildcards.region == 'multi':
		return f"external-eval/{wildcards.callset}/{wildcards.sample}/complex-bubbles.bed"
	elif wildcards.region == 'bi':
		return f"external-eval/{wildcards.callset}/{wildcards.sample}/biallelic-bubbles.bed"
	else:
		print("Region not in [all, multi, bi] ")


rule prepare_evaluation_beds:
	input:
		callable_regions = lambda wildcards: config['callsets'][wildcards.callset]["truthsets"][wildcards.sample]["callable_regions"],
		bed = define_bed_given_region
	output:
		"external-eval/{callset}/{sample}/bed-files/{sample}_{region}.bed"
	wildcard_constraints:
		callset = "|".join(callsets),
		sample = "|".join(truthsets),
		region = "all|bi|multi"
	resources:
		mem_total_mb=30000,
		runtime_hrs=1
	conda:
		"../envs/genotyping.yml"
	shell:
		"bedtools intersect -a {input.callable_regions} -b {input.bed} > {output}"

####################################################################################################
###  Compute Precision/Recall for both snp-indel vs. SVs
####################################################################################################

# compute precision/recall for small variants
rule vcfeval_callsets:
	input:
		callset="external-eval/{callset}/{sample}/{version}/{coverage}/callset-{filter}_vcfeval.vcf.gz",
		callset_tbi="external-eval/{callset}/{sample}/{version}/{coverage}/callset-{filter}_vcfeval.vcf.gz.tbi",
		baseline="external-eval/{callset}/{sample}/groundtruth-{filter}_vcfeval.vcf.gz",
		baseline_tbi="external-eval/{callset}/{sample}/groundtruth-{filter}_vcfeval.vcf.gz.tbi",
		regions= "external-eval/{callset}/{sample}/bed-files/{sample}_{region}.bed",
		reference=lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai'],
		sdf="external-eval/{callset}/SDF"
	output:
		fixed_vcf="external-eval/{callset}/{sample}/vcfeval_{version}_{coverage}_{filter}_region-{region}_fixed.vcf.gz",
		fixed_vcf_tbi="external-eval/{callset}/{sample}/vcfeval_{version}_{coverage}_{filter}_region-{region}_fixed.vcf.gz.tbi",
		summary="external-eval/{callset}/{sample}/vcfeval_{version}_{coverage}_{filter}_region-{region}/summary.txt"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		filter = "all|typable"
	log:
		"external-eval/{callset}/{sample}/vcfeval_{version}_{coverage}_{filter}_region-{region}.log"
	params:
		tmp = "external-eval/{callset}/{sample}/vcfeval_{version}_{coverage}_{filter}_region-{region}_tmp",
		outname = "external-eval/{callset}/{sample}/vcfeval_{version}_{coverage}_{filter}_region-{region}",
		bed = lambda wildcards: f"--evaluation-regions external-eval/{wildcards.callset}/{wildcards.sample}/bed-files/{wildcards.sample}_{wildcards.region}.bed" if wildcards.region != 'all' else ""
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
	threads: 24
	shell:
		"""
		bcftools view {input.callset} --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.fixed_vcf}
		tabix -p vcf {output.fixed_vcf}
		{rtg} vcfeval -b {input.baseline} -c {output.fixed_vcf} -t {input.sdf} -o {params.tmp} --evaluation-regions {input.regions} --threads {threads} > {output.summary}.tmp
		mv {params.tmp}/* {params.outname}/
		mv {output.summary}.tmp {output.summary}
		rm -r {params.tmp}
		"""

# compute precision/recall for SVs
rule truvari_callsets:
	input:
		callset="external-eval/{callset}/{sample}/{version}/{coverage}/callset-{filter}_truvari.vcf.gz",
		callset_tbi="external-eval/{callset}/{sample}/{version}/{coverage}/callset-{filter}_truvari.vcf.gz.tbi",
		baseline="external-eval/{callset}/{sample}/groundtruth-{filter}_truvari.vcf.gz",
		baseline_tbi="external-eval/{callset}/{sample}/groundtruth-{filter}_truvari.vcf.gz.tbi",
		regions= "external-eval/{callset}/{sample}/bed-files/{sample}_{region}.bed",
		reference=lambda wildcards: config['callsets'][wildcards.callset]['reference'],
		ref_index = lambda wildcards: config['callsets'][wildcards.callset]['reference_fai']
	output:
		fixed_vcf="external-eval/{callset}/{sample}/truvari_{version}_{coverage}_{filter}_region-{region}_fixed.vcf.gz",
		fixed_vcf_tbi="external-eval/{callset}/{sample}/truvari_{version}_{coverage}_{filter}_region-{region}_fixed.vcf.gz.tbi",
		summary="external-eval/{callset}/{sample}/truvari_{version}_{coverage}_{filter}_region-{region}/summary.json"
	conda:
		"../envs/genotyping.yml"
	wildcard_constraints:
		filter = "all|typable"
	log:
		"external-eval/{callset}/{sample}/truvari_{version}_{coverage}_{filter}_region-{region}.log"
	params:
		tmp = "external-eval/{callset}/{sample}/truvari_{version}_{coverage}_{filter}_region-{region}_tmp",
		outname = "external-eval/{callset}/{sample}/truvari_{version}_{coverage}_{filter}_region-{region}",
		bed = lambda wildcards: f"--includebed external-eval/{wildcards.callset}/{wildcards.sample}/bed-files/{wildcards.sample}_{wildcards.region}.bed" if wildcards.region != 'all' else ""
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=40
	threads: 1
	shell:
		"""
		bcftools view {input.callset} --min-ac 1 | python workflow/scripts/prepare-for-vcfeval.py {input.ref_index} | bgzip -c > {output.fixed_vcf}
		tabix -p vcf {output.fixed_vcf}
		{truvari} bench -b {input.baseline} -c {output.fixed_vcf} -f {input.reference} -o {params.tmp} --pick multi {params.bed} -r 2000 --no-ref a -C 2000 --passonly &> {log}
		mv {params.tmp}/* {params.outname}/
		rm -r {params.tmp}
		"""
