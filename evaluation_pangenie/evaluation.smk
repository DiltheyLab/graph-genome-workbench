configfile: "config.json"

rtg = config['rtg']
sample = config['sample']
variants = config['variants']
tools = config['tools']
eval_methods = config['evaluation_methods']

rule all:
    input:
        expand("evaluation/{sample}/{variant}/{tool}/{eval_method}/summary.txt", sample=sample, variant=variants, tool=tools, eval_method=eval_methods),
        #expand("downloaded/vcf/HGSVC-GRCh38/graph_annotated_onlyNA24385_biallelic_typables_{variant}.vcf.gz", variant=variants),
        #expand("genotyping/{sample}/genotyped-ids/{tool}_{variant}.tsv", sample=sample, variant=variants, tool=tools),
        #expand('evaluation/{sample}/{variant}/{tool}/vcfeval/fn.vcf.gz', sample=sample, variant=variants, tool=tools),
        #config['reference_sdf_dir'],
            

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



# prepare reference for vcfeval
rule rtg_format_reference_SDF:
	input:
		config['reference']
	output:
		directory(config['reference_sdf_dir'])
	log:
		'evaluation/prepare_reference.log'
	shell:
		#rm -r downloaded/SDF
		"""
		/usr/bin/time -v {rtg} format -o {output} {input} &> {log}
		"""

####### Remove untypable IDs
rule remove_untypable_truth:
    input:
        vcf="downloaded/vcf/HGSVC-GRCh38/graph_annotated_onlyNA24385_biallelic_{variant}.vcf.gz", 
        ids=config['untypable']
    output:
        vcf="downloaded/vcf/HGSVC-GRCh38/graph_annotated_onlyNA24385_biallelic_typables_{variant}.vcf.gz"
    shell:
        """
        zcat {input.vcf} | python3 scripts/skip-untypable.py {input.ids} | bgzip > {output.vcf}
        tabix -p vcf {output}
        """

rule remove_untypable_callset:
    input:
        vcf='genotyping/{sample}/{tool}/genotyping_{variant}.vcf.gz',
        ids=config['untypable']
    output:
        vcf='genotyping/{sample}/{tool}/genotyping_typables_{variant}.vcf.gz'
    shell:
        """
        zcat {input.vcf} | python3 scripts/skip-untypable.py {input.ids} | bgzip > {output.vcf}
        tabix -p vcf {output}
        """

################### assess performance: precision and recall
rule assess_precision_recall_per_tool_and_variantsize:
	input:
		truthvcf = "downloaded/vcf/HGSVC-GRCh38/graph_annotated_onlyNA24385_biallelic_typables_{variant}.vcf.gz",
		genotyped_sample_vcf = 'genotyping/{sample}/{tool}/genotyping_typables_{variant}.vcf.gz',
		reference_sdf_dir = config['reference_sdf_dir'], 
		#truthvcf = 'downloaded/vcf/giab/hg38/benchmark__biallelic_{variant}.vcf.gz',
	output:
		'evaluation/{sample}/{variant}/{tool}/vcfeval/summary.txt'
	params:
 		temp = 'evaluation/{sample}/{variant}/{tool}/vcfeval_temp',
 		outname = 'evaluation/{sample}/{variant}/{tool}/vcfeval'
	log:
		'evaluation/{sample}/{variant}/{tool}/vcfeval.log'
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	shell:
		"""
		/usr/bin/time -v {rtg} vcfeval -b {input.truthvcf} -c {input.genotyped_sample_vcf} -t {input.reference_sdf_dir} -o {params.temp} --squash-ploidy --ref-overlap --all-records  --Xmax-length 30000 &> {log}
		mv {params.temp}/* {params.outname}/
		rm -r {params.temp}
        """

################### assess performance: weighted Genotype Concordance
# determine the variants that went into re-typing per category
rule collected_typed_variants:
    input:
        genotyped_sample_vcf = 'genotyping/{sample}/{tool}/genotyping_typable_{variant}.vcf.gz'
    output:
        "genotyping/{sample}/genotyped-ids/{tool}_{variant}.tsv"
    shell:
        "zcat {input.genotyped_sample_vcf} | python3 scripts/get_ids.py > {output}"


# compute concordances
rule genotype_concordances:
	input:
		truthvcf = "downloaded/vcf/HGSVC-GRCh38/graph_annotated_onlyNA24385_biallelic_typables_{variant}.vcf.gz",
		genotyped_sample_vcf = 'genotyping/{sample}/{tool}/genotyping_typables_{variant}.vcf.gz',
		typed_ids = "genotyping/{sample}/genotyped-ids/{tool}_{variant}.tsv"
	output:
		summary="evaluation/{sample}/{variant}/{tool}/wGC/summary.txt"
	log:
		"evaluation/{sample}/{variant}/{tool}/genotype_concordances.log"
	resources:
		mem_total_mb=20000,
		runtime_hrs=0,
		runtime_min=20
	shell:
		"""
		python3 scripts/genotype-evaluation.py {input.truthvcf} {input.genotyped_sample_vcf} {input.typed_ids} 2> {log} 1> {output.summary}
		"""








