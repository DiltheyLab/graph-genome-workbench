configfile: "config.json"

rtg = config['rtg']
sample = config['sample']

rule all:
	input:
		expand('evaluation/{sample}/vcfeval/fn.vcf.gz', sample=sample)
        #config['reference_sdf_dir']
        

# prepare reference for vcfeval
rule rtg_format_reference_SDF:
	input:
		config['reference']
	output:
		directory(config['reference_sdf_dir'])
	log:
		'evaluation/prepare_reference.log'
	shell:
		"""
		/usr/bin/time -v {rtg} format -o {output} {input} &> {log}
		"""

# assess performance
rule assess_performance:
	input:
		truthvcf = config['truthvcf'],
		genotyped_sample_vcf = config['genotyped_sample_vcf'],
		reference_sdf_dir = config['reference_sdf_dir'] 
	output:
		'evaluation/{sample}/vcfeval/fn.vcf.gz'
	params:
 		temp = 'evaluation/{sample}/vcfeval_temp',
 		outname = 'evaluation/{sample}/vcfeval'
	log:
		'evaluation/{sample}_vcfeval.log'
	resources:
		mem_total_mb=30000,
		runtime_hrs=0,
		runtime_min=59
	shell:
		"""
		/usr/bin/time -v {rtg} vcfeval -b {input.truthvcf} -c {input.genotyped_sample_vcf} -t {input.reference_sdf_dir} -o {params.temp} --squash-ploidy --ref-overlap --all-records > {log} --Xmax-length 30000 &> {log}
		mv {params.temp}/* {params.outname}/
		rm -r {params.temp}
        """
