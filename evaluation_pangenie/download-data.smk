configfile: "config.json"

include: "rules/download-references.smk"
include: "rules/download-reads-NA24385.smk"
include: "rules/download-pangenome-graphs.smk"
include: "rules/download-giab-variants.smk"

sample = config['sample']
variants = config['variants']

rule all:
	input:
		# references
		'downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa',
		'downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai',
		
		# reads
		#expand("downloaded/reads/{sample}/raw_reads.fastq", sample=sample),
		
		# Input Pangenome Graph (PanGenie Github)
		'downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz',
		'downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze4_64haplotypes.vcf.gz',
		'downloaded/vcf/HPRC-GRCh38/Pangenome_graph_88haplotypes.vcf.gz',
		'downloaded/vcf/HPRC-CHM13/Pangenome_graph_88haplotypes.vcf.gz',
        #'downloaded/vcf/HGSVC-GRCh38/reduced_graph.vcf', 
        #'downloaded/vcf/HGSVC-GRCh38/graph_annotated_noNA24385.vcf.gz', 
        expand("downloaded/vcf/HGSVC-GRCh38/graph_annotated_{subset}_biallelic.vcf.gz", subset=['noNA24385', 'onlyNA24385']), 
        expand("downloaded/vcf/HGSVC-GRCh38/graph_annotated_{subset}_biallelic_{variant}.vcf.gz", subset=['onlyNA24385'], variant=variants), 
		
		# Callset (PanGenie Github)
		'downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz',
		'downloaded/vcf/HGSVC-GRCh38/Callset_freeze4_64haplotypes.vcf.gz',
		'downloaded/vcf/HPRC-GRCh38/Callset_88haplotypes.vcf.gz',
		'downloaded/vcf/HPRC-CHM13/Callset_88haplotypes.vcf.gz',
        "downloaded/vcf/HGSVC-GRCh38/graph_annotated.vcf.gz",
        #"downloaded/vcf/HGSVC-GRCh38/graph_annotated_biallelic.vcf.gz",
        #'downloaded/vcf/HGSVC-GRCh38/reduced_callset.vcf.gz', 
        #'downloaded/vcf/HGSVC-GRCh38/short_reduced_callset.vcf.gz', 

		# Benchmark set HG002_NA24385_son
		'downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz', 
		'downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi', 
		'downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.bed', 
		'downloaded/vcf/HGSVC-GRCh38/graph_biallelic.vcf.gz', 
		#'downloaded/vcf/HGSVC-GRCh38/benchmark_annotated_biallelic.vcf.gz', 
        expand('downloaded/vcf/giab/hg38/benchmark_annotated_biallelic_{variant}.vcf.gz', variant=variants), 
        
        # Create Alignments
        #"downloaded/BAM/aligned_reads.bam.bai",

        ## Other benchmarks
        #"downloaded/vcf/giab/hg38/HG002v11-align2-GRCh38.dip.vcf.gz", 
        #"downloaded/vcf/giab/hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
        #"downloaded/vcf/giab/hg38/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz",

        ## Download GRCh38 bundle for BayesTyper
        "downloaded/bayestyper/bayestyper_GRCh38_bundle.tar.gz",

        ## Annotate Callset and Pangenome Graph
        #"downloaded/vcf/HGSVC-GRCh38/reduced_graph_annotated.vcf",
        #"downloaded/vcf/HGSVC-GRCh38/reduced_graph_annotated_biallelic.vcf"
        #"downloaded/vcf/HGSVC-GRCh38/reduced_graph_annotated.vcf",
        #"downloaded/vcf/HGSVC-GRCh38/reduced_graph_annotated_biallelic.vcf"

        ## Get untypable ids 
        "downloaded/vcf/HGSVC-GRCh38/all-untypable-ids.tsv"