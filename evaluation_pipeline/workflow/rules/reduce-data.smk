##### We need to reduce the size of: Pangenome, biallelelic callset (two samples NA12878 and NA24285), truthset(HG002), reads 
##### as long as they contain of the cases: known & unknown alleles, snp-indel, sv (different types)
#configfile: "config/reduced.yaml"

rule all:
	input:
		config['reads_reduced'],
		config['multi_reduced']

rule reduce_reads:
	input:
		config['reads']
	output:
		config['reads_reduced']
	shell:
		"""
        zcat {input} | head -40000 | bgzip -c > {output}
		"""

rule reduce_pangenome:
	input:
		config['multi']
	output:
		config['multi_reduced']
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view {input} | head -120 | bgzip -c > {output}
		"""

rule reduce_pangenome:
	input:
		config['multi']
	output:
		config['multi_reduced']
	conda:
		"../envs/genotyping.yml"
	shell:
		"""
		bcftools view {input} | head -120 | bgzip -c > {output}
		"""

mkdir -p reduced-data
samtools faidx data/downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa chr17 > reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa > reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
zcat data/downloaded/reads/NA12878/NA12878_raw.fastq.gz | head -4000 | bgzip -c > reduced-data/NA12878_raw.fastq.gz
zcat data/downloaded/reads/NA24385/NA24385_raw.fastq.gz | head -4000 | bgzip -c > reduced-data/NA24385_raw.fastq.gz
bcftools view -r chr17 data/downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz | bgzip -c > reduced-data/Pangenome_graph_freeze3_64haplotypes.vcf.gz
bcftools view -r chr17 data/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz | bgzip -c > reduced-data/Callset_freeze3_64haplotypes.vcf.gz
bcftools view -r chr17 data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | python workflow/scripts/prepare-for-vcfeval.py reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | bgzip -c > reduced-data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

### Copy 
mkdir -p reduced-data
zcat data/downloaded/reads/NA12878/NA12878_raw.fastq.gz | head -4000 | bgzip -c > reduced-data/NA12878_raw.fastq.gz
zcat data/downloaded/reads/NA24385/NA24385_raw.fastq.gz | head -4000 | bgzip -c > reduced-data/NA24385_raw.fastq.gz
bcftools view -r chr17 data/downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz | bgzip -c > reduced-data/Pangenome_graph_freeze3_64haplotypes.vcf.gz
bcftools view -r chr17 data/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz | bgzip -c > reduced-data/Callset_freeze3_64haplotypes.vcf.gz
bcftools view -r chr17 data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | bgzip -c > reduced-data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
samtools faidx data/downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa chr17 > reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa > reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai





### Trash
bcftools view external-eval/HGSVC/pangenie.v3/NA12878/full/genotyping.vcf | bgzip -c > genotyping.vcf.gz
tabix -p vcf genotyping.vcf.gz
bcftools view -r chr17:65651680-65651681 genotyping.vcf.gz | bgzip -c > genotyping_special.vcf.gz
bcftools view genotyping_special.vcf.gz | python workflow/scripts/convert-to-biallelic.py reduced-data/Callset_freeze3_64haplotypes.vcf.gz



bcftools view data/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz | grep 65651681
bcftools view reduced-data/Callset_freeze3_64haplotypes.vcf.gz | grep 65651681
bcftools view reduced-data/Pangenome_graph_freeze3_64haplotypes.vcf.gz | grep 65651681
bcftools view external-eval/HGSVC/pangenie.v3/NA12878/full/genotyping.vcf | grep 65651681
bcftools view external-eval/HGSVC/pangenie.v3/NA12878/full/genotyping-biallelic.vcf | grep 65651681
bcftools norm -m +any -d all -f reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa reduced-data/Pangenome_graph_freeze3_64haplotypes.vcf.gz | bcftools sort | grep 65651681


bcftools view -H external-eval/HGSVC/NA12878/groundtruth-typable_truvari.vcf.gz | wc -l
bcftools view -H external-eval/HGSVC/NA12878/pangenie.v3/full/callset-typable_truvari.vcf.gz | wc -l 

bcftools view -H external-eval/HGSVC/NA12878/truvari_pangenie.v3_full_typable_region-bi_fixed.vcf.gz | wc -l
bcftools view -H external-eval/HGSVC/NA12878/groundtruth-typable_truvari.vcf.gz | wc -l
bcftools view -H external-eval/HGSVC/NA12878/groundtruth-typable_truvari.vcf.gz | cut -f1-5

