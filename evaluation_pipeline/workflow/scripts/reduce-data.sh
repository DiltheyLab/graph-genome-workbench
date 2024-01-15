#!/bin/bash

##### We need to reduce the size of: Pangenome, biallelelic callset (two samples NA12878 and NA24285), truthset(HG002), reads 
##### as long as they contain of the cases: known & unknown alleles, snp-indel, sv (different types)

chrom="chr17"

mkdir -p reduced-data
samtools faidx data/downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa $chrom > reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa > reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
zcat data/downloaded/reads/NA12878/NA12878_raw.fastq.gz | head -4000 | bgzip -c > reduced-data/NA12878_raw.fastq.gz
zcat data/downloaded/reads/NA24385/NA24385_raw.fastq.gz | head -4000 | bgzip -c > reduced-data/NA24385_raw.fastq.gz
bcftools view -r $chrom data/downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz | bgzip -c > reduced-data/Pangenome_graph_freeze3_64haplotypes.vcf.gz
bcftools view -r $chrom data/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz | bgzip -c > reduced-data/Callset_freeze3_64haplotypes.vcf.gz
bcftools view -r $chrom data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz | python workflow/scripts/prepare-for-vcfeval.py reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | bgzip -c > reduced-data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
bwa index reduced-data/GRCh38_full_analysis_set_plus_decoy_hla.fa