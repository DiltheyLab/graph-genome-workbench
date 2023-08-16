include: "rules/download-references.smk"
include: "rules/download-reads-NA24385.smk"
include: "rules/download-reads-1000G.smk"
include: "rules/download-assemblies.smk"
include: "rules/download-giab-variants.smk"
include: "rules/download-bed.smk"
include: "rules/download-reads-syndip.smk"

rule all:
    input:
        "downloaded/bed/hg38/ucsc-simple-repeats.merged.bed",

        # references
        ##'downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa',
        ##'downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai',
        ##'downloaded/fasta/no-alt/no-alt.fa',
        ##'downloaded/fasta/no-alt/no-alt.fa.fai',
        ##'downloaded/bed/hg38/ucsc-simple-repeats.merged.bed',

        # GIAB ground truth (small variants)
        ##'downloaded/vcf/giab/hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz',
        ##'downloaded/vcf/giab/hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed',

        # GIAB SVs in medically relevant genes
        ##'downloaded/vcf/giab/hg38/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz',
        ##'downloaded/vcf/giab/hg38/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.bed',

        # GIAB dipcall variants
        ##'downloaded/vcf/giab/hg38/HG002v11-align2-GRCh38.dip.vcf.gz',
        ##'downloaded/vcf/giab/hg38/HG002v11-align2-GRCh38.dip.bed',

        # reads
        ##expand("downloaded/reads/{sample}/{sample}_{number}.fastq.gz", number=[1,2], sample=["NA24385"])
