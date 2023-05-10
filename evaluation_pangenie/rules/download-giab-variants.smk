# download GIAB small variants (hg38)
rule download_GIAB_variants_hg38:
	output:
		vcf="downloaded/vcf/giab/hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz",
		tbi="downloaded/vcf/giab/hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi",
		bed="downloaded/vcf/giab/hg38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed"
	run:
		shell("wget -O {output.vcf} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz")
		shell("wget -O {output.tbi} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi")
		shell("wget -O {output.bed} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv3.3.2/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed")


# download GIAB SVs for med. relevant genes (hg38)
rule download_GIAB_medically_relevant:
	output:
		vcf="downloaded/vcf/giab/hg38/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz",
		tbi="downloaded/vcf/giab/hg38/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz.tbi",
		bed="downloaded/vcf/giab/hg38/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.bed"
	run:
		shell("wget -O {output.vcf} ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz")
		shell("wget -O {output.tbi} ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz.tbi")
		shell("wget -O {output.bed} ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_medical_genes_SV_benchmark_v0.01/HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.bed")


# download GIAB dipcall SVs for HG002
rule download_GIAB_dipcall_svs:
	output:
		vcf="downloaded/vcf/giab/hg38/HG002v11-align2-GRCh38.dip.vcf.gz",
		tbi="downloaded/vcf/giab/hg38/HG002v11-align2-GRCh38.dip.vcf.gz.tbi",
		bed="downloaded/vcf/giab/hg38/HG002v11-align2-GRCh38.dip.bed"
	run:
		shell("wget -O {output.vcf} ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SupplementaryFiles/HG002v11-align2-GRCh38/HG002v11-align2-GRCh38.dip.vcf.gz")
		shell("wget -O {output.tbi} ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SupplementaryFiles/HG002v11-align2-GRCh38/HG002v11-align2-GRCh38.dip.vcf.gz.tbi")
		shell("wget -O {output.bed} ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SupplementaryFiles/HG002v11-align2-GRCh38/HG002v11-align2-GRCh38.dip.bed")

# download GIAB benchmark VCF file for HG002_NA24385
rule download_GIAB_bechmark_HG002:
	output:
		vcf="downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
		tbi="downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
		bed="downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.bed"
	run:
		shell("wget -O {output.vcf} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
		shell("wget -O {output.tbi} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi")
		shell("wget -O {output.bed} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")
