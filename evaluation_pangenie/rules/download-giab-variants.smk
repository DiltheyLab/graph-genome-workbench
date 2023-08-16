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

rule download_BayesTyper_GRCh38_bundle:
    output:
        compressed_bundle="downloaded/bayestyper/bayestyper_GRCh38_bundle.tar.gz", 
        uncompressed_bundle=directory("downloaded/bayestyper")
    shell:
        """
        wget -O {output.compressed_bundle} http://people.binf.ku.dk/~lassemaretty/bayesTyper/bayestyper_GRCh38_bundle.tar.gz 
        tar -xvf {output.compressed_bundle} -C {output.uncompressed_bundle} 
        """

rule annotate_convert_to_biallelic_benchmark_set:
    input:
        truth="downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        callset="downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf"
    output:
        biallelic='downloaded/vcf/giab/hg38/benchmark_annotated_biallelic.vcf.gz'
    shell:
        """
        zcat {input.truth} | python3 scripts/annotate.py {input.callset} | python3 scripts/convert-to-biallelic.py {input.callset} | bgzip -c > {output.biallelic}
		tabix -f -p vcf {output.biallelic}
        """

rule split_variants_benchmark_set:
    input:
        vcf="downloaded/vcf/giab/hg38/benchmark_annotated_biallelic.vcf.gz"
    output:
        vcf="downloaded/vcf/giab/hg38/benchmark_annotated_biallelic_{variant}.vcf.gz"
    shell:
        """
        zcat {input.vcf} | python3 scripts/extract-varianttype.py {wildcards.variant} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
        """ 


##################################################           
#rule prepare_truth_set:
#    input:
#        truth="downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
#        callset="downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz"
#    output:
#        uncompressed_callset=temp("downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"),
#        biallelic='downloaded/vcf/giab/hg38/benchmark_biallelic.vcf.gz'
#    shell:
#        """
#        gunzip -c {input.truth} > {output.uncompressed_truth}
#        zcat {input.callset} | python3 scripts/ {output.uncompressed_truth} | bgzip > {output}
#		tabix -p vcf {output}
#        """
 
# rule convert_graph_to_biallelic:
#     input:
#         callset="downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz", 
#         graph="downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz"
#     output:
#         uncompressed_callset=temp("downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf"),
#         biallelic='downloaded/vcf/HGSVC-GRCh38/graph_biallelic.vcf.gz'
#     shell:
#         """
#         gunzip -c {input.callset} > {output.uncompressed_callset}
#         zcat {input.graph} | python3 scripts/convert-to-biallelic.py {output.uncompressed_callset} | bgzip > {output.biallelic}
# 		  tabix -p vcf {output.biallelic}
#         """


