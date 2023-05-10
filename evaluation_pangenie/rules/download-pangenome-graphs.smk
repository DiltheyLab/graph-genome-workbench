rule download_HGSVC_GRCh38_freeze3_64haplotypes:
	output:
		"downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz"
	shell:
		"wget -O {output} https://zenodo.org/record/7763717/files/pav-panel-freeze3.vcf.gz?download=1"

rule download_HGSVC_GRCh38_freeze4_64haplotypes:
        output:
                "downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze4_64haplotypes.vcf.gz"
        shell:
                "wget -O {output} https://zenodo.org/record/7763717/files/pav-panel-freeze4.vcf.gz?download=1"

rule download_HPRC_GRCh38_88haplotypes:
        output:
                "downloaded/vcf/HPRC-GRCh38/Pangenome_graph_88haplotypes.vcf.gz"
        shell:
                "wget -O {output} https://zenodo.org/record/6797328/files/cactus_filtered_ids.vcf.gz?download=1 "

rule download_HPRC_CHM13_88haplotypes:
        output:
                "downloaded/vcf/HPRC-CHM13/Pangenome_graph_88haplotypes.vcf.gz"
        shell:
                "wget -O {output} https://zenodo.org/record/7839719/files/chm13_cactus_filtered_ids.vcf.gz?download=1 "
