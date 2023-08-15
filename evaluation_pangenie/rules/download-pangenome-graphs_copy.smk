### This downloads both pangenome graphs as well as the known set of variants derived called from the assemblies, i.e. callsets
scripts = config['scripts']
sample = config['sample']

subset_to_samples = {
    'only' + sample: [sample],
    'no' + sample: ["^" + sample]
}

rule download_HGSVC_GRCh38_freeze3_64haplotypes_pangenome_graph:
    output:
        vcf="downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz",
        tbi="downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz.tbi",
    shell:
        """
        wget -O {output.vcf} https://zenodo.org/record/7763717/files/pav-panel-freeze3.vcf.gz?download=1
        tabix -f -p vcf {output.vcf} > {output.tbi}
        """

rule download_HGSVC_GRCh38_freeze3_64haplotypes_callset:
    output:
        callset="downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz",
        callset_tbi="downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz.tbi"
    shell:
        """
        wget -O {output.callset} https://zenodo.org/record/7763717/files/pav-calls-freeze3.vcf.gz?download=1
        tabix -f -p vcf {output.callset} > {output.callset_tbi}
        """


rule download_HGSVC_GRCh38_freeze4_64haplotypes_pangenome_graph:
    output:
        vcf="downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze4_64haplotypes.vcf.gz",
        tbi="downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze4_64haplotypes.vcf.gz.tbi",
    shell:
        """
        wget -O {output.vcf} https://zenodo.org/record/7763717/files/pav-panel-freeze4.vcf.gz?download=1
        tabix -f -p vcf {output.vcf} > {output.tbi}
        """

rule download_HGSVC_GRCh38_freeze4_64haplotypes_callset:
    output:
        callset="downloaded/vcf/HGSVC-GRCh38/Callset_freeze4_64haplotypes.vcf.gz",
        callset_tbi="downloaded/vcf/HGSVC-GRCh38/Callset_freeze4_64haplotypes.vcf.gz.tbi"
    shell:
        """
        wget -O {output.callset} https://zenodo.org/record/7763717/files/pav-calls-freeze4.vcf.gz?download=1
        tabix -f -p vcf {output.callset} > {output.callset_tbi}
        """

rule download_HPRC_GRCh38_88haplotypes_pangenome_graph:
    output:
        vcf="downloaded/vcf/HPRC-GRCh38/Pangenome_graph_88haplotypes.vcf.gz", 
        tbi="downloaded/vcf/HPRC-GRCh38/Pangenome_graph_88haplotypes.vcf.gz.tbi",
    shell:
        """
        wget -O {output.vcf} https://zenodo.org/record/6797328/files/cactus_filtered_ids.vcf.gz?download=1 
        tabix -f -p vcf {output.vcf} > {output.tbi}
        """

rule download_HPRC_GRCh38_88haplotypes_callset:
    output:
        callset="downloaded/vcf/HPRC-GRCh38/Callset_88haplotypes.vcf.gz", 
        callset_tbi="downloaded/vcf/HPRC-GRCh38/Callset_88haplotypes.vcf.gz.tbi"
    shell:
        """
        wget -O {output.callset} https://zenodo.org/record/6797328/files/cactus_filtered_ids_biallelic.vcf.gz?download=1
        tabix -f -p vcf {output.callset} > {output.callset_tbi}
        """

rule download_HPRC_CHM13_88haplotypes_pangenome_graph:
    output:
        vcf="downloaded/vcf/HPRC-CHM13/Pangenome_graph_88haplotypes.vcf.gz",
        tbi="downloaded/vcf/HPRC-CHM13/Pangenome_graph_88haplotypes.vcf.gz.tbi",
    shell:
        """
        wget -O {output.vcf} https://zenodo.org/record/7839719/files/chm13_cactus_filtered_ids.vcf.gz?download=1 
        tabix -f -p vcf {output.vcf} > {output.tbi}
        """

rule download_HPRC_CHM13_88haplotypes_callset:
    output:
        callset="downloaded/vcf/HPRC-CHM13/Callset_88haplotypes.vcf.gz",
        callset_tbi="downloaded/vcf/HPRC-CHM13/Callset_88haplotypes.vcf.gz.tbi"
    shell:
        """
        wget -O {output.callset} https://zenodo.org/record/7839719/files/chm13_cactus_filtered_ids_biallelic.vcf.gz?download=1
        tabix -f -p vcf {output.callset} > {output.callset_tbi}
        """


rule annotate_and_biallelic_graph:
    input:
        callset="downzloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf", 
        graph="downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf"
    output:
        annotated="downloaded/vcf/HGSVC-GRCh38/graph_annotated.vcf.gz",
        biallelic="downloaded/vcf/HGSVC-GRCh38/graph_annotated_biallelic.vcf.gz"
    shell: 
        """
        cat {input.graph} | python3 {scripts}/annotate.py {input.callset} | bgzip -c > {output.annotated}
        tabix -f -p vcf {output.annotated}
        cat {output.annotated} | python3 scripts/convert-to-biallelic.py {input.callset} | bgzip -c > {output.biallelic}
		tabix -f -p vcf {output.biallelic}
        """

rule extract_sample:
    input:
        vcf=lambda wildcards: config['graph'] if wildcards.mode == 'multiallelic' else config['biallelic']
        #vcf="downloaded/vcf/HGSVC-GRCh38/graph_annotated.vcf.gz"
    output:
        vcf="downloaded/vcf/HGSVC-GRCh38/graph_annotated_{mode}_{subset}.vcf.gz",
    #wildcard_constraints:
    #    subset = "|".join([s for s in subset_to_samples.keys()]), 
    #    mode = "multiallelic|biallelic"
    params:
        samples = lambda wildcards : subset_to_samples[wildcards.subset]
    log:
        "downloaded/vcf/HGSVC-GRCh38/extract_{subset}_{mode}.log"
    shell:
        """
        bcftools view -s {params.samples} {input.vcf} | bcftools view --min-ac 1 | bgzip -c > {output.vcf}
        tabix -f -p vcf {output.vcf}
        """

########################################################################################
rule reduced_graph:
    input: 
        graph = "downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf"
    output:
        reduced_graph_vcf = "downloaded/vcf/HGSVC-GRCh38/reduced_graph.vcf.gz"
    shell:
        """
        cat {input.graph} | head -10000 | bgzip -c > {output.reduced_graph_vcf}
        tabix -f -p vcf {output.reduced_graph_vcf} 
        """

rule annotate_and_biallelic_reduced_graph:
    input:
        callset="downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf", 
        graph="downloaded/vcf/HGSVC-GRCh38/reduced_graph.vcf"
    output:
        annotated="downloaded/vcf/HGSVC-GRCh38/reduced_graph_annotated.vcf.gz",
        biallelic="downloaded/vcf/HGSVC-GRCh38/reduced_graph_annotated_biallelic.vcf.gz"
    shell: 
        """
        cat {input.graph} | python3 {scripts}/annotate.py {input.callset} | bgzip -c > {output.annotated}
        tabix -f -p vcf {output.annotated}
        zcat {output.annotated} | python3 scripts/convert-to-biallelic.py {input.callset} | bgzip -c > {output.biallelic}
		tabix -f -p vcf {output.biallelic}
        """

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


############################################################################
############################################################################
############################################################################
#rule extract_vcf:
#    input:
#        graph="downloaded/vcf/HGSVC-GRCh38/reduced_graph.vcf.gz"
#    output:


#rule reduce_callset:
#    input: 
#        input_callset = "downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz"
#    output:
#        reduced_callset_vcf = "downloaded/vcf/HGSVC-GRCh38/reduced_callset.vcf.gz",
#        reduced_callset_tbi = "downloaded/vcf/HGSVC-GRCh38/reduced_callset.vcf.gz.tbi"
#    shell:
#        """
#        zcat {input.input_callset} | cut -f1-10  | bgzip > {output.reduced_callset_vcf}
#        tabix -f -p vcf {output.reduced_callset_vcf} > {output.reduced_callset_tbi}
#        """

#rule short_reduce_graph:
#    input: 
#        reduced_graph = "downloaded/vcf/HGSVC-GRCh38/reduced_graph.vcf.gz"
#    output:
#        short_reduced_graph_vcf = "downloaded/vcf/HGSVC-GRCh38/short_reduced_graph.vcf.gz",
#        short_reduced_graph_tbi = "downloaded/vcf/HGSVC-GRCh38/short_reduced_graph.vcf.gz.tbi"
#    shell:
#        """
#        zcat {input.reduced_graph} | head -10000 | bgzip > {output.short_reduced_graph_vcf}
#        tabix -f -p vcf {output.short_reduced_graph_vcf} > {output.short_reduced_graph_tbi}
#        """

#rule short_reduce_input_callset:
#    input: 
#        reduced_callset = "downloaded/vcf/HGSVC-GRCh38/reduced_callset.vcf.gz"
#    output:
#        short_reduced_callset_vcf = "downloaded/vcf/HGSVC-GRCh38/short_reduced_callset.vcf.gz",
#        short_reduced_callset_tbi = "downloaded/vcf/HGSVC-GRCh38/short_reduced_callset.vcf.gz.tbi"
#    shell:
#        """
#        zcat {input.reduced_callset} | head -10000 | bgzip > {output.short_reduced_callset_vcf}
#        tabix -f -p vcf {output.short_reduced_callset_vcf} > {output.short_reduced_callset_tbi}
#        """


#rule convert_graph_to_biallelic:
#    input:
#        callset="downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf", 
#        graph="downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf"
#    output:
#        biallelic='downloaded/vcf/HGSVC-GRCh38/graph_biallelic.vcf.gz'
#    shell:
#        """
#        cat {input.graph} | python3 scripts/convert-to-biallelic.py {input.callset} | bgzip > {output.biallelic}
#		tabix -p vcf {output.biallelic}
#        """

#rule convert_reduced_graph_to_biallelic:
#    input:
#        callset="downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf", 
#        graph="downloaded/vcf/HGSVC-GRCh38/reduced_graph.vcf"
#    output:
#        biallelic='downloaded/vcf/HGSVC-GRCh38/graph_biallelic.vcf.gz'
#    shell:
#        """
#        cat {input.graph} | python3 scripts/convert-to-biallelic.py {input.callset} | bgzip > {output.biallelic}
#		tabix -p vcf {output.biallelic}
#        """



#rule annotate_callset:
#    input:
#        callset=config['input_callset'],
#        full_callset=config['full_callset']
#    output:
#        callset="downloaded/vcf/HGSVC-GRCh38/short_reduced_callset_annotated.vcf.gz",
#    shell: 
#        """
#        cat {input.callset} | python3 {scripts}/annotate.py {input.full_callset} | bgzip -c > {output.callset}
#        tabix -p vcf {output.callset}
#        """

#rule annotate_graph:
#    input:
#        graph=config['input_pangenome_graph'],
#        full_graph=config['full_graph']
#    output:
#        graph="downloaded/vcf/HGSVC-GRCh38/short_reduced_graph_annotated.vcf.gz"
#    shell: 
#        """
#        cat {input.graph} | python3 {scripts}/annotate.py {input.full_graph} | bgzip -c > {output.graph}
#        tabix -p vcf {output.graph}
#        """
