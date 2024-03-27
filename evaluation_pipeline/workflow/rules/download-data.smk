# stores paths to reads
gdown = config['programs']['gdown']
samtools = config['programs']['samtools']
bedtools = config['programs']['bedtools']
read_paths = {}

for line in open(config['reads'], 'r'):
	if line.startswith('#'):
		continue
	fields = line.strip().split()
	sample_id = fields[1]
	read_path = fields[7]
	read_paths[sample_id] = read_path


rule download_data:
    input:
        expand("data/downloaded/bed/hg38/{callset}/ucsc-simple-repeats.merged.bed", callset = [c for c in config['callsets'].keys()]),

        # references + indexing
        'data/downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa',
        'data/downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai',
        "data/downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann",
        
        # Input Pangenome Graph (PanGenie Github)
        'data/downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz',
        
        # Callset (PanGenie Github)
        'data/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz',
         
        # reads
        expand("{sample_paths}", sample_paths=list(read_paths.values())),

        ## Download GRCh38 bundle for BayesTyper
        "data/downloaded/bayestyper_utils/bayestyper_GRCh38_bundle.tar.gz",

        ## Download HG002/NA24385 benchmark external validation dataset
        "data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
        "data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.bed"


        

####################################################################
#########################   Download bed     #######################
####################################################################

# bed is downloaded already, sort and merge it
rule ucsc_repeats:
	input:
		lambda wildcards: config['callsets'][wildcards.callset]['repeat_regions']
	output:
		"data/downloaded/bed/hg38/{callset}/ucsc-simple-repeats.merged.bed"
	shell:
		"{bedtools} sort -i {input} | {bedtools} merge -i - > {output}"



####################################################################
##########   Download Pangenome graphs and Callsets     ############
####################################################################

rule download_HGSVC_GRCh38_freeze3_64haplotypes_pangenome_graph:
    output:
        vcf="data/downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz",
        tbi="data/downloaded/vcf/HGSVC-GRCh38/Pangenome_graph_freeze3_64haplotypes.vcf.gz.tbi",
    shell:
        """
        wget -O {output.vcf} https://zenodo.org/record/7763717/files/pav-panel-freeze3.vcf.gz?download=1
        tabix -f -p vcf {output.vcf} > {output.tbi}
        """

rule download_HGSVC_GRCh38_freeze3_64haplotypes_callset:
    output:
        callset="data/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz",
        callset_tbi="data/downloaded/vcf/HGSVC-GRCh38/Callset_freeze3_64haplotypes.vcf.gz.tbi"
    shell:
        """
        wget -O {output.callset} https://zenodo.org/record/7763717/files/pav-calls-freeze3.vcf.gz?download=1
        tabix -f -p vcf {output.callset} > {output.callset_tbi}
        """

rule download_HPRC_GRCh38_88haplotypes_pangenome_graph:
    output:
        vcf="data/downloaded/vcf/HPRC-GRCh38/Pangenome_graph_88haplotypes.vcf.gz", 
        tbi="data/downloaded/vcf/HPRC-GRCh38/Pangenome_graph_88haplotypes.vcf.gz.tbi",
    shell:
        """
        wget -O {output.vcf} https://zenodo.org/record/6797328/files/cactus_filtered_ids.vcf.gz?download=1 
        tabix -f -p vcf {output.vcf} > {output.tbi}
        """

rule download_HPRC_GRCh38_88haplotypes_callset:
    output:
        callset="data/downloaded/vcf/HPRC-GRCh38/Callset_88haplotypes.vcf.gz", 
        callset_tbi="data/downloaded/vcf/HPRC-GRCh38/Callset_88haplotypes.vcf.gz.tbi"
    shell:
        """
        wget -O {output.callset} https://zenodo.org/record/6797328/files/cactus_filtered_ids_biallelic.vcf.gz?download=1
        tabix -f -p vcf {output.callset} > {output.callset_tbi}
        """


####################################################################
######################   Download reads  GIAB   ####################
####################################################################

### SAMPLE 24385
numbers = ['00' + str(i) for i in range(1,10)] + ['0' + str(i) for i in range(10,18)]

rule download_giab_fastq:
	output:
		temp("data/downloaded/reads/NA24385/raw/D1_S1_{l}_R{r}_{n}.fastq.gz")
	shell:
		"wget -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/D1_S1_{wildcards.l}_R{wildcards.r}_{wildcards.n}.fastq.gz"

rule combine_giab_fastq:
	input:
		expand("data/downloaded/reads/NA24385/raw/D1_S1_{l}_R{r}_{n}.fastq.gz", l = ['L001', 'L002'], r = ['1', '2'], n = numbers)
	output:
		"data/downloaded/reads/NA24385/NA24385_raw.fastq.gz"
	shell:
		"cat {input} > {output}"


####################################################################
#####################   Download reads  1000G   ####################
####################################################################

sample_to_path = {
	'HG00731': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/004/ERR3241754/ERR3241754'],
	'HG00732': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/005/ERR3241755/ERR3241755'],
	'NA19238': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/003/ERR3239453/ERR3239453'],
	'NA19239': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/004/ERR3239454/ERR3239454'],
	'NA19240': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR398/000/ERR3989410/ERR3989410'],
	'HG00512': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR398/000/ERR3988780/ERR3988780'],
	'HG00513': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/004/ERR3241684/ERR3241684'],
	'HG00514': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR398/001/ERR3988781/ERR3988781'],
	'NA20847': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/009/ERR3239999/ERR3239999'],
	'NA19036': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/005/ERR3239695/ERR3239695'],
	'HG00171': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/000/ERR3240160/ERR3240160'],
	'HG01571': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/006/ERR3241986/ERR3241986'],
	'NA12878': ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR323/004/ERR3239334/ERR3239334']
}

samples = [s for s in sample_to_path.keys()]

rule download_fastq:
	output:
		'data/downloaded/reads/{sample}/{sample}_{number}.fastq.gz'
	wildcard_constraints:
		sample="HG00731|HG00732|NA19238|NA19239|NA19240|HG00512|HG00513|HG00514|NA20847|NA19036|HG00171|HG01571|NA12878",
		number="1|2"
	params:
		path = lambda wildcards: sample_to_path[wildcards.sample]
	shell:
		'wget -O {output} ftp://{params.path}_{wildcards.number}.fastq.gz'


rule combine_fastq:
	input:
		lambda wildcards: expand("data/downloaded/reads/{{sample}}/{{sample}}_{number}.fastq.gz", number = ['1', '2'])
	output:
		"data/downloaded/reads/{sample}/{sample}_raw.fastq.gz"
	shell:
		"cat {input} > {output}"


####################################################################
################   Download references + indexing    ###############
####################################################################

bwa=config['programs']['bwa']

rule download_reference_hg38:
	output:
		fasta="data/downloaded/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa"
	shell:
		"wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

rule index_hg38_reference:
	input:
		"data/downloaded/fasta/{filename}.fa"
	output:
		"data/downloaded/fasta/{filename}.fa.fai"
	shell:
		"{samtools} faidx {input}"

rule bwa_index_fasta:
	input:
		"data/downloaded/fasta/{filename}.fa"
	output:
		"data/downloaded/fasta/{filename}.fa" + ".ann"
	log:
		"logs/data/downloaded/fasta/{filename}-indexing.log"
	resources:
		mem_total_mb=5000
	shell:
		"(/usr/bin/time -v {bwa} index {input}) &> {log}"


####################################################################
#######################   Download utils     #######################
####################################################################

### Following data is downloaded from Google Drive, since the original link doesn't work anymore. 
### --> See issue: https://github.com/bioinformatics-centre/BayesTyper/issues/48

rule download_BayesTyper_GRCh38_bundle:
    output:
        compressed_bundle="data/downloaded/bayestyper_utils/bayestyper_GRCh38_bundle.tar.gz", 
        uncompressed_bundle=directory("data/downloaded/bayestyper_utils")
    params:
        file_id = "1ioTjLFkfmvOMsXubJS5_rwpfajPv5G1Q"
    shell:
        """
        {gdown} {params.file_id} -O {output.compressed_bundle}
        tar -xvf {output.compressed_bundle} -C {output.uncompressed_bundle} 
        """


####################################################################
#####   Download benchmark external validation dataset     #########
####################################################################

# download GIAB benchmark VCF file for HG002_NA24385
rule download_GIAB_bechmark_HG002:
	output:
		vcf="data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
		tbi="data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
		bed="data/downloaded/vcf/giab/hg38/HG002_GRCh38_1_22_v4.2.1_benchmark.bed"
	run:
		shell("wget -O {output.vcf} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
		shell("wget -O {output.tbi} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi")
		shell("wget -O {output.bed} https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")
