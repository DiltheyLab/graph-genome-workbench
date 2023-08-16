
################################ Download HGSVC small variants and structural variants ################################
#
# things to be fixed:
#
#  1.) SV vcfs are malformatted in multiple ways (e.g. not sorted, fields missing, disagree with reference sequence)
#  2.) In SV vcfs, alternative alleles are not given explicitly
#  3.) small variants vcfs are given for each trio, sample needs to be extracted
#
# output:
#
# SV and snp vcfs for each sample
#
#######################################################################################################################

configfile: 'config.json'

# HGSVC data
sample_to_fam = {'HG00733': 'PR05', 'HG00514': 'SH032', 'NA19240': 'Y117'}
sample_to_pop = {'HG00733': 'PUR', 'HG00514': 'CHS' , 'NA19240': 'YRI'}
# SV files
sv_vcfs = ['HG00733.merged_nonredundant.vcf', 'HG00514.merged_nonredundant.vcf', 'NA19240.merged_nonredundant.vcf']
# snp files
snp_vcfs = ["SH032.wgs.whatshap.strandseq-10X.20160704.phased-genotypes.vcf.gz", "PR05.wgs.whatshap.strandseq-10X.20160704.phased-genotypes.vcf", "Y117.wgs.whatshap.strandseq-10X.20160704.phased-genotypes.vcf"]

# reference sequence
reference = config['reference']
samples = sample_to_fam.keys()
scripts = config['scripts']

rule download_sv_vcf:
	output:
		'generated/vcf/hgsvc/original-vcfs/sv/{sample}.merged_nonredundant.vcf'
	shell:
		'wget -O {output} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180627_PanTechnologyIntegrationSet/{wildcards.sample}.merged_nonredundant.vcf'

rule download_snp_vcf:
	output:
		'generated/vcf/hgsvc/sample-vcfs/snp/{sample}.wgs.whatshap.strandseq-10X.20160704.phased-genotypes.vcf'
	run:
		family = sample_to_fam[wildcards.sample]
		population = sample_to_pop[wildcards.sample]
		shell('wget -O - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160704_whatshap_strandseq_10X_phased_SNPs/{population}/{family}.wgs.whatshap.strandseq-10X.20160704.phased-genotypes.vcf.gz | bcftools view -s {wildcards.sample} | python3 {scripts}/hom_to_phased.py > {output}')

# fix format of SV vcfs
rule fix_vcf:
	input:
		vcf='generated/vcf/hgsvc/original-vcfs/sv/{sample}.merged_nonredundant.vcf',
		fasta=reference
	output:
		'generated/vcf/hgsvc/sample-vcfs/sv/{sample}.merged_nonredundant.vcf'
	log:
		'generated/vcf/hgsvc/sample-vcfs/sv/{sample}.merged_nonredundant.log'
	run:
		shell("grep '^#' {input.vcf} > {output}.tmp && grep -v '^#' {input.vcf} | sort -k1,1 -k2,2n >> {output}.tmp")
		shell('python3 {scripts}/fix-vcf.py {output}.tmp {input.fasta} > {output} 2> {log}')
		shell('rm {output}.tmp')


# ignore nested variants (e.g. a SNP within a SV, as they might lead to conflicts)
rule merge_samples_variation:
	input:
		snp_vcf = "generated/vcf/hgsvc/sample-vcfs/snp/{sample}.wgs.whatshap.strandseq-10X.20160704.phased-genotypes.vcf",
		sv_vcf = "generated/vcf/hgsvc/sample-vcfs/sv/{sample}.merged_nonredundant.vcf"
	output:
		'generated/vcf/hgsvc/sample-vcfs/all/{sample}-all-variants.vcf'
	params:
		tmp = 'generated/vcf/hgsvc/sample-vcfs/all/{sample}-overlap-removed.vcf.gz',
		scripts = config['scripts']
	run:
		shell('bedtools subtract -header -a {input.snp_vcf} -b {input.sv_vcf} | bgzip > {params.tmp}')
		shell('tabix -p vcf {params.tmp}')
		shell('bgzip -c {input.sv_vcf} > {input.sv_vcf}.gz')
		shell('tabix -p vcf {input.sv_vcf}.gz')
		shell('rtg vcfmerge {params.tmp} {input.sv_vcf}.gz -o - | python3 {params.scripts}/remove_clusters.py 1 > {output}')

# merge single-sample vcfs into multisample vcfs
rule merge_vcfs:
	input:
		vcfs=expand("generated/vcf/hgsvc/sample-vcfs/all/{sample}-haplotypes.vcf", sample=samples),
		reference = config['reference']
	output:
		"generated/vcf/hgsvc/multisample-vcfs/hgsvc-all-variants.vcf"
	shell:
		"{scripts}/merge_vcfs.py merge -vcf {input.vcfs} -r {input.reference} -ploidy 2 > {output}"

