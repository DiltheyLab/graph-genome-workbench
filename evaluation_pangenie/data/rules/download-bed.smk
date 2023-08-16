configfile: "config.json"

repeats = config['repeats']


# bed is downloaded already, sort and merge it
rule ucsc_repeats:
	input:
		repeats
	output:
		"downloaded/bed/hg38/ucsc-simple-repeats.merged.bed"
	shell:
		"bedtools sort -i {input} | bedtools merge -i - > {output}"
