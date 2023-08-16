rule download_syndip:
	output:
		'downloaded/reads/syndip/syndip_{number}.fastq.gz'
	shell:
		'wget -O {output} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/006/ERR1341796/ERR1341796_{wildcards.number}.fastq.gz'
