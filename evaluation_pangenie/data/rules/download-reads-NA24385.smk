numbers = ['00' + str(i) for i in range(1,10)] + ['0' + str(i) for i in range(10,18)]

rule download_giab_fastq:
	output:
		"downloaded/reads/NA24385/raw/D1_S1_{l}_R{r}_{n}.fastq.gz"
	shell:
		"wget -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/D1_S1_{wildcards.l}_R{wildcards.r}_{wildcards.n}.fastq.gz"

rule combine_giab_fastq:
	input:
		expand("downloaded/reads/NA24385/raw/D1_S1_{l}_R{{r}}_{n}.fastq.gz", l = ['L001', 'L002'], n = numbers)
	output:
		"downloaded/reads/NA24385/NA24385_{r}.fastq.gz"
	shell:
		"cat {input} > {output}"
