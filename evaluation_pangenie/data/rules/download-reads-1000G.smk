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
print(samples)

rule download_fastq:
	output:
		'downloaded/reads/{sample, HG00731|HG00732|NA19238|NA19239|NA19240|HG00512|HG00513|HG00514|NA20847|NA19036|HG00171|HG01571|NA12878}/{sample}_{number}.fastq.gz'
	params:
		path = lambda wildcards: sample_to_path[wildcards.sample]
	shell:
		'wget -O {output} ftp://{params.path}_{wildcards.number}.fastq.gz'
