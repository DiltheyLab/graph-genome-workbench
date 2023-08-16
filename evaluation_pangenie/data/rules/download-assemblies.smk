## Shilpa's assemblies

rule download_NA12878:
	output:
		"downloaded/assemblies/NA12878/NA12878-denovo-H{number}.fa.gz"
	shell:
		"wget -O {output} ftp://ftp.dfci.harvard.edu/pub/hli/whdenovo/asm/NA12878-denovo-H{wildcards.number}.fa.gz"

rule download_NA24385:
	output:
		"downloaded/assemblies/NA24385/NA24385-denovo-H{number}.fa.gz"
	shell:
		"wget -O {output} ftp://ftp.dfci.harvard.edu/pub/hli/whdenovo/asm/NA24385-denovo-H{wildcards.number}.fa.gz"

rule download_PGP1:
	output:
		"downloaded/assemblies/PGP1/PGP1-denovo-H{number}.fa.gz"
	shell:
		"wget -O {output} ftp://ftp.dfci.harvard.edu/pub/hli/whdenovo/asm/PGP1-denovo-H{wildcards.number}.fa.gz"


## Peter's assemblies

rule download_trio_assemblies:
	output:
		fasta=temp("downloaded/assemblies/{sample, HG00731|HG00732|HG00733}/{sample}_hgsvc_pbsq2-ccs_1000-pereg.h{number}-un.racon-p2.fasta"),
		gz="downloaded/assemblies/{sample, HG00731|HG00732|HG00733}/{sample}_hgsvc_pbsq2-ccs_1000-pereg.h{number}-un.racon-p2.fasta.gz"
	run:
		shell('wget -O {output.fasta} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200417_Marschall-Eichler_NBT_hap-assm/{wildcards.sample}_hgsvc_pbsq2-ccs_1000-pereg.h{wildcards.number}-un.racon-p2.fasta')
		shell('bgzip -c {output.fasta} > {output.gz}')
		
