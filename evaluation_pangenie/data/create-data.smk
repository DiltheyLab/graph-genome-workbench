#configfile: 'config-no-alt.json'

include: "rules/assembly-vcfs.smk"
include: "rules/prepare-external-truthsets.smk"

subsets = ['all']
for sample in config['dataset']['assemblies']['leave_out_sample']:
	subsets.append('only' + sample)
	subsets.append('no' + sample)

outdir=config['outdir']

rule all:
    input:
        expand(outdir + "{truthset}/{truthset}-unique.tsv", truthset=['GIAB']),
        #outdir + "statistics/untypable-ids.tsv",
        #expand(config['outdir'] + '{truthset}/samples/{sample}.vcf.gz', truthset=['giab-small', 'giab-sv'], sample=['NA24385']),
        #expand(outdir + "multisample-vcfs/assemblies-{subset}-filtered.normalized.vcf.gz", subset = subsets),
        #expand(outdir + "multisample-vcfs/assemblies-{subset}-filtered.vcf.gz", subset = subsets),
        #expand(outdir + "multisample-vcfs/assemblies-{subset}-biallelic-filtered.vcf.gz", subset = subsets),
        #	expand(outdir + "statistics/bcftools-plots/summary-{subset}.pdf", subset=subsets),
        #expand(outdir + "statistics/vcftools-plots/indel-histogram-{subset}.pdf", subset=subsets),
        #expand(outdir + "statistics/vcfstats-{subset}.txt", subset=subsets),
        #expand(outdir + "statistics/vcfstats-plots/{plots}-{subset}.pdf", plots=["variant-numbers", "het-hom"], subset=subsets),
        #expand(outdir + "statistics/untypable-{mode}.tsv", mode=['ids', 'bubbles']),
        #expand(outdir + "statistics/raw-{subset}.tsv", subset=subsets),
        #expand(outdir + "statistics/raw-callable-{subset}.tsv", subset=subsets),
        #expand(outdir + "statistics/filtered-{subset}.tsv", subset=subsets),
        #expand(outdir + "statistics/merged-filtered-{subset}.tsv", subset=subsets),
        #expand(outdir + "statistics/alleles-per-bubble-{subset}.pdf", subset=subsets),
        #expand(outdir + "statistics/complex-bubbles-{subset}.bed", subset=subsets),
        #expand(outdir + "bed/callable-regions-missing.bed"),
        #expand(outdir + "bed/callable-regions.bed"),
        #expand(outdir + "bed/{sample}_callable.bed", sample=config['dataset']['assemblies']['leave_out_sample']),
        #expand(outdir + "{truthset}/{truthset}-unique.tsv", truthset=['giab-small', 'giab-sv']),
        #expand(outdir + "{truthset}/{truthset}-annotated-sv.vcf", truthset=['giab-sv']),
        #expand(outdir + "{truthset}/{truthset}-unique-sv.vcf", truthset=['giab-sv']),
        #outdir + "giab-small/giab-small-annotated-indel.vcf",
        #outdir + "giab-small/giab-small-unique-indel.vcf",
        #outdir + "giab-sv/comparison/giab-sv-calls-sv/summary.txt",
        #outdir + "giab-small/comparison/giab-small-calls-indel/summary.txt"
