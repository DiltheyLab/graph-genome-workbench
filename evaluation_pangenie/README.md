### Working pipeline for HGSVC data

However, still inconsistent as some data is taken from other folders. For instance: 

From evaluation_pangenie I’m using:
config-data:
- downloaded/fasta  “reference”
- downloaded/vcf/HGSVC  “untypable_input”
- downloaded/vcf/HGSVC  “input_graph”
- downloaded/vcf/GIAB  “truth_bed”
- downloaded/vcf/GIAB  “GIAB”
config-genotyping:
- downloaded/vcf/HGSVC  “full_callset”
- downloaded/vcf/HGSVC  “graph”
- downloaded/vcf/HGSVC  “biallelic” 
- downloaded/vcf/HGSVC  “truth”
- downloaded/vcf/GIAB    “external”
- downloaded/vcf/GIAB    “external_bed”
- downloaded/bayestyper    “bayestyper_reference_canon”
- downloaded/bayestyper    “bayestyper_reference_decoy”
From genotyping-experiments I’m using:
config-genotyping:
- genotyping/reads   “reads”
- data/callset/…/bed “bed”
- data/downloaded/fasta   “reference”


Next steps are to run the pipeline from scratch and debug the files that are missing in the main folder. After that, the branch can be merged.

