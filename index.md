---
layout: default
---

![GGW_logo](/logo_trans.png){: width="450"}

## Software

This is a list of software that is being installed when you use this image. If you have any issues using this playbook on a deNBI virtual machine or any other system, please write an [issue](https://github.com/DiltheyLab/graph-genome-workbench/issues).

| software      | command line |
|:--------------|:-------------|
| [Miniconda](https://docs.conda.io/en/latest/miniconda.html) | `conda` |
| [Pangenie](https://github.com/eblerjana/PanGenie) | `PanGenie, Pangenie-graph` |
| [Minimap2](https://github.com/lh3/minimap2) | `minimap2` |
| [Miniasm](https://github.com/lh3/miniasm) | `miniasm` |
| [Minigraph](https://github.com/lh3/minigraph) | `minigraph` |
| [vg](https://github.com/vgteam/vg) | `vg` |
| [odgi](https://github.com/pangenome/odgi) | `odgi` |
| [graphtyper2](https://github.com/DecodeGenetics/graphtyper) | `graphtyper` |
| [Bandage](https://github.com/rrwick/Bandage) | `Bandage` |
| [Seqwish](https://github.com/ekg/seqwish) | `seqwish` |
| [Gfaffix](https://github.com/marschall-lab/GFAffix) | `gfaffix` |
| [Smoothxg](https://github.com/pangenome/smoothxg) | `smoothxg` |
| [Wfmash](https://github.com/waveygang/wfmash) | `wfmash` |
| [Pggb](https://github.com/pangeonme/pggb) | `pggb` |
| [IGV](https://software.broadinstitute.org/software/igv/) | `IGV` |
| [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) | see below |
| [RTG Tools](https://www.realtimegenomics.com/products/rtg-tools) | `/home/ubuntu/rtg-tools-3.12.1/rtg` |
| [SnpEff](https://netcologne.dl.sourceforge.net/project/snpeff) | `java -jar /home/ubuntu/snpEff/snpEff.jar` |
| [ClinEff](https://netcologne.dl.sourceforge.net/project/snpeff) | `java -jar /home/ubuntu/clinEff/ClinEff.jar`  |
| [Gfaffix](https://github.com/marschall-lab/GFAffix) | `gfaffix` |
| [Novograph](https://github.com/NCBI-Hackathons/NovoGraph) |  |

All software with exception of cactus is available directly from the command line. You can see all executables in `/usr/local/bin/`. Cactus has its own environment because it comes packaged with its own version of some of the software in the list. In order to use cactus and the versions of the software that cactus relies on, activate the conda environment `conda activate cactus`.

Bandage and IGV are GUI software. This playbook also installs an X-server which can be used via ssh when you provide the `-X` parameter.

