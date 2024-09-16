## Graph genomes workbench repository

The folder [denbi_gg](denbi_gg/) contains the development space for Ansible playbooks and conda environment files. These are meant to be used to install graph genome software on virtual machines quickly. 

 ![Workbench logo](/logo_trans.png) [Information page](https://diltheylab.github.io/graph-genome-workbench/) 

Additionally, in the folder [evaluation_pipeline](evaluation_pipeline/) there is a pipeline to assess the performance of three state-of-the-art genotyping algorithms (PanGenie, BayesTyper and GraphTyper). The pipeline uses genomes graphs provided by the "Human Genome Structural Variant Consortium" (HGSVC) and it is applied to a standard benchmarking sample NA24385/HG002 of the "Genome In A Bottle" (GIAB).


### Requirements

- ansible (in your local machine)

We tested the installation with `ansible [core 2.14.3]`, however, other versions should work, as well. 

### Tested VMs

The playbook was tested in deNBI (Deutsches Netzwerk für Bioinformatik-Infrastruktur) VMs with the following OS:

- Ubuntu 20.04
- Ubuntu 22.04

## Acknowledgements
Funding provided by the BMBF - Förderkennzeichen 031L0184B.


