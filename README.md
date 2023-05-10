## Graph genomes workbench repository

The folder [denbi_gg](denbi_gg/) contains the development space for Ansible playbooks and conda environment files. These are meant to be used to install graph genome software on virtual machines quickly. 

 ![Workbench logo](/denbi_gg/logo_trans.png) [Information page](https://diltheylab.github.io/graph-genome-workbench/) 

Additionally, in the folder [evaluation_pangenie](evaluation_pangenie/) there is a pipeline showing how to use the tool PanGenie in order to genotype for varaints of a sample using short-read data as well as to measure the performance of the resulting genotyped data and a benchmark assembly of such sample.


### Requirements

- ansible (in your local machine)

We tested the installation with `ansible [core 2.14.3]`, however, other versions should work, as well. 

### Tested VMs

The playbook was tested in deNBI (Deutsches Netzwerk für Bioinformatik-Infrastruktur) VMs with the following OS:

- Ubuntu 20.04
- Ubuntu 22.04


