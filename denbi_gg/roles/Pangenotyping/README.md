Pangenotyping
=========

Bioninformatic tools to be installed in any Linux VM in order to run the following graph-genome-based inference pipeline: https://github.com/DiltheyLab/graph-genome-workbench. For more information about the pipeline, see the GitHub page.

Requirements
------------

Role and inference pipeline were used within the de.NBI Cloud (German Network for Bioinformatics Infrastructure) in a VM with Ubuntu 20.04 and 22.04, but it should work in other cloud suppliers. 

Role Variables
--------------

Adapt the software versions to your preferences. The given variables were used to run the pipeline.   

Dependencies
------------



Example Playbook
----------------

  - hosts: denbi_vms
  roles:
    - Pangenotyping

License
-------

BSD

Author Information
------------------

Email: Alberto.Descalzo.Garcia@hhu.de