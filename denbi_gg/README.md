## Deployment

In order to deploy the playbook into one/several VM/s, we need to set up the following 4 steps:

#### Ansible Installation 
- Install Ansible in your local machine following "Installing Ansible" section from the [official Ansible documentation](https://docs.ansible.com/ansible/latest/installation_guide/intro_installation.html).

#### SSH Connection: 
- Create an entry in your SSH configuration file (located at `~/.ssh/config`) in order to register the hostnames of the VMs in the deNBI Cloud and access them more quickly through SSH (without specifying passwords at every log in). An example for such an entry goes as follows:

```
Host nickname_to_access_this_host
    HostName VM_IP_address
    User ubuntu
    AddKeysToAgent yes
    IdentityFile /path/to/private_key/for/deNBI_Cloud
    Port port_to_connect_to_VM
```

Please feel free to adapt the values `nickname_to_access_this_host`, `VM_IP_address`, `/path/to/private_key/for/deNBI_Cloud` and `port_to_connect_to_VM` to your specific situation.   

#### Target VM(s) Specification
- In the inventory file, add the host nicknames of the VMs, that you provided in the previous step. There is an example of inventory file in this folder. In there, `preciseHopper` indicates the group naming of our VMs, `preciseHopper1` and `preciseHopper2` are the host nicknames we assigned to our VMs and `denbi_vms:children` show the host groups which we want to install the software in.

#### Role Installation & Execution
- Download the desired role to install the required software. There are two possibilities: *pangenomics* or *pangenotyping* role. The former includes a much wider range of pangenomic tools, whereas the latter is useful to run the evaluation pipeline of this project. To do so, run: `ansible-galaxy role install albertodescalzo.pangenotyping` 
- .
- Finally, after navigating to the [denbi_gg](https://github.com/DiltheyLab/graph-genome-workbench/tree/master/denbi_gg) folder, uncomment the line corresponding to your role in `main.yaml` and execute the playbook with the command:

`ansible-playbook main.yaml`

## Testing

The roles are tested through Github Actions before being uploaded to Ansible-Galaxy, but it is a good practice to test them locally after installation. To do so, we can test whether the software was properly installed, once the installation is over. To do so, follow the next steps:
- Access your VM: use the command `ssh "nickname_to_access_this_host"`.
- Clone the repository into your directory of preference with the command: `git clone https://github.com/DiltheyLab/graph-genome-workbench.git`
- Navigate to [denbi_gg](https://github.com/DiltheyLab/graph-genome-workbench/tree/master/denbi_gg) folder.
- Run: `./test_software.sh`. You might make the file executable first with `chmod +x test_software.sh`  This should show the different help commands for each tool. Additionally, in the first execution you might want to accept/reject the logging message for RTG Tools so that the message is not displayed during evaluation pipeline execution.


Now, your VM is ready to go. Check the [evaluation](https://github.com/DiltheyLab/graph-genome-workbench/tree/master/evaluation_pangenie) folder to run the genotype inference pipeline.


### Advance feature (for de.NBI Users)

There is a possibility to execute the pipeline in a SLURM cluster. After being given access to creation of cluster in de.NBI, and in fact creating a cluster, it is possible to share the volume you mount among the worker nodes so that the computations do not run out of memory.

To do so and have a ready cluster, follow these steps:
- Create the SSH connections to your master and workers nodes
- Modify the categories `[masternode] and [workers]` in `inventory` with your hostnames
- Download the roles that you want to install from Ansible-Galaxy. E.g. `albertodescalzo.kirgenotyping` and `albertodescalzo.pangenotyping`
- Create & Attach volumes from the de.NBI portal website and mount them **in the master node**. To do so, log in into the terminal's VM and run `sudo mkdir -p /vol/whopper && sudo mount /dev/vdb /vol/whopper`
- Change volume paths in both [variable files master](https://github.com/DiltheyLab/graph-genome-workbench/tree/master/denbi_gg/roles/master-mount-volumes/vars/main.yaml) and [variable files workers](https://github.com/DiltheyLab/graph-genome-workbench/tree/master/denbi_gg/roles/workers-mount-volumes/vars/main.yaml). Additionally, adapt the IP address of master node (this can be found by running: `cat playbook/vars/instances.yml` in the master-node's terminal).
- Uncomment the line corresponding to your role in `cluster.yaml`
- Run: `ansible-playbook cluster.yaml`
