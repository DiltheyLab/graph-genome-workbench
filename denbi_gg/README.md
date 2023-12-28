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

#### Playbook Execution
- Finally, after navigating to the [denbi_gg](https://github.com/DiltheyLab/graph-genome-workbench/tree/master/denbi_gg) folder, execute the playbook with the command:

`ansible-playbook main.yaml`

## Testing

Additionally, we can test whether the software was properly installed, once the installation is over. To do so, follow the next steps:
- Access your VM: use the command `ssh "nickname_to_access_this_host"`.
- Clone the repository into your directory of preference with the command: `git clone https://github.com/DiltheyLab/graph-genome-workbench.git`
- Navigate to [denbi_gg](https://github.com/DiltheyLab/graph-genome-workbench/tree/master/denbi_gg) folder.
- Run: `./test_software.sh`. You might make the file executable first with `chmod +x test_software.sh`  This should show the different help commands for each tool. Additionally, in the first execution you might want to accept/reject the logging message for RTG Tools so that the message is not displayed during evaluation pipeline execution.


Now, your VM is ready to go. Check the [evaluation](https://github.com/DiltheyLab/graph-genome-workbench/tree/master/evaluation_pangenie) folder to run the genotype inference pipeline.
