## About your "omica" cluster 
 
 "omica" is a small cluster, with master and login node 'omica',
30 processing nodes, named nodo1, nodo2, ..., nodo30
and a distributed file storage with lustre filesystem, mounted
on /LUSTRE mount point.

Jobs management and scheduler is slurm software.


For this course, we will be installing all software into:

```
 /LUSTRE/apps/workshop
```

 and will keep test data in:

```
/LUSTRE/bioinformatica_data/bioinformatica2018
```


## Logging In


If you're on a windows machine, open `PUTTY` and log into the `omica` cluster:

``
ip: 158.97.9.9
port: 22
```
You will need to enter your password. 


If you have a mac, log in via ssh instead:

```
ssh <usuario>@omica.cicese.mx
```
 You will need to enter your password. 



## Using the Installed Software

To use the software we have previously installed you'll need to execute the following:

```
echo 'export PATH=/LUSTRE/apps/workshop/miniconda3/bin:$PATH' >> ~/.bashrc
echo 'export PATH=/LUSTRE/apps/workshop/transrate-1.0.3-linux-x86_64:$PATH' >> ~/.bashrc
source ~/.bashrc
```


## Working on the cluster

After logging in, make sure to switch to a node:

```
~/.works18

```
Your prompt should now look like this: `[<usuario>@nodo11]`

Now, to enter our `conda` environment, run:

```
source activate tara
```


## Future Cluster Use
 
Normally, apps are installed on: `/LUSTRE/apps` and `/LUSTRE/bioinformatica`

Working directories are created by each user on;
`/LUSTRE/bioinformatica_data/"group_subdirectory"`


