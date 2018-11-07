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

 To log in via ssh:

 ```
 ssh <usuario>@omica.cicese.mx
 ```
 You will need to enter your password. 



## Future Cluster Use
 
Normally, apps are installed on: `/LUSTRE/apps` and `/LUSTRE/bioinformatica`

Working directories are created by each user on;
`/LUSTRE/bioinformatica_data/"group_subdirectory"`


