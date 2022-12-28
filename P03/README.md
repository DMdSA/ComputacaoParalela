# MPI - Lloyd's k-means simple algorithm (clusters)

> An attempt to improve the `Lloyd's k-means algorithm`, with the parallelization of processes [`SPMD`]


```bash

# compile
mpicc -g -Wall -O2 lloyd.c


# run
mpirun -np <n_processes> ./a.out <n_points> <n_clusters>

```