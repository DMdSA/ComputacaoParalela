#include <stdio.h>

#include "../include/cluster.h"


void print_cluster(struct cluster c) {

    printf("Cluster:{.dimension=%d, .centroid={.x=%.5lf, .y=%.5f, .k=%d}}",
                        c.dimension, (c.centroid).x, (c.centroid).y, (c.centroid).k);
}


