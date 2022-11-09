#include <stdio.h>

#include "../include/cluster.h"


void print_cluster(struct cluster c) {

    printf("Cluster:{.dimension=%d, .centroid={.x=%.3lf, .y=%.3f}",
                        c.dimension, c.x, c.y);
}


