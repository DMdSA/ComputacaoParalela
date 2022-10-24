#include <stdio.h>

#include "../include/cluster.h"
const int CSIZE = 300;

int cluster_tostring(struct cluster c, char* s) {

   return snprintf(s, CSIZE, "Cluster:{.dimension=%d, .centroid={.x=%.3f, .y=%.3f, .k=%d}}",
                        c.dimension, (c.centroid).x, (c.centroid).y, (c.centroid).k);
}

void print_cluster(struct cluster c) {

    printf("Cluster:{.dimension=%d, .centroid={.x=%.3lf, .y=%.3f, .k=%d}}",
                        c.dimension, (c.centroid).x, (c.centroid).y, (c.centroid).k);
}


