#ifndef CENTROID_HEADER_FILE
#define CENTROID_HEADER_FILE

#include "point.h"
extern const int CSIZE;

struct cluster {

    float x, y;
    int dimension;
};

typedef struct cluster *Cluster;

void print_cluster(struct cluster c);

#endif