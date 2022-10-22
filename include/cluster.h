#ifndef CENTROID_HEADER_FILE
#define CENTROID_HEADER_FILE

#include "point.h"
extern const int CSIZE;

struct cluster {

    int dimension;
    struct spoint centroid;
};

typedef struct cluster *Cluster;

int cluster_tostring(struct cluster c, char* s);
void print_cluster(struct cluster c);

#endif