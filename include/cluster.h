#ifndef CENTROID_HEADER_FILE
#define CENTROID_HEADER_FILE

#include "point.h"
extern const int CSIZE;

struct __attribute__ ((__packed__, __aligned__(8))) cluster {

    int dimension;
    float x, y;
};

typedef struct cluster *Cluster;

int cluster_tostring(struct cluster c, char* s);
void print_cluster(struct cluster c);

#endif