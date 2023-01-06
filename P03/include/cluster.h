#ifndef CLUSTER_FILE_HEADER
#define CLUSTER_FILE_HEADER

/**
 * @brief cluster struct
 * 
 */
struct cluster {

    float x, y;             // (x, y) from space
    int dimension;          // number of points associated
};

typedef struct cluster *Cluster;

#endif