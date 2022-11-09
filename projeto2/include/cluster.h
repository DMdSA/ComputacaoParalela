#ifndef CENTROID_HEADER_FILE
#define CENTROID_HEADER_FILE

/**
 * @brief cluster struct
 * 
 */
struct cluster {

    float x, y;             // (x, y) from space
    int dimension;          // number of points associated
};

typedef struct cluster *Cluster;

/**
 * @brief default print of cluster's info
 * 
 * @param c cluter's struct
 */
void print_cluster(struct cluster c);

#endif