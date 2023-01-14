/**
 * @file initialization.c
 * @author group a93313_a93180
 * @brief testing lloyd's initialization technique for random points and their clusters - mpi approach
 * @version 1.0
 * @date 2022-12-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <stddef.h>
#include <mpi.h>
#include <assert.h>

#include "../include/point.h"
#include "../include/cluster.h"

// pointer to array of sample's points and clusters
struct point* RANDOM_SAMPLE;
struct cluster* CLUSTERS;


/**
 * @brief populate points and clusters' arrays with random values
 * 
 * @param n_points number of points
 * @param n_clusters number of clusters
 */
void populate(int n_points, int n_clusters) {

   // default rand seed
    srand(10);
    
    // random points + first k clusters
    int i = 0;
    for (; i < n_clusters; i++) {

        float x = (float) rand() / RAND_MAX, y = (float) rand() / RAND_MAX;
        RANDOM_SAMPLE[i] = (struct point) {.x = x, .y = y, .cluster = -1};
        CLUSTERS[i] = (struct cluster) {.x = x, .y = y, .dimension = 0};
    }

    // remaining points
    for (int j = i; j < n_points; j++) {

        float x = (float) rand() / RAND_MAX, y = (float) rand() / RAND_MAX;
        RANDOM_SAMPLE[j] = (struct point) {.x = x, .y = y, .cluster = -1};
    }
}



/**
 * @brief gets the distance between two points, without square rooting.
 * As this distance is only used as a means of comparison, it is not necessary to use the square root.
 * 
 * @param a point a
 * @param b cluster b
 * @return float distance
 */
float distanceComparator(struct point a, struct cluster b) {

    float x = (b.x - a.x);
    x *= x;
    float y = (b.y - a.y);
    y *= y;

    return y+x;
}



/**
 * @brief associates each sample's point with its corresponding cluster
 * 
 * @param SAMPLE_SUBSET set of points
 * @param changes auxiliar array to save dimension changes on clusters
 * @param number_points_per_proccess number of points being processed
 * @param n_clusters number of clusters do evaluate
 * @return int if any change has been made
 */
int update_samples(struct point* SAMPLE_SUBSET, float* changes, int number_points_per_proccess, int n_clusters ){

    int changes_flag = 0;
    
    for(int i = 0; i < number_points_per_proccess; i++){
            
            float minDist = FLT_MAX;
            int minK = 0, lastK = 0;

            struct point p = SAMPLE_SUBSET[i];
            minK = lastK = p.cluster;
            struct cluster c = {0};
            
            // for each point, evaluate all clusters
            for(int d = 0; d < n_clusters; d++){
                
                c = CLUSTERS[d];
                // distance measure
                float dist = distanceComparator(p,c);
                
                // if a nearest cluster is found, get its info
                if(minDist > dist ){
                    
                    minDist = dist;
                    minK = d;
                }
            }

            // if a nearest cluster was in fact found (int relation to the previous one)
            if (minK != lastK){ 
                
                // update point's info
                SAMPLE_SUBSET[i].cluster = minK;
                changes_flag = 1;
            }

            // update auxiliar array by adding coordinates of the point to its cluster and incrementing dimension
            // {d, x, y} , {d, x, y} , ...
            ++changes[minK * 3];
            changes[minK * 3 + 1] += p.x; 
            changes[minK * 3 + 2] += p.y; 
        }

    return changes_flag;
}



/**
 * @brief same as update_samples, only needed to process the first iteration of points, as their clusters are initialized as (-1)
 * 
 * @param SAMPLE_SUBSET set of points
 * @param changes auxiliar array to save dimension changes on clusters
 * @param number_points_per_proccess number of points being processed
 * @param n_clusters number of clusters do evaluate
 */
void first_update_samples(struct point* SAMPLE_SUBSET, float* changes, int number_points_per_proccess, int n_clusters ){

    //First cicle because of the -1 in the points' default klusters
    for(int i = 0; i < number_points_per_proccess; i++){
            
        float minDist = FLT_MAX;
        int minK = 0, lastK = 0;

        struct point p = SAMPLE_SUBSET[i];
        minK = lastK = p.cluster;
        struct cluster c = {0};
        
        // for each point, evaluate all clusters
        for(int d = 0; d < n_clusters; d++){
            
            c = CLUSTERS[d];
            
            // distance measure
            float dist = distanceComparator(p, c);
            
            // if a nearest cluster is found, get its info
            if(minDist > dist){
                
                minDist = dist;
                minK = d;
            }
        }

        // if a nearest cluster was in fact found (int relation to the previous one)
        if(minK != lastK){
        
            SAMPLE_SUBSET[i].cluster = minK;
        }
        // {d, x, y}, {d, x, y}
        ++changes[minK * 3];
        changes[minK * 3 + 1] += p.x;   //Somar o X e o Y 
        changes[minK * 3 + 2] += p.y; 
    }
}



/**
 * @brief update clusters' centroids
 * 
 * @param changes clusters' points info on auxiliar array
 * @param nclusters number of clusters to process
 */
void updateCentroids(float* changes, int nclusters){
        
    for (int i = 0; i < nclusters; i++) {
        
        int dimension = (CLUSTERS[i].dimension) = (int) changes[i*3];

        (CLUSTERS[i]).x = changes[i*3 + 1] / (float)dimension;
        (CLUSTERS[i]).y = changes[i*3 + 2] / (float)dimension;
    }
}



/**
 * @brief Create a Mpi Point object for POINT struct
 * 
 * @return MPI_Datatype 
 */
MPI_Datatype createMpiPoint() {
    
    const int p_nitems=3;
    int          p_blocklengths[3] = {1,1,1};
    MPI_Datatype p_types[3] = {MPI_FLOAT, MPI_FLOAT, MPI_INT};
    MPI_Datatype mpi_point;
    MPI_Aint     p_offsets[3];

    p_offsets[0] = offsetof(struct point, x);
    p_offsets[1] = offsetof(struct point, y);
    p_offsets[2] = offsetof(struct point, cluster);

    MPI_Type_create_struct(p_nitems, p_blocklengths, p_offsets, p_types, &mpi_point);
    MPI_Type_commit(&mpi_point);
    return mpi_point;
}



/**
 * @brief Create a Mpi Cluster object for CLUSTER struct
 * 
 * @return MPI_Datatype 
 */
MPI_Datatype createMpiCluster() {

    const int c_nitems=3;
    int          c_blocklengths[3] = {1,1,1};
    MPI_Datatype c_types[3] = {MPI_FLOAT, MPI_FLOAT, MPI_INT};
    MPI_Datatype mpi_cluster;
    MPI_Aint     c_offsets[3];

    c_offsets[0] = offsetof(struct cluster, x);
    c_offsets[1] = offsetof(struct cluster, y);
    c_offsets[2] = offsetof(struct cluster, dimension);

    MPI_Type_create_struct(c_nitems, c_blocklengths, c_offsets, c_types, &mpi_cluster);
    MPI_Type_commit(&mpi_cluster);
    return mpi_cluster;
}



/**
 * @brief checks correct usage of program
 * 
 * @param argc number of arguments
 * @param argv arguments
 */
void main_arg_control(int argc, char* argv[]) {

    if (argc != 3) {

        printf("\n#> [%s:error]\n\t./%s <number of points> <number of clusters>\n\n", argv[0],argv[0]);
        exit(1);
    }
}



int main(int argc, char* argv[]) {

    // controll main arguments
    //main_arg_control(argc, argv);
    
    int n_points = 0, n_clusters = 0;
    sscanf(argv[1], "%d", &n_points);
    sscanf(argv[2], "%d", &n_clusters);

    if (n_clusters > n_points) {

        fprintf(stderr, "\n#> Fatal: you can't have more clusters than points.\n");
        exit(3);
    }

    // ---------------------------------------------------------


    int size = 0, rank = 0;
    
    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    
    // commit struct dataypes for mpi
    /* MPI_POINT */
    MPI_Datatype mpi_point = createMpiPoint();
    
    /* MPI_CLUSTER */
    MPI_Datatype mpi_cluster = createMpiCluster();
    
    // get processes' ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // clusters space must be allocated for all processes
    CLUSTERS = (struct cluster*) malloc(sizeof(struct cluster) * n_clusters);
    
    if (!CLUSTERS) {
            fprintf(stderr, "\n#> Fatal: failed to allocate requested bytes for clusters.\n");
            exit(2);
    }
    
    if (rank == 0) {
        
        // only "master"/root process will initialize and populate the
        // random numbers sample. The clusters will only be populated

        RANDOM_SAMPLE = (struct point*) malloc(sizeof(struct point) * n_points);
    
        if (!RANDOM_SAMPLE) {
            fprintf(stderr, "\n#> Fatal: failed to allocate requested bytes for samples.\n");
            exit(2);
        }

        populate(n_points, n_clusters);
    }

    /**
     * for each process, create a buffer that will hold a subset/chunk of the entire
     * points sample
    */
    int number_points_per_proccess = n_points / size;

    struct point* SAMPLE_SUBSET = (Point) malloc(sizeof(struct point) * number_points_per_proccess);
    assert(SAMPLE_SUBSET != NULL);

    /**
     * for each process, receive the clusters via broadcast
     * 
     * "for MPI collective communicators, _everyone_ has to participate; everyone
     * has to call the Bcast"
    */

    MPI_Bcast(CLUSTERS, n_clusters, mpi_cluster, 0, MPI_COMM_WORLD);

    /**
     * scatter the random numbers from the root process to all processes in the
     * MPI World, via communicator
     * 
     * int MPI_Scatter(const void *sendbuf, int sendcount
     *                  , MPI_Datatype sendtype, void *recvbuf
     *                  , int recvcount, MPI_Datatype recvtype
     *                  , int root, MPI_Comm comm)
    */
   
    MPI_Scatter(RANDOM_SAMPLE, number_points_per_proccess
                , mpi_point, SAMPLE_SUBSET
                , number_points_per_proccess, mpi_point
                , 0, MPI_COMM_WORLD);


    // {{d0, x0, y0},  {d1, x1, y1},  {d2, x2, y1},  {d3, x3, y3}}
    const int arraysize = n_clusters*3;
    float changes[arraysize];
    
    for(int i = 0; i < arraysize; i++)
        changes[i] = 0.0f;

    
    //First cicle because of the -1 in the points' default klusters
    first_update_samples(SAMPLE_SUBSET, changes, number_points_per_proccess, n_clusters);

    // update each cluster with new dimensions (after first iteration)
    for(int i = 0; i < n_clusters; i++)
        CLUSTERS[i].dimension = (int) changes[i*3];
    

    // all processes get the same array with changes from all processes
    MPI_Allreduce(MPI_IN_PLACE, changes, arraysize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    // all processes update their centroids (it's the exact same work for everyone...)
    updateCentroids(changes,n_clusters);


    int change_flag = 1, iterations = 1;

    // while changes are being made
    while(change_flag) {

        // reset change's array
        for(int i = 0; i < arraysize; i++) {
            changes[i] = 0.0f;
        }

        // update points' clusters
        change_flag = update_samples(SAMPLE_SUBSET, changes, number_points_per_proccess, n_clusters);
        
        // has any change been made by any of the clusters?
        MPI_Allreduce(MPI_IN_PLACE, &change_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // if so, allreduce the changes arrays from all processes and then update the centroids
        if (change_flag != 0){

            MPI_Allreduce(MPI_IN_PLACE, changes, arraysize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            updateCentroids(changes, n_clusters);
            iterations++;
        }
    }
    
    // when the algorithm converges, get the point's info from all processes back to the root
    MPI_Gather(SAMPLE_SUBSET, number_points_per_proccess, mpi_point,
                        RANDOM_SAMPLE,number_points_per_proccess,
                        mpi_point, 0, MPI_COMM_WORLD);
    

    if (rank == 0) {
        for (int i = 0; i < n_clusters; i++) {
    
            printf("\nCenter: (%.3f, %.3f) : Size: %d", CLUSTERS[i].x, CLUSTERS[i].y, CLUSTERS[i].dimension);
        }
        printf("\n#> n iterations: %d\n", iterations);
        free(RANDOM_SAMPLE);
    }
    free(SAMPLE_SUBSET);
    free(CLUSTERS);
    MPI_Finalize();

    return 0;
}
