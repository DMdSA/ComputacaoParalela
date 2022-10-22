#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>

#include "../include/point.h"
#include "../include/cluster.h"

#include <time.h>
#include <math.h>

#define N 10000000
#define K 4


struct spoint* RANDOM_SAMPLE;
struct cluster* CLUSTERS; 


/**
 * @brief inicializa N amostras, no espaço, com valores random, e inicializa K clusters com random centroids
 * 
 * @param n_points número de amostras a utilizar
 * @param n_clusters número de clusters
 */
void inicializa(int n_points, int n_clusters) {

    // malloc for vector of size N
    RANDOM_SAMPLE = (struct spoint*) malloc(sizeof(struct spoint) * N);

    // malloc for N clusters
    CLUSTERS = (struct cluster*) malloc(sizeof(struct cluster) * K);

    
    // if malloc fails
    if (!RANDOM_SAMPLE || !CLUSTERS) {
        fprintf(stderr, "Fatal: failed to allocate bytes.\n");
        exit(1);
    }

    // rand seed
    srand(10);

    // random sample generator
    for (int i = 0; i < N; i++) {

        RANDOM_SAMPLE[i] = (struct spoint) {.x = (float) rand()/RAND_MAX, .y = (float) rand()/RAND_MAX, .k = -1};
    }

    // initialize each of K clusters and assign their first centroid
    for (int i = 0; i < K; i++) {

        CLUSTERS[i] = (struct cluster) {.dimension = 0, .centroid = {.x = RANDOM_SAMPLE[i].x, .y = RANDOM_SAMPLE[i].y, i}};
    }
}
 


/**
 * @brief goes through all samples and updates their cluster
 * 
 * @return int 0 if there was no changes were made
 */
int update_samples() {

    // auxiliar variables
    float minDist = FLT_MAX, dist = 0.0f;
    int minK = INT_MAX, changes_flag = 0;
    struct spoint p;

    // for each of the samples
    for (int i = 0; i < N; i++) {

        // get current point
        p = RANDOM_SAMPLE[i];

        // default values for minimum calculation
        
        minDist = INT_MAX;
        minK = p.k;

        // for each of the *other* clusters
        for (int j = 0; j < K; j++) {

            // calculate the euclidian distance
            dist = euclidianDistance(p, (CLUSTERS[j]).centroid);
            // if a new minimum is found
            if (dist < minDist) {

                // update the minimum distance and it's associated cluster
                minDist = dist;
                minK = j;
            }
        }


        if (minK != (p.k)) {
            
            // assuming k will always be >= 0
            // update previous cluster
            if (p.k != -1)
                (CLUSTERS[p.k]).dimension--;
            // update newer cluster
            (CLUSTERS[minK]).dimension++;
            // update point's cluster
            (RANDOM_SAMPLE[i]).k = minK;
            changes_flag = 1;
        }
    }

    return changes_flag;
}



void update_centroids() {

    /*
        k = 4, #tam = 20 bytes (cada ponto), chacheL1 = 64bytes     => 3.2
        k = 4, #tam = 16 (cada ponto), cacheL1 = 64                 => 4
    */

    // initialize means array
    int means_size = K*2;
    float means_array[means_size];
    for (int i = 0; i < means_size; i++) {
        means_array[i] = 0.0f;
    }

    struct spoint p;
    int index = 0;

    // for each of the samples
    for (int i = 0; i < N; i++) {

        p = RANDOM_SAMPLE[i];
        index = p.k * 2;
        means_array[index] += p.x;
        means_array[index + 1] += p.y;

        // k = 0, 0 1
        // k = 1, 2 3
        // k = 2, 4 5
    }

    struct cluster c;
    int dimension = 0;

    // for each cluster, calculate the new centroid
    for (int i = 0; i < K; i++) {
        
        // aux variables
        index = i*2;
        dimension = (CLUSTERS[i]).dimension;
        // mean calculation
        //printf("\n#> dimension %d, index %d, %.5f, %.5f", dimension, index, (CLUSTERS[i].centroid).x, (CLUSTERS[i].centroid).y);
        (CLUSTERS[i].centroid).x = means_array[index] / dimension; 
        (CLUSTERS[i].centroid).y = means_array[index + 1] / dimension;

    }

}



int main() {

    // debug purposes
    char spoint_string[SPSIZE];
    char cluster_string[CSIZE];
    
    // initialize random samples + clusters
    inicializa(N, K);
    int end_flag = 1, n_loops = 0;

    do {
        end_flag = update_samples();
        
        if (end_flag) {
            
            n_loops++;
            update_centroids();
            
        }

    } while (end_flag);
    

    

    printf("\n#> n_loops: %d\n", n_loops);
    for (int i = 0; i < K; i++) {

        print_cluster(CLUSTERS[i]);
        printf("\n");
    }

    printf("\n#> done!\n\n");
    return 0;
}

/* usage example:



    clock_t begin1 = clock();
    printf("\n#> %.5f pow(97,2)\n", pow(97, 2));
    clock_t end1 = clock();

    clock_t begin2 = clock();
    printf("#> %.5f  pow(97,2)\n", square(97));
    clock_t end2 = clock();

    double time_spent1 = (double)(end1 - begin1) / CLOCKS_PER_SEC;
    double time_spent2 = (double)(end2 - begin2) / CLOCKS_PER_SEC;

    printf("\n#> time1:%.5f \t time2:%.5f\n\n", time_spent1, time_spent2);

    --------------------------

    struct spoint p1 = {.x = 1.0f, .y = 1.0f};
    print_spoint(p1);

    SPoint sp1 = NULL;
    sp1 = initSPoint(sp1, 2.0f, 2.0f);

    print_SPoint(sp1);

    float distance = 0.0f;
    distance = euclidianDistance(p1, *(sp1));
    printf("\n#> distance = %.5f\n", distance);

*/