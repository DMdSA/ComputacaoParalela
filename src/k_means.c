#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>

#include "../include/point.h"
#include "../include/cluster.h"

#include <time.h>
#include <math.h>
#include <malloc.h>

#define N 10000000
#define K 4


struct spoint* RANDOM_SAMPLE;
struct cluster* CLUSTERS; 
int i = 0, j = 0, means_size = K*2;
float MEANS_ARRAY[K*2];

float euclidian_distance(float x1, float y1, float x2, float y2) {

    return sqrt((x2-x1) * (x2-x1) + (y2-y1) * (y2-y1));
}

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
    for (i = 0; i < N; i++) {

        RANDOM_SAMPLE[i] = (struct spoint) {.x = (float) rand()*(1.0f/RAND_MAX), .y = (float) rand()*(1.0f/RAND_MAX), .k = -1};
    }

    // initialize each of K clusters and assign their first centroid
    for (i = 0; i < K; i++) {

        CLUSTERS[i] = (struct cluster) {.dimension = 0, .x = RANDOM_SAMPLE[i].x, .y = RANDOM_SAMPLE[i].y};
    }
}
 

/**
 * @brief goes through all samples and updates their cluster
 * 
 * @return int 0 if there was no changes were made
 */
int update_samples() {

    // auxiliar variables
    float minDist = FLT_MAX, dist = 0.0f, x=0.0f, y=0.0f;
    int minK = INT_MAX, changes_flag = 0;
    struct spoint p = {0}, c= {0};


    // for each of the samples
    #pragma omp simd
    for (i = 0; i < N; i++) {

        minDist = INT_MAX;
        // get current point
        p = *(RANDOM_SAMPLE+i);
        minK = p.k;
        // default values for minimum calculation

        // for each of the *other* clusters
        for (j = 0; j < K; j++) {
            // calculate the euclidian distance
            //dist = euclidianDistance(p, (CLUSTERS[j]).centroid);
            dist = euclidian_distance(p.x, p.y, CLUSTERS[j].x, CLUSTERS[j].y);
            
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
    int index = 0, dimension = 0;
    for (i = 0; i < means_size; i++) {
        MEANS_ARRAY[i] = 0.0f;
    }

    struct spoint p;

    // for each of the samples
    for (i = 0; i < N; i++) {

        p=RANDOM_SAMPLE[i];
        index = p.k * 2;
        MEANS_ARRAY[index] += p.x;
        MEANS_ARRAY[index+1] += p.y;

        // k = 0, 0 1
        // k = 1, 2 3
        // k = 2, 4 5
    }

    // for each cluster, calculate the new centroid
    for (i = 0; i < K; i++) {
        
        // aux variables
        index = i*2;
        dimension = (CLUSTERS[i]).dimension;

        // mean calculation
        //printf("\n#> dimension %d, index %d, %.5f, %.5f", dimension, index, (CLUSTERS[i].centroid).x, (CLUSTERS[i].centroid).y);
        (CLUSTERS[i].x) = MEANS_ARRAY[index] * 1.0f/(float)dimension; 
        (CLUSTERS[i].y) = MEANS_ARRAY[index + 1] * 1.0f/(float)dimension;

    }

}



int main() {

    // debug purposes
    //printf("It's %ld bit system /n", sizeof(void *) * 8); 
    //                4           4            8
    //printf("\n#> int(%ld), float(%ld), double(%ld)\n\n", sizeof(int), sizeof(float), sizeof(double));
    //char spoint_string[SPSIZE];
    //char cluster_string[CSIZE];
    


    // initialize random samples + clusters
    int end_flag = 1, n_loops = 0;
    clock_t begin = clock();
    inicializa(N, K);

    do {
        end_flag = update_samples();
        
        if (end_flag) {
            
            n_loops++;
            update_centroids();
            
        }

    } while (end_flag);
    
    
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\n#> time:%.5f\n", time_spent);
    
    printf("\n#> n_loops: %d\n", n_loops);
    for (int i = 0; i < K; i++) {

        print_cluster(CLUSTERS[i]);
        printf("\n");
    }

    free(RANDOM_SAMPLE);
    free(CLUSTERS);
    printf("\n#> done!\n\n");
    
    return 0;
}
