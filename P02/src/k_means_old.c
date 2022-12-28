#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../include/point.h"
#include "../include/cluster.h"

#include <time.h>
#include <math.h>

#define N 10000000
#define K 4


struct spoint* RANDOM_SAMPLE;
struct cluster* CLUSTERS;
float CENTR_MEANS[K*2];

/**
 * @brief inicializa N amostras, no espaço, com valores random, e inicializa K clusters com random centroids
 * 
 * @param n_points número de amostras a utilizar
 * @param n_clusters número de clusters
 */
void inicializa(const int n_points, const int n_clusters) {

    // malloc for K clusters
    CLUSTERS = (struct cluster*) malloc(sizeof(struct cluster) * n_clusters);
    
    // malloc for vector of size N
    RANDOM_SAMPLE = (struct spoint*) malloc(sizeof(struct spoint) * n_points);

    // if malloc fails
    if (!RANDOM_SAMPLE || !CLUSTERS) {
        fprintf(stderr, "Fatal: failed to allocate bytes.\n");
        exit(1);
    }

    int i = 0;
    float x = 0.0f, y = 0.0f;
    // rand seed
    srand(10);
    // random sample generator
    
    do {
        x = (float) rand()/RAND_MAX;
        y = (float) rand()/RAND_MAX;
        RANDOM_SAMPLE[i++] = (struct spoint) {.x = x, .y = y, .k = -1};
    } while (i < n_points);

    // initialize each of K clusters and assign their first centroid
    i = 0;
    do {
        x = RANDOM_SAMPLE[i].x;
        y = RANDOM_SAMPLE[i].y;
        CLUSTERS[i] = (struct cluster) {.x =x, .y = y, .dimension = 0};
        i++;
    } while (i < n_clusters);
}
 

int update_samples(const int samples, const int klusters) {

    //Auxiliar Variables
    struct spoint sp = {0};
    struct cluster* restrict c = NULL;
    float dist = 0.0f, minDist = FLT_MAX;
    int lastK = 0, minK = 0, changes = 0;

    float x = 0.0f, y = 0.0f;

    for (int p = 0; p < samples; p++) {

        sp = RANDOM_SAMPLE[p];
        x = sp.x; y = sp.y;
        minK = lastK = sp.k;
        minDist = FLT_MAX;
        
        for (int k = 0; k < klusters; k++) {

            c = &(CLUSTERS[k]);
            dist = sqrtf((x - c->x) * (x - c->x) + (y - c->y) * (y - c->y));
            if (dist < minDist) {
                
                minDist = dist;
                minK = k;
            }
        }

        if (minK != lastK) {
            if (lastK != -1) ((CLUSTERS[lastK]).dimension)--;
            (CLUSTERS[minK].dimension)++;
            RANDOM_SAMPLE[p].k = minK;
            if (!changes) changes = 1;
        }
    }
    return changes;
}



void update_centroids(const int samples, const int klusters) {

    const int arraysize = klusters * 2;
    for (int i = 0; i < arraysize; i++) {

        *(CENTR_MEANS+i) = 0.0f;
    }

    int index = 0, dimension = 0;
    struct spoint p = RANDOM_SAMPLE[0];

    for(int i = 0; i < samples; i++) {

        p = *(RANDOM_SAMPLE+i);
        index = p.k * 2;
        CENTR_MEANS[index] += p.x;
        CENTR_MEANS[index+1] += p.y;
    }

    for (int i = 0; i < klusters; i++) {

        index = i+i;
        dimension = CLUSTERS[i].dimension;
        (CLUSTERS[i]).x = CENTR_MEANS[index++] / (float)dimension;
        (CLUSTERS[i]).y = CENTR_MEANS[index] / (float)dimension;
    }

}



int main() {

    // initialize random samples + clusters
    int end_flag = 1, n_loops = 0;
    const int n=N, k=K;
 
    inicializa(n, k);

    do {
        end_flag = update_samples(n,k);
        
        if (end_flag) {
            
            n_loops++;
            update_centroids(n,k);    
        }
    
    } while (end_flag);
    
    
    printf("N = %d, K = %d",N,K);
    for (int i = 0; i < K; i++) {
    
       printf("\nCenter: (%.3f, %.3f) : Size: %d", CLUSTERS[i].x, CLUSTERS[i].y, CLUSTERS[i].dimension);
    }
    printf("\nIterations: %d\n", n_loops);
    fflush(stdout);
    free(CLUSTERS);
    free(RANDOM_SAMPLE);
    
    return 0;
}

