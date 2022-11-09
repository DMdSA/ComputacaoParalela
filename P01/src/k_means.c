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

/*
float euclidian_distance(float x1, float y1, float x2, float y2) {

    return sqrt((x2-x1) * (x2-x1) + (y2-y1) * (y2-y1));
}
*/

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

    // debug purposes
    //printf("It's %ld bit system /n", sizeof(void *) * 8); 
    //                4           4            8
    //printf("\n#> int(%ld), float(%ld), double(%ld)\n\n", sizeof(int), sizeof(float), sizeof(double));
    //char spoint_string[SPSIZE];
    //char cluster_string[CSIZE];
    

    //struct spoint ex[2] = {0};
    //ex[0] = (struct spoint) {.x = 1.0f, .y = 0.0f, .k = 1};
    //ex[1] = (struct spoint) {.x = 4.0f, .y = 2.0f, .k = 3};
//
    //printf("\n#> [0]: x:%p, y:%p, k:%p", &ex[0].x, &ex[0].y, &ex[0].k);
    //printf("\n#> [1]: x:%p, y:%p, k:%p", &ex[1].x, &ex[1].y, &ex[1].k);
    //float exx[6] = {0};
    //exx[0] = 1.0f;exx[1] = 1.0f;exx[2] = 1.0f;exx[3] = 1.0f;exx[4] = 1.0f;exx[5] = 1.0f;
//
    //printf("\n#> struct size: %ld, float size: %ld", sizeof(ex), sizeof(exx));

    // initialize random samples + clusters
    int end_flag = 1, n_loops = 0;
    const int n=N, k=K;
    //clock_t begin = clock();
    inicializa(n, k);

    do {
        end_flag = update_samples(n,k);
        
        if (end_flag) {
            
            n_loops++;
            update_centroids(n,k);    
        }
    
    } while (end_flag);
    
    
    //clock_t end = clock();
    //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("\n#> time:%.5f\n", time_spent);
    
    printf("\n#> n_loops: %d\n", n_loops);
    for (int i = 0; i < K; i++) {
    
       printf("\n#> Center: (%.3f, %.3f), Size: %d", CLUSTERS[i].x, CLUSTERS[i].y, CLUSTERS[i].dimension);
    }
    fflush(stdout);
    free(CLUSTERS);
    free(RANDOM_SAMPLE);
    printf("\n#> done!\n\n");
    
    return 0;
}

