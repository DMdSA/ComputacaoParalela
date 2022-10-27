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
    // rand seed
    srand(10);
    // random sample generator
    
    do {
        RANDOM_SAMPLE[i] = (struct spoint) {.x = (float)rand()/RAND_MAX, .y = (float)rand()/RAND_MAX, .dist = FLT_MAX, .k = -1};
        i++;
    } while (i < n_points);

    // initialize each of K clusters and assign their first centroid
    i = 0;
    do {
        CLUSTERS[i] = (struct cluster) {.x = RANDOM_SAMPLE[i].x, .y = RANDOM_SAMPLE[i].y, .dimension = 0};
        i++;
    } while (i < n_clusters);
}
 

int update_samples(const int samples, const int klusters) {

    struct spoint* restrict sp = NULL;
    struct cluster* restrict c = NULL;
    float dist = 0.0f, sqred = 0.0f, minDist = FLT_MAX;
    int lastK = 0, minK = 0, changes = 0;

    for (int p = 0; p < samples; p++) {

        sp = &(RANDOM_SAMPLE[p]);
        lastK = sp->k;
        minK = lastK;
        minDist = FLT_MAX;

        for (int k = 0; k < klusters; k++) {

            c = &(CLUSTERS[k]);
            sqred = (sp->x - c->x) * (sp->x - c->x) + (sp->y - c->y) * (sp->y - c->y);
            dist = sqrtf(sqred);

            if (dist < minDist) {
                
                minDist = dist;
                //sp->k = k;
                minK = k;
            }
        }

        if (minK != lastK) {
            if (lastK != -1) ((CLUSTERS[lastK]).dimension)--;
            (CLUSTERS[minK].dimension)++;
            sp->dist = minDist;
            sp->k = minK;
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

    #pragma omp simd
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
    


    // initialize random samples + clusters
    int end_flag = 1, n_loops = 0;
    const int n=N, k=K;
    clock_t begin = clock();
    inicializa(n, k);

    do {
        end_flag = update_samples(n,k);
        
        if (end_flag) {
            
            n_loops++;
            update_centroids(n,k);    
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

    free(CLUSTERS);
    free(RANDOM_SAMPLE);
    printf("\n#> done!\n\n");
    
    return 0;
}
