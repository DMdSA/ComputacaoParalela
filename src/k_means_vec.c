#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>

#include <time.h>
#include <math.h>
#include <malloc.h>

#define N 10000000
#define K 4

float* RANDOM_SAMPLES __attribute__((aligned (32))), *CLUSTERS __attribute__((aligned (32)));;

float MEANS_ARRAY[K*2];


/**
 * @brief inicializa N amostras, no espaço, com valores random, e inicializa K clusters com random centroids
 * 
 * @param n_points número de amostras a utilizar
 * @param n_clusters número de clusters
 */
void inicializa(int n_points, int n_clusters) {

    // malloc for vector of size N
    RANDOM_SAMPLES = (float *) malloc(sizeof(float) * n_points);

    // malloc for N clusters
    CLUSTERS = (float *) malloc(sizeof(float) * n_clusters);

    
    // if malloc fails
    if (!RANDOM_SAMPLES || !CLUSTERS) {
        fprintf(stderr, "Fatal: failed to allocate bytes.\n");
        exit(1);
    }

    // rand seed
    srand(10);

    // random sample generator
    for (int i = 0; i < n_points;) {

        RANDOM_SAMPLES[i++] = (float) rand()/RAND_MAX;          // point's X
        RANDOM_SAMPLES[i++] = (float) rand()/RAND_MAX;          // point's Y
        RANDOM_SAMPLES[i++] = -1;                               // point's default K
    }

    // initialize each of K clusters and assign their first centroid
    for (int i = 0; i < n_clusters; i+=3) {

        CLUSTERS[i] = RANDOM_SAMPLES[i];                        // cluster's centroid's X
        CLUSTERS[i+1] = RANDOM_SAMPLES[i+1];                    // cluster's centroid's Y
        CLUSTERS[i+2] = 0;                                      // cluster's centroid's dimension
    }
}
 

/**
 * @brief goes through all samples and updates their cluster
 * 
 * @return int 0 if there was no changes were made
 */
int update_samples(const int n_points, const int n_clusters) {

    // auxiliar variables
    float minDist = FLT_MAX, dist = 0.0f;
    int minK = INT_MAX, pk = 0, changes = 0;    
    float px = 0.0f, py = 0.0f, kx = 0.0f, ky = 0.0f;

    // for each of the samples
    for (int i = 0; i < n_points; i+=3) {

        // get current point
        px = RANDOM_SAMPLES[i];                 // 0
        py = RANDOM_SAMPLES[i+1];               // 1
        pk = (int)RANDOM_SAMPLES[i+2];          // 2
        minK = pk;
        
        // default values for minimum calculation
        
        minDist = FLT_MAX;

        // for each of the *other* clusters
        for (int k = 0; k < n_clusters; k+=3) {

            kx = CLUSTERS[k];
            ky = CLUSTERS[k+1];

            dist = sqrtf((px-kx) * (px-kx) + (py-ky) * (py-ky));
            
            if (dist < minDist){
                minDist = dist;
                minK = k/3;
            }
        }

        if (minK != pk) {
            if (pk != -1) (CLUSTERS[(pk*3)+2])--;
            (CLUSTERS[(minK*3)+2])++;
            RANDOM_SAMPLES[i+2] = (float)minK;
            if (!changes) changes = 1;
        }   
    }

    return changes;
}



void update_centroids(const int n_points, const int n_clusters) {


    // initialize means array
    float px = 0.0f, py = 0.0f, dimension = 0.0f;
    int index = 0, means_size = K*2, pk = 0;
    
    for (int i = 0; i < means_size; MEANS_ARRAY[i] = 0.0f, i++);

    // for each of the samples
    for (int i = 0; i < n_points;) {

        px = RANDOM_SAMPLES[i++];
        py = RANDOM_SAMPLES[i++];
        pk = (float)RANDOM_SAMPLES[i++];
        index = pk * 2;
        MEANS_ARRAY[index] += px;
        MEANS_ARRAY[index+1] += py;

        // k = 0, 0 1
        // k = 1, 2 3
        // k = 2, 4 5
    }

    int auxIndex = 0;
    // for each cluster, calculate the new centroid
    for (int i = 0; i < n_clusters; i+=3) {
        
        // aux variables
        index = auxIndex*2;
        dimension = CLUSTERS[i+2];

        // mean calculation
        //printf("\n#> dimension %d, index %d, %.5f, %.5f", dimension, index, (CLUSTERS[i].centroid).x, (CLUSTERS[i].centroid).y);
        (CLUSTERS[i]) = MEANS_ARRAY[index] /dimension; 
        (CLUSTERS[i+1]) = MEANS_ARRAY[index + 1] /dimension;
        auxIndex++;
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
    const int new_n = N*3, new_k = K*3;
    
    clock_t begin = clock();

    inicializa(new_n, new_k);
    do {
        end_flag = update_samples(new_n, new_k);
        
        if (end_flag) {
            
            n_loops++;
            update_centroids(new_n, new_k);
            
        }

    } while (end_flag);
    
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\n#> time:%.5f\n", time_spent);
    
    printf("\n#> n_loops: %d\n", n_loops);
    for (int i = 0; i < new_k; i+=3) {

        printf("\n#> Center: (%.3f, %.3f) : Size : %.0f", CLUSTERS[i], CLUSTERS[i+1], CLUSTERS[i+2]);
    }

    free(RANDOM_SAMPLES);
    free(CLUSTERS);
    printf("\n#> done!\n\n");
    
    return 0;
}
