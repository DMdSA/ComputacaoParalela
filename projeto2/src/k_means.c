#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <float.h>

#include <time.h>

#include "../include/point.h"
#include "../include/cluster.h"


/**
 * @brief array of random points
 * 
 */
struct point* RANDOM_SAMPLE;

/**
 * @brief array of clusters
 * 
 */
struct cluster* CLUSTERS;


/**
 * @brief initialize points and clusters' arrays with necessary memory
 * 
 * @param n_points 
 * @param n_clusters 
 * @return int 0, if successfully allocated
 */
int initialize(const int n_points, const int n_clusters) {

    CLUSTERS = (struct cluster*) malloc(sizeof(struct cluster) * n_clusters);
    RANDOM_SAMPLE = (struct point*) malloc(sizeof(struct point) * n_points);

    if (!RANDOM_SAMPLE || !CLUSTERS) {
        fprintf(stderr, "\n#> Fatal: failed to allocate requested bytes.\n");
        exit(2);
    }
    return 0;
}

/**
 * @brief populate each points and clusters' arrays with random values
 * 
 * @param n_points number of points
 * @param n_clusters number of clusters
 */
void populate(int n_points, int n_clusters) {

    // default rand seed
    srand(10);

    // random points
    for (int i = 0; i < n_points; i++) {
        
        float x = (float) rand() / RAND_MAX, y = (float) rand() / RAND_MAX;
        RANDOM_SAMPLE[i] = (struct point) {.x = x, .y = y, .k = -1};
    }

    // clusters
    for (int i = 0; i < n_clusters; i++) {

        float x = (RANDOM_SAMPLE[i]).x, y = (RANDOM_SAMPLE[i]).y;
        CLUSTERS[i] = (struct cluster) {.x = x, .y = y, .dimension = 0};
    }
}

/**
 * @brief populate with threads
 * 
 * @param n_points 
 * @param n_clusters 
 */
void populate_thread(int n_points, int n_clusters) {

    // default rand seed
    srand(10);

    #pragma omp parallel num_threads(2)
    #pragma omp for
    // random points
    for (int i = 0; i < n_points; i++) {
        
        float x = (float) rand() / RAND_MAX, y = (float) rand() / RAND_MAX;
        RANDOM_SAMPLE[i] = (struct point) {.x = x, .y = y, .k = -1};
    }

    #pragma omp for
    // clusters
    for (int i = 0; i < n_clusters; i++) {

        float x = (RANDOM_SAMPLE[i]).x, y = (RANDOM_SAMPLE[i]).y;
        CLUSTERS[i] = (struct cluster) {.x = x, .y = y, .dimension = 0};
    }
}

float euclidian_distance(struct point point, struct cluster cluster) {

    float x = point.x-cluster.x, y = point.y-cluster.y;
    return (sqrtf(x*x + y*y));
}

float euclidian_distance_conditional(struct point point, struct cluster cluster) {

    float x1 = point.x, y1 = point.y, x2 = cluster.x, y2 = cluster.y;
    
    if (x1 != x2 && y1 != y2) {

        float x = x2-x1, y = y2-y1;
        return (sqrtf(x*x + y*y));
    }

    if (x1 == x2 && y1 == y2) {

        return 0.0f;
    }

    if (x1 == x2) {

        float y = y2-y1;
        return (y >= 0) ? y : -y;
    }

    float x = x2-x1;
    return (x >= 0) ? x : -x;
}

int update_samples(const int samples, const int klusters, const int n_threads) {


    //Auxiliar Variables
    struct point point = {0};
    struct cluster cluster = {0};
    float dist = 0.0f, minDist = FLT_MAX;
    int lastK = 0, minK = 0, changes = 0;

    float x = 0.0f, y = 0.0f;

    
    #pragma omp for schedule(static, 5)
    for (int p = 0; p < samples; p++) {

        point = RANDOM_SAMPLE[p];
        minK = lastK = point.k;
        minDist = FLT_MAX;


        //#pragma omp critical
        for (int k = 0; k < klusters; k++) {

            cluster = CLUSTERS[k];
            dist = euclidian_distance(point, cluster);
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


void update_centroids(float* CENTR_MEANS, const int samples, const int klusters, const int n_threads) {

    /*
        - os valores de schedule fazem variar os tempos de execução
        - cuidado: ao dividir, se o numero for inferior o schedule fica a 0...
    */

    const int arraysize = klusters * 2;
    
    //#pragma omp parallel num_threads(n_threads)
    #pragma omp for schedule (static, arraysize/8)
    for (int i = 0; i < arraysize; i++) {

        CENTR_MEANS[i] = 0.0f;
    }

    struct point p = RANDOM_SAMPLE[0];

    int index = 0, i = 0;
    #pragma omp for private(i, index) schedule(static, samples/8)
    for(i = 0; i < samples; i++) {

        p = RANDOM_SAMPLE[i];
        index = p.k + p.k;
        CENTR_MEANS[index] += p.x;
        CENTR_MEANS[index+1] += p.y;
    }

    // aqui, se houver divisao de schedule, não termina
    #pragma omp for private(i, index) schedule(static,5)
    for (i = 0; i < klusters; i++) {

        index = i+i;
        int dimension = CLUSTERS[i].dimension;
        (CLUSTERS[i]).x = CENTR_MEANS[index++] / (float)dimension;
        (CLUSTERS[i]).y = CENTR_MEANS[index] / (float)dimension;
    }
}



/**
 * @brief checks correct usage of program
 * 
 * @param argc number of arguments
 * @param argv arguments
 */
void main_arg_control(int argc, char* argv[]) {

    if (argc != 4) {

        printf("\n#> [k_means:error]\n\t%s <number of points> <number of clusters> <number of threads>\n\n"
        , argv[0]);
        
        exit(1);
    }
}



int main(int argc, char* argv[]) {

    // controll main arguments
    main_arg_control(argc, argv);
    
    int n_points = 0, n_clusters = 0, n_threads = 0;
    sscanf(argv[1], "%d", &n_points);
    sscanf(argv[2], "%d", &n_clusters);
    sscanf(argv[3], "%d", &n_threads);

    float CENTR_MEANS[n_clusters * 2];

    if (n_clusters > n_points) {

        // k_means 1 1 and 2 2 are aborted...?

        fprintf(stderr, "\n#> Fatal: you can't have more clusters than points.\n");
        exit(3);
    }
    

    int end_flag = 1, n_loops = 0;

    clock_t begin = clock();
    #pragma omp parallel num_threads(n_threads)
    initialize(n_points, n_clusters);
    populate(n_points, n_clusters);

    
    for (end_flag = 1; end_flag; ) {

        end_flag = update_samples(n_points, n_clusters, n_threads);
        
        if (end_flag) {
            
            n_loops++;
            update_centroids(CENTR_MEANS, n_points, n_clusters, n_threads);    
        }
    }
    clock_t end = clock();
    double time_spent = (double) (end-begin) / CLOCKS_PER_SEC;
    printf("\n#> exec: %lf\n\n", time_spent);
    
    printf("N = %d, K = %d, T = %d", n_points, n_clusters, n_threads);
    for (int i = 0; i < n_clusters; i++) {
    
       printf("\nCenter: (%.3f, %.3f) : Size: %d", CLUSTERS[i].x, CLUSTERS[i].y, CLUSTERS[i].dimension);
    }
    printf("\nIterations: %d\n", n_loops);


    /**
     * if n_clusters <= 2, free aborts... (?)
     * 
     */
    free(CLUSTERS);
    free(RANDOM_SAMPLE);
    return 0;
}