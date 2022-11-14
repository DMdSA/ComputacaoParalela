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
int initialize(int n_points, int n_clusters, int n_threads) {

    #pragma omp parallel sections num_threads(n_threads) 
    
    {
    
        #pragma omp section 
        {
            CLUSTERS = (struct cluster*) malloc(sizeof(struct cluster) * n_clusters);
        }
        #pragma omp section 
        {
            RANDOM_SAMPLE = (struct point*) malloc(sizeof(struct point) * n_points);
        }
    }
    
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
    #pragma omp for schedule(static)
    for (int i = 0; i < n_points; i++) {
        
        float x = (float) rand() / RAND_MAX, y = (float) rand() / RAND_MAX;
        RANDOM_SAMPLE[i] = (struct point) {.x = x, .y = y, .k = -1};
    }

    // clusters
    #pragma omp for schedule(static)
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


void updateSamples(int samples, int klusters, float* centroid_mean_array) {

    float min_dist = 0.0f;
    int min_index = 0;

    #pragma omp default(shared) private(min_dist, min_index) firstprivate(samples, klusters)
    {

        #pragma omp for schedule (static)

            for (int i = 0; i < samples; i++) {

                struct point point = RANDOM_SAMPLE[i];
                min_dist = euclidian_distance(point, CLUSTERS[0]);
                min_index = 0;

                for (int j = 1; j < klusters; j++) {

                    struct cluster cluster = CLUSTERS[j];

                    float dist = euclidian_distance(point, cluster);

                    if (dist < min_dist) {

                        min_dist = dist;
                        min_index = j;
                    }
                }

                //point.k = min_index;
                CLUSTERS[min_index].dimension++;

                centroid_mean_array[min_index*2] += point.x;
                centroid_mean_array[min_index*2+1] += point.y; 
            }
    }
}


int update_cluster(struct cluster* cluster, int cluster_index, float* centroid_mean_array) {

    float x = (*cluster).x, y = (*cluster).y;
    int index = cluster_index*2;

    (*cluster).x = centroid_mean_array[index] / (float)(*cluster).dimension;
    (*cluster).y = centroid_mean_array[index+1] / (float)(*cluster).dimension;

    // clear the auxiliar values
    centroid_mean_array[index] = 0.0f;
    centroid_mean_array[index+1] = 0.0f;

    // clear cluster dimension for next iteration
    //(*cluster).dimension = 0;

    if (cluster->x == x && cluster->y == y) return 0;
    return 1;
}

int update_clusters(int samples, int klusters, int n_threads, float* centroid_mean_array) {

    int flag = 0;

    #pragma omp parallel for simd num_threads(n_threads) private(klusters) reduction(+:flag) schedule(static)
    for (int i = 0; i < klusters; i++) {
        
        flag += update_cluster(&CLUSTERS[i], i, centroid_mean_array);
        int id = omp_get_thread_num();
        printf("T%d:i%d\n", id, i);
    }
    printf("\n\n");
    return flag;
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

    if (n_clusters > n_points) {

        // k_means 1 1 and 2 2 are aborted...?

        fprintf(stderr, "\n#> Fatal: you can't have more clusters than points.\n");
        exit(3);
    }
    
    int centr_means_size = n_clusters + n_clusters;
    float CENTR_MEANS[centr_means_size];
    
    #pragma omp parallel num_threads(n_threads)
    
    #pragma omp for
    for (int i = 0; i < centr_means_size; i++) CENTR_MEANS[i]=0.0f;

    int end_flag = 1, n_loops = 0;
    float begin = omp_get_wtime();

    initialize(n_points, n_clusters, n_threads);
    populate(n_points, n_clusters);

    
    for (end_flag = 1; (end_flag && n_loops <= 20); n_loops++) {

        for (int i = 0; i < n_clusters; CLUSTERS[i].dimension = 0, i++);
        updateSamples(n_points, n_clusters, CENTR_MEANS);
        end_flag = update_clusters(n_points, n_clusters, n_threads, CENTR_MEANS);
    }
    --n_loops;
    
    float end = omp_get_wtime();
    float time_spent = (end-begin);
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