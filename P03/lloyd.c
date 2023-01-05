/**
 * @file initialization.c
 * @author dma
 * @brief testing lloyd's initialization technique for random points and their clusters
 * @version 0.1
 * @date 2022-12-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <stddef.h>
#include <mpi.h>
#include <assert.h>

/**
 * @brief struct of a point in space
 * 
 */
struct point {

    float x, y, cluster_dist;          // (x, y) values on space
    int cluster;               // current associated cluster
} ; 
typedef struct point *Point;


/**
 * @brief cluster struct
 * 
 */
struct cluster {

    float x, y;             // (x, y) from space
    int dimension;          // number of points associated
};

typedef struct cluster *Cluster;



struct point* RANDOM_SAMPLE;
struct cluster* CLUSTERS;

/**
 * @brief initialize points and clusters' arrays with necessary memory
 * 
 * @param n_points 
 * @param n_clusters 
 * @return int 0, if successfully allocated
 */
int initialize(int n_points, int n_clusters) {

    // clusters must be initialized for all threads..
    //CLUSTERS = (struct cluster*) malloc(sizeof(struct cluster) * n_clusters);
    RANDOM_SAMPLE = (struct point*) malloc(sizeof(struct point) * n_points);
    
    if (!RANDOM_SAMPLE) {
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
    
    // random points + first k clusters
    int i = 0;
    for (; i < n_clusters; i++) {

        float x = (float) rand() / RAND_MAX, y = (float) rand() / RAND_MAX;
        RANDOM_SAMPLE[i] = (struct point) {.x = x, .y = y, .cluster_dist = FLT_MAX, .cluster = -1};
        CLUSTERS[i] = (struct cluster) {.x = x, .y = y, .dimension = 0};
    }
    // remaining points
    for (int j = i; j < n_points; j++) {

        float x = (float) rand() / RAND_MAX, y = (float) rand() / RAND_MAX;
        RANDOM_SAMPLE[j] = (struct point) {.x = x, .y = y, .cluster_dist = FLT_MAX, .cluster = -1};
    }
    
}

/**
 * @brief gets the distance between two points, without square rooting
 * 
 * @param a point a
 * @param b point b
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


void showPoint(struct point p, int* rank){

    if (rank != NULL)
        printf("\n#> Process[%d] - P(%.2f,%.2f,%.2f,%d)", *(rank), p.x, p.y, p.cluster_dist, p.cluster);
    else
        printf("\n#> P(%.2f,%.2f,%.2f,%d)", p.x, p.y, p.cluster_dist, p.cluster);
}

void showCluster(struct cluster c, int* rank) {
    
    if (rank != NULL)
        printf("\n#> Process[%d] - C(%.2f, %.2f, %d)", *rank, c.x, c.y, c.dimension);
    else
        printf("\n#> C(%.2f, %.2f, %d)", c.x, c.y, c.dimension);

}

int update_samples(struct point* SAMPLE_SUBSET, float* changes, int number_points_per_proccess, int n_clusters ){

    struct point p = {0};
    struct cluster c = {0};
    float minDist, dist; 
    int minK, lastK;
    int change_flag = 1;
    for(int i = 0; i < number_points_per_proccess; i++){
            p = SAMPLE_SUBSET[i];
            minDist = p.cluster_dist;
            minK = lastK = p.cluster;

            // 0 1 2 3, se o ultimo kluster for 1, começa no 2 e vai até 4 (4%4 = 0)
            for(int d = lastK+1;d < lastK + n_clusters; d++){
                c = CLUSTERS[d%n_clusters];
                dist = distanceComparator(p,c); //Calculo da distância
                //Se encontrar uma distância inferior à registada, substitui-a
                if(minDist > dist ){
                    minDist = dist;
                    minK = d%n_clusters;
                }
            }

            if(minK != lastK){ //Se o cluster mais proximo for diferente do atual cluster do ponto
                changes[lastK]--;  //Alterações na dimensão
                changes[minK]++;

                SAMPLE_SUBSET[i].cluster = minK;    //Atualizar os valores do ponto
                SAMPLE_SUBSET[i].cluster_dist = minDist;
               change_flag = 1;
            }
            changes[minK+n_clusters] += p.x;   //Somar o X e o Y 
            changes[minK+n_clusters+1] += p.y; 

        }
    return change_flag;
}

void updateCentroids(float* CENTER_MEANS, int klusters){
        int index = 0, i = 0;
        for (i = 0; i < klusters; i++) {

        index = i+i;
        CLUSTERS[i].dimension = CENTER_MEANS[i];
        int dimension = CENTER_MEANS[i];
        (CLUSTERS[i]).x = CENTER_MEANS[i+klusters] / (float)dimension;
        (CLUSTERS[i]).y = CENTER_MEANS[i+(klusters*2)] / (float)dimension;
    }
}

int main(int argc, char* argv[]) {

    // controll main arguments
    main_arg_control(argc, argv);
    
    int n_points = 0, n_clusters = 0;
    sscanf(argv[1], "%d", &n_points);
    sscanf(argv[2], "%d", &n_clusters);

    printf("\n#> %d %d\n", n_points, n_clusters);
    const int arraysize = n_clusters*2;
    float CENTR_MEANS[arraysize];
    float changes[n_clusters];

    if (n_clusters > n_points) {

        fprintf(stderr, "\n#> Fatal: you can't have more clusters than points.\n");
        exit(3);
    }

    // ---------------------------------------------------------

    int size = 0, rank = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // estudar a melhor alternativa
    int number_points_per_proccess = n_points / size; 
    
    // commit struct dataypes for mpi
    /* MPI_POINT */
    const int p_nitems=4;
    int          p_blocklengths[4] = {1,1,1,1};
    MPI_Datatype p_types[4] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
    MPI_Datatype mpi_point;
    MPI_Aint     p_offsets[4];

    p_offsets[0] = offsetof(struct point, x);
    p_offsets[1] = offsetof(struct point, y);
    p_offsets[2] = offsetof(struct point, cluster_dist);
    p_offsets[3] = offsetof(struct point, cluster);

    MPI_Type_create_struct(p_nitems, p_blocklengths, p_offsets, p_types, &mpi_point);
    MPI_Type_commit(&mpi_point);
    
    /* MPI_CLUSTER */
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


    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // clusters must be initialized for all processes
    CLUSTERS = (struct cluster*) malloc(sizeof(struct cluster) * n_clusters);
    if (rank == 0) {
        
        // only "master"/root process will initialize and populate the
        // random numbers sample and the clusters

        initialize(n_points, n_clusters);
        populate(n_points, n_clusters);
        
    }

    /**
     * for each process, create a buffer that will hold a subset/chunk of the entire
     * points sample
    */

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

    struct point p = {0};
    struct cluster c = {0};
    float minDist, dist; 
    int minK, lastK;
    int change_flag = 1;  

     for(int i = 0; i < arraysize; i++){  //Reiniciar o array com os centros
            CENTR_MEANS[i] = 0.0f;
        }
        for(int i = 0; i < n_clusters*3; i++){
            changes[i] = (float) 0;
        }
    //First cicle because of the -1 in the points' default klusters
    for(int i = 0; i < number_points_per_proccess; i++){
            p = SAMPLE_SUBSET[i];
            minDist = FLT_MAX;
            minK = lastK = p.cluster;

            // 0 1 2 3, se o ultimo kluster for 1, começa no 2 e vai até 4 (4%4 = 0)
            for(int d = lastK+1;d < lastK + n_clusters; d++){
                c = CLUSTERS[d%n_clusters];
                dist = distanceComparator(p,c); //Calculo da distância
                //Se encontrar uma distância inferior à registada, substitui-a
                if(minDist > dist ){
                    minDist = dist;
                    minK = d%n_clusters;
                }
            }

            if(minK != lastK){ //Se o cluster mais proximo for diferente do atual cluster do ponto
                changes[minK]++;

                SAMPLE_SUBSET[i].cluster = minK;    //Atualizar os valores do ponto
                SAMPLE_SUBSET[i].cluster_dist = minDist;
               change_flag = 1;
            }
            changes[minK+n_clusters] += p.x;   //Somar o X e o Y 
            changes[minK+n_clusters+1] += p.y; 

        }
        /*
    MPI_Gather(SAMPLE_SUBSET, number_points_per_proccess, mpi_point,
                        RANDOM_SAMPLE,number_points_per_proccess,
                        mpi_point,0,MPI_COMM_WORLD);*/
    for(int i = 0; i < n_clusters; i++){
            CLUSTERS[i].dimension = (int) changes[i];
        }
        
    MPI_Allreduce(MPI_IN_PLACE,changes,n_clusters,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    
    if(rank==0){
        updateCentroids(changes,n_clusters);
    }


    while(change_flag){  //Enquanto houverem alterações nos clusters

        for(int i = 0; i < arraysize; i++){  //Reiniciar o array com os centros
            CENTR_MEANS[i] = 0.0f;
        }
        for(int i = 0; i < n_clusters*3; i++){
            changes[i] = (float) CLUSTERS[i].dimension;
        }
        change_flag = 0;

        update_samples(SAMPLE_SUBSET, changes, number_points_per_proccess, n_clusters);

        for(int i = 0; i < n_clusters; i++){
            CLUSTERS[i].dimension = (int) changes[i];
        }
        
        MPI_Allreduce(MPI_IN_PLACE,changes,n_clusters,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

        if(rank==0){
            updateCentroids(changes,n_clusters);
        }
    }
    
    for (int i = 0; i < number_points_per_proccess; i++) {
        showPoint(SAMPLE_SUBSET[i], &rank);
    }

    MPI_Gather(SAMPLE_SUBSET, number_points_per_proccess, mpi_point,
                        RANDOM_SAMPLE,number_points_per_proccess,
                        mpi_point,0,MPI_COMM_WORLD);


    if (rank == 0) {
        free(RANDOM_SAMPLE);
    }
    free(CLUSTERS);
    MPI_Finalize();

    return 0;
}


/*
void sumCentroids(float* CENTER_MEANS, const int samples, const int klusters){
    const int arraysize = klusters*2;

    for(int i = 0; i < arraysize; i++){
        CENTER_MEANS[i] = 0.0f;
    }
    struct point p = {0};

    int index = 0, i = 0;
    for(i = 0; i < samples; i++) {

        p = RANDOM_SAMPLE[i];
        index = p.cluster + p.cluster;
        CENTER_MEANS[index] += p.x;
        CENTER_MEANS[index+1] += p.y;
    }

}*/

    /**
    printf("\n#> n_clusters from rank[%d] = %d\n", rank, n_clusters);
    for (int i = 0; i < n_clusters; i++) {

        showCluster(CLUSTERS[i], &rank);
    }
    */

            /*
            MPI_Allreduce(
                    void* send_data,
                    void* recv_data,
                    int count,
                    MPI_Datatype datatype,
                    MPI_Op op,
                    MPI_Comm communicator)*/

    /**
     * Gather all partial data down to the root process
     * MPI_Gather( void* send_data, int send_count,
                    MPI_Datatype send_datatype,
                    void* recv_data,
                    int recv_count,
                    MPI_Datatype recv_datatype,
                    int root,
                    MPI_Comm communicator)
    */

/*

    -Inicio:                                                                            - done
        - scatter dos pontos
        - bcast dos clusters
    
    - Funcao Inicial (subset_pontos, clusters):                                         - @todo
        for (i=0, i < nclusters; i++)
            eachpoint : mindist

    - Funcao lidar clusters (clusters, &arrayDistanciasAuxiliar):                       - @todo
        - ReduceAll (clusters.dimensao e arrayDistancias) (sum)
        - aqui, ja todos vao ter toda a info
        - divisão e atualizar clusters
    
    - Funcao loop "otimizado" (subset_pontos, clusters) [int::flag]:                    - @todo
        for (optimized):
            min optimized

    - Gather                                                                            - @todo
        - resposta


*/



/*

for (int i = (.cluster+1) % nclusters; i != .cluster; i = (i+1)%4) {
}

for (int i = .cluster+1; i < .cluster+nclusters; i++) {
}

.cluster = 0,  : 1, 2, 3
.cluster = 1,  : 2, 3, 4%4=0

.cluster = -1, : 0, 1, 2, 3 

.cluster = 0,  : 1, 2, 3
.cluster = 1,  : 2, 3, 0
.cluster = 2,  : 3, 0, 1
.cluster = 3,  : 0, 1, 2
*/
