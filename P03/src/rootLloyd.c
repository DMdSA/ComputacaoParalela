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

#include "../include/point.h"
#include "../include/cluster.h"

struct point* RANDOM_SAMPLE;
struct cluster* CLUSTERS;

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


int update_samples(struct point* SAMPLE_SUBSET, float* changes, int number_points_per_proccess, int n_clusters ){

    int change_flag = 0;
    
    for(int i = 0; i < number_points_per_proccess; i++){
            
            struct point p = SAMPLE_SUBSET[i];
            float minDist = FLT_MAX;
            int minK = 0, lastK = 0;
            minK = lastK = p.cluster;

            struct cluster c = {0};
            // 0 1 2 3, se o ultimo kluster for 1, começa no 2 e vai até 4 (4%4 = 0)
            for(int d = 0; d < n_clusters; d++){
                
                c = CLUSTERS[d];
                float dist = distanceComparator(p,c); //Calculo da distância
                
                //Se encontrar uma distância inferior à registada, substitui-a
                if(minDist > dist ){
                    
                    minDist = dist;
                    minK = d;
                }
            }

            if (minK != lastK){ //Se o cluster mais proximo for diferente do atual cluster do ponto
                
                //--changes[lastK * 3];  //Alterações na dimensão

                SAMPLE_SUBSET[i].cluster = minK;    //Atualizar os valores do ponto
                change_flag = 1;
            }
            // {d, x, y}
            ++changes[minK * 3];
            changes[minK * 3 + 1] += p.x;   //Somar o X e o Y 
            changes[minK * 3 + 2] += p.y; 
        }

    return change_flag;
}


void first_update_samples(struct point* SAMPLE_SUBSET, float* changes, int number_points_per_proccess, int n_clusters ){

    //First cicle because of the -1 in the points' default klusters
    for(int i = 0; i < number_points_per_proccess; i++){
            
        struct point p = SAMPLE_SUBSET[i];
        float minDist = FLT_MAX;
        int minK = 0, lastK = 0;
        minK = lastK = p.cluster;
    
        struct cluster c = {0};
        // 0 1 2 3, se o ultimo kluster for 1, começa no 2 e vai até 4 (4%4 = 0)
        for(int d = 0; d < n_clusters; d++){
            
            c = CLUSTERS[d];
            float dist = distanceComparator(p, c); //Calculo da distância
            
            //Se encontrar uma distância inferior à registada, substitui-a
            if(minDist > dist){
                
                minDist = dist;
                minK = d;
            }
        }

        if(minK != lastK){ //Se o cluster mais proximo for diferente do atual cluster do ponto
        
            ++changes[minK * 3];
            SAMPLE_SUBSET[i].cluster = minK;    //Atualizar os valores do ponto
        }
        // {d, x, yf}
        changes[minK * 3 + 1] += p.x;   //Somar o X e o Y 
        changes[minK * 3 + 2] += p.y; 
    }
}


void updateCentroids(float* changes, int klusters){
        
    for (int i = 0; i < klusters; i++) {
        
        // (CLUSTERS[i].dimension) =
        int dimension = (CLUSTERS[i].dimension) = (int) changes[i*3];
        //printf("\n#> %f %f %f\n", changes[i*3], changes[i*3+1], changes[i*3+2]);
        (CLUSTERS[i]).x = changes[i*3 + 1] / (float)dimension;
        (CLUSTERS[i]).y = changes[i*3 + 2] / (float)dimension;
    }
}

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
    // inicializar MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // estudar a melhor alternativa
    int number_points_per_proccess = n_points / size; 
    
    // commit struct dataypes for mpi
    /* MPI_POINT */
    MPI_Datatype mpi_point = createMpiPoint();
    /* MPI_CLUSTER */
    MPI_Datatype mpi_cluster = createMpiCluster();
    
    // atribuir ranks
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


    // {d0, d1, d2, d3,  x0, x1, x2, x3,  y0, y1, y2, y3}
    // {d0, x0, y0,   d1, x1, y1,    d2, x2, y1,   d3, x3, y3}
    const int arraysize = n_clusters*3;
    float changesSend[arraysize], changesRecv[arraysize];
    
    for(int i = 0; i < arraysize; i++) {
        changesSend[i] = 0.0f;
    }

    if (rank == 0) {
        for (int i = 0; i < arraysize; i++){
            changesRecv[i] = 0.0f;
        }
    }
    
    //First cicle because of the -1 in the points' default klusters
    first_update_samples(SAMPLE_SUBSET, changesSend, number_points_per_proccess, n_clusters);

    // atualizar a dimensão de cada cluster, após a primeira iteração
    float previousDimensions[n_clusters];
    for(int i = 0; i < n_clusters; i++) {
        previousDimensions[i] = CLUSTERS[i].dimension = (int) changesSend[i*3];
    }
    
    //MPI_Allreduce(MPI_IN_PLACE, changes, arraysize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    // apenas o rank 0 recebe o resultado do reduce
    MPI_Reduce(changesSend, changesRecv, arraysize, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    // todos atualizam os centroids

    if (rank == 0) {
        // atualizar os centroids
        updateCentroids(changesRecv,n_clusters);
    }

    // partilhar os novos clusters com todos os processos
    // obs: a dimensão correta de cada processo está guardada em "previousDimensions"!
    MPI_Bcast(CLUSTERS, n_clusters, mpi_cluster, 0, MPI_COMM_WORLD);

    // recuperar as dimensões corretas para cada processo
    for(int i = 0; i < n_clusters; i++) {
        CLUSTERS[i].dimension = previousDimensions[i];
    }

    int change_flag = 1, iterations = 1;
    
    while(change_flag) {  //Enquanto houver alterações nos clusters

        // reiniciar o array de changes
        for(int i = 0; i < arraysize; i++) {
            changesSend[i] = 0.0f;
        }
        change_flag = update_samples(SAMPLE_SUBSET, changesSend, number_points_per_proccess, n_clusters);
        
        MPI_Allreduce(MPI_IN_PLACE, &change_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (change_flag != 0){

            for(int i = 0; i < n_clusters; i++)
                previousDimensions[i] = CLUSTERS[i].dimension = (int) changesSend[i*3];
        
            MPI_Reduce(changesSend, changesRecv, arraysize, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

            if (rank == 0)
                updateCentroids(changesRecv,n_clusters);
            
            MPI_Bcast(CLUSTERS, n_clusters, mpi_cluster, 0, MPI_COMM_WORLD);

            iterations++;
        }
    }
    
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