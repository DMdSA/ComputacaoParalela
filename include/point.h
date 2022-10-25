#ifndef POINT_HEADER_FILE
#define POINT_HEADER_FILE

struct spoint {
    // struct of a space point

    float x, y, dist;   // x, y and distance to it's cluster
    int k;              // cluster
} ; 
typedef struct spoint *SPoint;

void print_spoint(struct spoint a);

//__attribute__ ((packed, aligned(32)))
#endif