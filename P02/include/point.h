#ifndef POINT_HEADER_FILE
#define POINT_HEADER_FILE

/**
 * @brief struct of a point in space
 * 
 */
struct point {

    float x, y;          // (x, y) values on space
    int k;               // current associated cluster
} ; 
typedef struct point *Point;

/**
 * @brief default print of a point info
 * 
 * @param a 
 */
void print_point(struct point a);

//__attribute__ ((packed, aligned(32)))
#endif