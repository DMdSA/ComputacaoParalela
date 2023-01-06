#ifndef POINT_FILE_HEADER
#define POINT_FILE_HEADER

/**
 * @brief struct of a point in space
 * 
 */
struct point {

    float x, y;          // (x, y) values on space
    int cluster;               // current associated cluster
} ; 
typedef struct point *Point;

#endif