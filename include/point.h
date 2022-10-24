#ifndef POINT_HEADER_FILE
#define POINT_HEADER_FILE

extern const int SPSIZE;

struct __attribute__ ((__packed__, __aligned__(8))) spoint {
    // struct of a space point

    float x;            // x value
    float y;            // y value

    int k;            // cluster
};

typedef struct spoint *SPoint;

float sqrt1(const float x);
float squared(float f);
int minFloat(float a, float b);

SPoint initSPoint(SPoint sp, float x, float y);

float spoint_euclidianDistance(struct spoint a, struct spoint b);
float SPoint_euclidianDistance(SPoint a, SPoint b);

int spoint_tostring(struct spoint a, char* s);
void print_spoint(struct spoint a);
void print_SPoint(SPoint a);

#endif