#ifndef POINT_HEADER_FILE
#define POINT_HEADER_FILE

extern const int SPSIZE;

struct spoint {
    // struct of a space point

    float x;            // x value
    float y;            // y value

    int k;            // cluster
};

typedef struct spoint *SPoint;

float sqrt1(const float x);
float squared(float f);

SPoint initSPoint(SPoint sp, float x, float y);

float euclidianDistance(struct spoint a, struct spoint b);
float euclidianDistanceP(SPoint a, SPoint b);

int spoint_tostring(struct spoint a, char* s);
void print_spoint(struct spoint a);
void print_SPoint(SPoint a);

#endif