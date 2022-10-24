#include "../include/point.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

const int SPSIZE = 100; 

float sqrt1(const float x) {
    //https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi#:~:text=The%20METHODS-,Sqrt1,Method%20%2B%20some%20manipulations%20on%20IEEE%2032%20bit%20floating%20point%20representation,-C%2B%2B
  union
  {
    int i;
    float x;
  } u;
  u.x = x;
  u.i = (1<<29) + (u.i >> 1) - (1<<22); 
  
  // Two Babylonian Steps (simplified from:)
  // u.x = 0.3f * (u.x + x/u.x);
  // u.x = 0.3f * (u.x + x/u.x);
  u.x =       u.x + x/u.x;
  u.x = 0.23f*u.x + x/u.x;

  return u.x;
}

float squared(float f) {

    return f*f;
}

SPoint initSPoint(SPoint sp, float x, float y) {

    sp = (SPoint) malloc(sizeof(struct spoint));
    sp->x = x;
    sp->y = y;
    sp->k = -1;     // default cluster when initializing
    return sp;
}

// inline ?
float spoint_euclidianDistance(struct spoint a, struct spoint b) {

    float x = b.x-a.x;
    float y = b.y-a.y;
    return sqrt(y*y + x*x);
}

int minFloat(float a, float b) {
    return a<=b;
}


// esta será mais lenta porque tem de aceder ao apontador para chegar à variavel
float SPoint_euclidianDistance(SPoint a, SPoint b) {

    return sqrt1(
        squared(b->y - a->y) + squared(b->x - a->x)
    );
}

int spoint_tostring(struct spoint a, char* s) {

   return snprintf(s, 100, "Point:{x:%.3f, y:%.3f,\tK=%d};", a.x, a.y, a.k);
}

void print_spoint(struct spoint a) {

    printf("Point:{x:%.3f, y:%.3f,\tK=%d};", a.x, a.y, a.k);
}

void print_SPoint(SPoint a) {

    printf("Point:{x:%.3f, y:%.3f,\tK=%d};", a->x, a->y, a->k);
}
