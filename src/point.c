#include "../include/point.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>

const int SPSIZE = 100; 

void print_spoint(struct spoint a) {

    printf("Point:{x:%.3f, y:%.3f,\tK=%d};", a.x, a.y, a.k);
}
