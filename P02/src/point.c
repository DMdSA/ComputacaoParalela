#include "../include/point.h"
#include "math.h"
#include <stdlib.h>
#include <stdio.h>


void print_point(struct point a) {

    printf("Point:{x:%.3f, y:%.3f,\tK=%d};", a.x, a.y, a.k);
}
