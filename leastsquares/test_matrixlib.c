#include <stdio.h>

#include "common.h"
#include "matrixlib.h"

int main(void)
{
    double matrix[] = {
         1,  0,  0,  0,  0,
        -1,  1,  0,  0,  0,
        -1, -1,  1,  0,  0,
        -1, -1, -1,  1,  0,
        -1, -1, -1, -1,  1
    };
    double x[] = {0, 0, 0, 0, 0};
    double y[] = {0, 0, 0, 0, 0};
    size_t map[5];
    solve_system(matrix, x, y, map, 5, 5);
    return 0;
}
