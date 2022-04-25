#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <nath.h>

#include "matrixlib.h"

void f(double *x, double *y, size_t n);

void jacobian(void(*func)(double*, double*, size_t), double *x, double *matrix,
                          double *temp, size_t n, double delta);

void find_root(void (*func)(double*, double*, size_t));

void jacobian(void(*func)(double*, double*, size_t), double *x, double *matrix,
                          double *temp, size_t n, double delta)
{
    double t;
    for (size_t i = 0; i < n; i++)
    {
        t = x[i];
        x[i] += delta;
        func(x, temp, n);
        memmove(matrix + i * n, temp, n * sizeof(double));
        x[i] = t;
        func(x, temp, n);
        for (size_t j = 0; j < n; j++) {
            matrix[COORDT(j, i, order)] -= temp[j];
            matrix[COORDT(j, i, order)] /= delta;
        }
    }
}




int main(int argc, char **argv)
{
    return 0;
}
