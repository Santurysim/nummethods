#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "common.h"
#include "matrixlib.h"

void f(double *x, double *y, size_t n);

void f1(double *x, double *y, size_t n);

void generate_jacobian_t(void(*func)(double*, double*, size_t), double *x,
                         double *matrix, double *temp, size_t n, double delta);

void find_root(void (*func)(double*, double*, size_t), double *x0,
                            double *matrix, double *x, size_t n);

double norm2(double const *x, size_t n);

void generate_jacobian_t(void(*func)(double*, double*, size_t), double *x,
                         double *matrix, double *temp, size_t n, double delta)
{
    double t;
    for (size_t i = 0; i < n; i++) {
        t = x[i];
        x[i] += delta;
        func(x, temp, n);
        memmove(matrix + i * n, temp, n * sizeof(double));
        x[i] = t - delta;
        func(x, temp, n);
        x[i] = t;
        for (size_t j = 0; j < n; j++) {
            matrix[COORDT(j, i, n)] -= temp[j];
            matrix[COORDT(j, i, n)] /= (2.0  * delta);
        }
    }
}

double norm2(double const *x, size_t n)
{
    double result = 0.0;
    for (size_t i = 0; i < n; i++) {
        result += SQUARE(x[i]);
    }
    return sqrt(result);
}

void find_root(void (*func)(double*, double*, size_t), double *x0,
               double *matrix, double *x, size_t n)
{
    do {
        generate_jacobian_t(func, x0, matrix, x, n, 1e-4); // XXX
        func(x0, x, n);
        for (size_t i = 0; i < n; i++) x[i] = -x[i];
        solve_system(matrix, x, n);
        for (size_t i = 0; i < n; i++) x0[i] += x[i];
        func(x0, x, n);
    } while (norm2(x, n) > EPS);
}

// Sample functions

void f(double *x, double *y, size_t n)
{
    assert(n == 10);
    y[0] = sin(x[0]);
    y[1] = sin(x[1]);
    y[2] = sin(x[2]);
    y[3] = sin(x[3]);
    y[4] = sin(x[4]);
    y[5] = sin(x[5]);
    y[6] = sin(x[6]);
    y[7] = sin(x[7]);
    y[8] = sin(x[8]);
    y[9] = sin(x[9]);

}

void f1(double *x, double *y, size_t n)
{
    assert(n == 10);
    y[0] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9]
           + x[0] * x[0] - 11.0;
    y[1] = x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9]
           + x[1] * x[1] - 10.0;
    y[2] = x[2] + x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[2] * x[2]
           - 9.0;
    y[3] = x[3] + x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[3] * x[3] - 8.0;
    y[4] = x[4] + x[5] + x[6] + x[7] + x[8] + x[9] + x[4] * x[4] - 7.0;
    y[5] = x[5] + x[6] + x[7] + x[8] + x[9] + x[5] * x[5] - 6.0;
    y[6] = x[6] + x[7] + x[8] + x[9] + x[6] * x[6] - 5.0;
    y[7] = x[7] + x[8] + x[9] + x[7] * x[7] - 4.0;
    y[8] = x[8] + x[9] + x[8] * x[8] - 3.0;
    y[9] = x[9] + x[9] * x[9] - 2.0;
}

int main(int argc, char **argv)
{
    size_t n;
    double *matrix, *x, *x0;
    FILE *in;

    if (argc != 2) return 2;
    
    n = 10;

    matrix = (double*)malloc(n * n * sizeof(double));
    if (!matrix)
        return 2;

    x = (double*)malloc(n * sizeof(double));
    if (!x)
        return 2;

    x0 = (double*)malloc(n * sizeof(double));
    if (!x0)
        return 2;

    in = fopen(argv[1], "r");
    if (!in)
        return 3;

    for (size_t i = 0; i < n; i++)
        if ((fscanf(in, "%lf", x0 + i)) != 1)
            return 4;

    find_root(f1, x0, matrix, x, n);

    for (size_t i = 0; i < n; i++)
        printf("%lf ", x0[i]);

    printf("\n");

    return 0;
}
