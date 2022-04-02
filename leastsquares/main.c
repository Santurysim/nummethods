#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "common.h"
#include "matrixlib.h"

typedef double (*function_t)(double);
typedef double (*basis_t)(double, size_t, size_t);

double ipow(double x, size_t n);
double basis(double x, size_t n, size_t N);

double approximation(double x, double *solution, basis_t phi, size_t n,
                     size_t N);

void generate_matrix(double *matrix, double *x, size_t n, size_t m,
                     basis_t phi);

void solve_and_print_summary(double *matrix, double *x, double *y, size_t *map,
                             double *solution, size_t n, size_t N,
                             function_t target, basis_t phi);

double F(double);

double ipow(double x, size_t n)
{
    double b = 1.0, c = x;
    while (n > 0) {
        if(n % 2) {
            b *= c;
            n--;
        }
        else {
            c *= c;
            n /= 2;
        }
    }
    return b;
}

double basis(double x, size_t n, size_t N)
{
    (void)N;
    return ipow(x, n);
}

double sin_basis(double x, size_t n, size_t N)
{
    (void)N;
    return sin(M_PI * (double)(n + 1) * x);
}

double F(double x)
{
    return x * sin(M_PI * x);
}

double runge(double x)
{
    return 1.0 / (1.0 + 25.0 * x * x);
}

void generate_matrix(double *matrix, double *x, size_t n, size_t N,
                     basis_t phi)
{
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < N; j++)
            matrix[COORD(i, j, N)] = phi(x[j], i, N);
}

double approximation(double x, double *solution, basis_t phi, size_t n,
                     size_t N)
{
    double result = 0.0;
    for (size_t i = 0; i < n; i++)
        result += solution[i] * phi(x, i, N);
    return result;
}

void solve_and_print_summary(double *matrix, double *x, double *y, size_t *map,
                             double *solution, size_t n, size_t N,
                             function_t target, basis_t phi)
{
    double t, h;
    generate_matrix(matrix, x, n, N, phi);
    for (size_t i = 0; i < N; i++) {
        y[i] = target(x[i]);
    }

    solve_system(matrix, solution, y, map, N, n);

    for (size_t i = 0; i < N - 1; i++) {
        h = (x[i + 1] - x[i]) / 3;
        t = approximation(x[i], solution, phi, n, N);
        printf("%e\t%e\t%e\t%e\n", x[i], target(x[i]), t,
               fabs(target(x[i]) - t));

        t = approximation(x[i] + h, solution, phi, n, N);
        printf("%e\t%e\t%e\t%e\n", x[i] + h, target(x[i] + h), t,
               fabs(target(x[i] + h) - t));
 
        t = approximation(x[i] + 2 * h, solution, phi, n, N);
        printf("%e\t%e\t%e\t%e\n", x[i] + 2 * h, target(x[i] + 2 * h), t,
               fabs(target(x[i] + 2 * h) - t));
    }
}

int main(int argc, char **argv)
{
    int result;
    size_t n, N, *map;
    double *x, *y, *matrix, *solution, t1;
    FILE *in;

    if (argc != 3)
        return 1;

    in = fopen(argv[1], "r");
    if (!in) {
        perror("fatal: fopen");
        return 1;
    }

    if (sscanf(argv[2], "%zu", &n) != 1)
        return 1;

    N = 0;
    while ((result = fscanf(in, "%lf", &t1) != EOF)) {
        if (result != 1) {
            fprintf(stderr, "fatal: malformed data\n");
            fclose(in);
            return 1;
        }
        N++;
    }

    fclose(in);

    matrix = (double*) calloc(n * N, sizeof(double));
    if (!matrix) {
        perror("fatal: malloc");
        return 2;
    }

    x = (double*)malloc(N * sizeof(double));
    if (!x) {
        perror("fatal: malloc");
        return 2;
    }

    y = (double*)malloc(N * sizeof(double));
    if (!y) {
        perror("fatal: malloc");
        return 2;
    }

    map = (size_t*)malloc(n * sizeof(size_t));
    if (!map) {
        perror("fatal: malloc");
        return 2;
    }

    solution = (double*)malloc(n * sizeof(size_t));
    if (!solution) {
        perror("fatal: malloc");
        return 2;
    }

    in = fopen(argv[1], "r");
    if (!in) {
        perror("fatal: fopen");
        return 1;
    }

    for (size_t i = 0; i < N; i++) {
        result = fscanf(in, "%lf", x + i);
        if (result != 1) {
            fprintf(stderr, "fatal: malformed data\n");
            return 1;
        }
    }

    fclose(in);

    solve_and_print_summary(matrix, x, y, map, solution, n, N, runge, basis);

    return 0;
}
