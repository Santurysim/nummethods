#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "matrixlib.h"

typedef double (*function_t)(double, unsigned int);

void usage(void);

double ipow(double, unsigned int);
double f(double, unsigned int);
double runge_func(double, unsigned int);
double myabs(double, unsigned int);

double polynomial(double*, size_t, double);
double lagrange_basis_polynomial(double*, size_t, size_t, double);
double lagrange_polynomial(double*, double*, size_t, double);

void solve_and_print_summary(double*, function_t, size_t, double*, double*,
                             double*);

void correctness(size_t n, double*, double*, double*, double*);
void instability(size_t n, double*, double*, double*, double*);
void runge_convergence(size_t n, double*, double*, double*, double*);
void abs_convergence(size_t n, double*, double*, double*, double*);

void generate_uniform_mesh(double*, size_t, double, double);
void generate_chebyshev_mesh(double*, size_t, double, double);

void usage(void) {
    fputs("Usage: ./a.out task_number mesh_size\n", stderr);
}

// functions for approximations

double ipow(double x, unsigned int n) {
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

double f(double x, unsigned int n) {
    double temp = ipow(x, n - 2);
    return x * temp + 2.0 * temp + 3;
}

double runge_func(double x, unsigned int n) {
    (void)n;
    return 1.0 / (25.0 * x * x + 1.0);
}

double myabs(double x, unsigned int n){
    (void)n;
    return fabs(x);
}

// Explicit interpolation polynomial. Requires linear system to be solved
double polynomial(double *a, size_t n, double x) {
    double result = a[n - 1];
    for(size_t i = n - 1; i > 0; i--) {
        result = result * x + a[i - 1];
    }
    return result;
}

double lagrange_basis_polynomial(double *mesh, size_t n, size_t i, double x) {
    double result = 1.0;
    for(size_t j = 0; j < n; j++) {
        if(j == i) continue;
        result *= (x - mesh[j]) / (mesh[i] - mesh[j]);
    }
    return result;
}

// Interpolation polynomial in Lagrange form
double lagrange_polynomial(double *mesh, double *y, size_t n, double x) {
    double result = 0.0;
    for(size_t i = 0; i < n; i++)
        result += y[i] * lagrange_basis_polynomial(mesh, n, i, x);
    return result;
}

void solve_and_print_summary(double *mesh, function_t func, size_t n,
                             double *matrix, double *coeffs, double *vals) {
    double h;
    // Solve in explicit way
    // Write matrix in trasposed form
    for(size_t j = 0; j < n; j++)
        matrix[COORD(0, j, n)] = 1.0;

    for(size_t i = 1; i < n; i++)
        for(size_t j = 0; j < n; j++)
            matrix[COORD(i, j, n)] = matrix[COORD(i - 1, j, n)] * mesh[j];

    for(size_t i = 0; i < n; i++)
        vals[i] = coeffs[i] = func(mesh[i], n);

    // solve system when respective function is ready
    solve_system(matrix, coeffs, n);

    // time to print out results
    // TODO make a groff table?
    for(size_t i = 0; i < n - 1; i++) {
        h = (mesh[i + 1] - mesh[i]) / 3.0;
        printf("%e\t%e\t%e\t%e\n", mesh[i], func(mesh[i], n),
               polynomial(coeffs, n, mesh[i]),
               lagrange_polynomial(mesh, vals, n, mesh[i]));
        printf("%e\t%e\t%e\t%e\n", mesh[i] + h,  func(mesh[i] + h, n),
               polynomial(coeffs, n, mesh[i] + h),
               lagrange_polynomial(mesh, vals, n, mesh[i] + h));
        printf("%e\t%e\t%e\t%e\n", mesh[i] + 2 * h, func(mesh[i] + 2 * h, n),
               polynomial(coeffs, n, mesh[i] + 2 * h),
               lagrange_polynomial(mesh, vals, n, mesh[i] + 2 * h));
    }
    printf("%e\t%e\t%e\t%e\n", mesh[n - 1], func(mesh[n - 1], n),
           polynomial(coeffs, n, mesh[n - 1]),
           lagrange_polynomial(mesh, vals, n, mesh[n - 1]));

}

// mesh generators
void generate_uniform_mesh(double *mesh, size_t n, double left, double right) {
    for(size_t i = 0; i < n; i++)
        mesh[i] = left + ((right - left) * i) / n;
}

void generate_chebyshev_mesh(double *mesh, size_t n, double left,
                             double right) {
    // TODO
}

void correctness(size_t n, double* matrix, double* mesh, double* vals,
                 double *coeffs) {
    generate_uniform_mesh(mesh, n, -1.0, 1.0);
    solve_and_print_summary(mesh, f, n, matrix, coeffs, vals);

    printf("\n\n");

    generate_uniform_mesh(mesh, n, -1.0, 2.0);
    solve_and_print_summary(mesh, f, n, matrix, coeffs, vals);
}

void instability(size_t n, double* matrix, double* mesh, double* vals,
                 double *coeffs) {
    generate_uniform_mesh(mesh, n, -1.0, 2.0);
    solve_and_print_summary(mesh, f, n, matrix, coeffs, vals);
}

void runge_convergence(size_t n, double* matrix, double* mesh, double* vals,
                       double *coeffs) {
    generate_uniform_mesh(mesh, n, -1.0, 1.0);
    solve_and_print_summary(mesh, runge_func, n, matrix, coeffs, vals);

    printf("\n\n");

    generate_chebyshev_mesh(mesh, n, -1.0, 1.0);
    solve_and_print_summary(mesh, runge_func, n, matrix, coeffs, vals);
}

void abs_convergence(size_t n, double* matrix, double* mesh, double* vals,
                     double *coeffs) {
    generate_uniform_mesh(mesh, n, -1.0, 2.0);
    solve_and_print_summary(mesh, myabs, n, matrix, coeffs, vals); 
}

int main(int argc, char **argv) {
    int task_number;
    size_t n;
    double *matrix, *mesh, *vals, *coeffs;

    if(argc != 3) {
        usage();
        return 1;
    }

    if(sscanf(argv[1], "%d", &task_number) != 1) {
        usage();
        return 1;
    }

    if(sscanf(argv[2], "%zu", &n) != 1) {
        usage();
        return 1;
    }

    matrix = (double*)malloc(n * n * sizeof(double));
    if(!matrix) {
        perror("Failed to allocate memory");
        return 2;
    }

    mesh = (double*)malloc(n * sizeof(double));
    if(!mesh) { 
        perror("Failed to allocate memory");
        free(matrix);
        return 2;
    }

    vals = (double*)malloc(n * sizeof(double));
    if(!vals) { 
        perror("Failed to allocate memory");
        free(matrix);
        free(mesh);
        return 2;
    }

    coeffs = (double*)malloc(n * sizeof(double));
    if(!coeffs) { 
        perror("Failed to allocate memory");
        free(matrix);
        free(mesh);
        free(vals);
        return 2;
    }

    switch(task_number) {
        case 1:
            correctness(n, matrix, mesh, vals, coeffs);
            break;
        case 2:
            instability(n, matrix, mesh, vals, coeffs);
            break;
        case 3:
            runge_convergence(n, matrix, mesh, vals, coeffs);
            break;
        case 4:
            abs_convergence(n, matrix, mesh, vals, coeffs);
            break;
        default:
            usage();
            return 1;
    }
    return 0;

}
