#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double (*function_t)(double, size_t);

void usage(void) {
    fputs("Usage: ./a.out task_number mesh_size", stderr);
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
        for(int j = 0; j < n; j++)
            matrix[COORD(i, j, n)] = matrix[COORD(i, j, n)] * mesh[j];

    for(size_t i = 0; i < n; i++)
        vals[i] = coeffs[i] = func(mesh[i], n);

    // solve system when respective function is ready

    // time to print out results
    // TODO make a groff table?
    for(size_t i = 0; i < n - 1; i++) {
        h = (mesh[i + 1] - mesh[i]) / 3.0;
        printf("%e\t%e\t%e\n", func(mesh[i], n),
               polynomial(coeffs, n, mesh[i]),
               lagrange_polynomial(mesh, vals, n, mesh[i]));
        printf("%e\t%e\t%e\n", func(mesh[i] + h, n),
               polynomial(coeffs, n, mesh[i] + h),
               lagrange_polynomial(mesh, vals, n, mesh[i] + h));
        printf("%e\t%e\t%e\n", func(mesh[i] + 2 * h, n),
               polynomial(coeffs, n, mesh[i] + 2 * h),
               lagrange_polynomial(mesh, vals, n, mesh[i] + 2 * h));
    }
    printf("%e\t%e\t%e\n", func(mesh[n - 1], n),
           polynomial(coeffs, n, mesh[n - 1]),
           lagrange_polynomial(mesh, vals, n, mesh[n - 1]));

}

void correctness(size_t n) {}

void instability(size_t n) {}

void runge_convergence(size_t n) {}

void abs_convergence(size_t n) {}

int main(int argc, char **argv) {
    int task_number;
    size_t n;

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

    switch(task_number) {
        case 1:
            correctness(n);
            break;
        case 2:
            instability(n);
            break;
        case 3:
            runge_convergence(n);
            break;
        case 4:
            abs_convergence(n);
            break;
        default:
            usage();
            return 1;
    }
    return 0;

}
