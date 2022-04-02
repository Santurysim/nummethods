#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libintegral.h"

#define POWERFUNC_DECL(n) double pow##n (double x);

#define POWERFUNC(n) \
double pow##n (double x) \
{ \
    return ipow(x, n); \
} \

double ipow_derivative_max(double a, double b, size_t n, size_t k);

double ipow_true_integral(double a, double b, size_t n);

double ipow(double x, size_t n);

POWERFUNC_DECL(0)
POWERFUNC_DECL(1)
POWERFUNC_DECL(2)
POWERFUNC_DECL(3)
POWERFUNC_DECL(5)
POWERFUNC_DECL(9)

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

POWERFUNC(0)
POWERFUNC(1)
POWERFUNC(2)
POWERFUNC(3)
POWERFUNC(5)
POWERFUNC(9)

double ipow_true_integral(double a, double b, size_t n)
{
    return (ipow(b, n + 1) - ipow(a, n + 1)) / (double)(n + 1);
}

double ipow_derivative_max(double a, double b, size_t n, size_t k)
{
    double result = 1.0;
    (void) a;
    if (k > n)
        return 0.0;
    for (size_t i = 0; i < k; i++)
        result *= (double)(n - i);
    result *= ipow(b, n - k);
    return result;
}

int main(int argc, char **argv)
{
    double a, b, integ, resid;
    if (argc != 3)
        return 1;

    if(sscanf(argv[1], "%lf", &a) != 1)
        return 1;
    if(sscanf(argv[2], "%lf", &b) != 1)
        return 1;

    printf("Gauss formula:\n");
    integ = gauss_integral(a, b, pow0);
    resid = fabs(integ - ipow_true_integral(a, b, 0));
    printf("%d\t%e\t%e\t%e\n", 0, integ, resid,
            ipow_derivative_max(a, b, 0, 6) * ipow(b - a, 7) / 2016000.0);
    integ = gauss_integral(a, b, pow1);
    resid = fabs(integ - ipow_true_integral(a, b, 1));
    printf("%d\t%e\t%e\t%e\n", 1, integ, resid,
            ipow_derivative_max(a, b, 1, 6) * ipow(b - a, 7) / 2016000.0);
    integ = gauss_integral(a, b, pow2);
    resid = fabs(integ - ipow_true_integral(a, b, 2));
    printf("%d\t%e\t%e\t%e\n", 2, integ, resid,
            ipow_derivative_max(a, b, 2, 6) * ipow(b - a, 7) / 2016000.0);
    integ = gauss_integral(a, b, pow3);
    resid = fabs(integ - ipow_true_integral(a, b, 3));
    printf("%d\t%e\t%e\t%e\n", 3, integ, resid,
            ipow_derivative_max(a, b, 3, 6) * ipow(b - a, 7) / 2016000.0);
    integ = gauss_integral(a, b, pow5);
    resid = fabs(integ - ipow_true_integral(a, b, 5));
    printf("%d\t%e\t%e\t%e\n", 5, integ, resid,
            ipow_derivative_max(a, b, 5, 6) * ipow(b - a, 7) / 2016000.0);
    integ = gauss_integral(a, b, pow9);
    resid = fabs(integ - ipow_true_integral(a, b, 9));
    printf("%d\t%e\t%e\t%e\n", 9, integ, resid,
            ipow_derivative_max(a, b, 9, 6) * ipow(b - a, 7) / 2016000.0);

    printf("\nSimpson formula:\n");
    integ = simpson_integral(a, b, pow0);
    resid = fabs(integ - ipow_true_integral(a, b, 0));
    printf("%d\t%e\t%e\t%e\n", 0, integ, resid,
            ipow_derivative_max(a, b, 0, 4) * ipow(b - a, 5) / 2880.0);
    integ = simpson_integral(a, b, pow1);
    resid = fabs(integ - ipow_true_integral(a, b, 1));
    printf("%d\t%e\t%e\t%e\n", 1, integ, resid,
            ipow_derivative_max(a, b, 1, 4) * ipow(b - a, 5) / 2880.0);
    integ = simpson_integral(a, b, pow2);
    resid = fabs(integ - ipow_true_integral(a, b, 2));
    printf("%d\t%e\t%e\t%e\n", 2, integ, resid,
            ipow_derivative_max(a, b, 2, 4) * ipow(b - a, 5) / 2880.0);
    integ = simpson_integral(a, b, pow3);
    resid = fabs(integ - ipow_true_integral(a, b, 3));
    printf("%d\t%e\t%e\t%e\n", 3, integ, resid,
            ipow_derivative_max(a, b, 3, 4) * ipow(b - a, 5) / 2880.0);
    integ = simpson_integral(a, b, pow5);
    resid = fabs(integ - ipow_true_integral(a, b, 5));
    printf("%d\t%e\t%e\t%e\n", 5, integ, resid,
            ipow_derivative_max(a, b, 5, 4) * ipow(b - a, 5) / 2880.0);
    integ = simpson_integral(a, b, pow9);
    resid = fabs(integ - ipow_true_integral(a, b, 9));
    printf("%d\t%e\t%e\t%e\n", 9, integ, resid,
            ipow_derivative_max(a, b, 9, 4) * ipow(b - a, 5) / 2880.0);

    return 0;
}
