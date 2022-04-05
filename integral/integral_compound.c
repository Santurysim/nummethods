#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libintegral.h"

double f1(double x);
double f2(double x);
double f3(double x);

double f1(double x)
{
    return cos(100.0 * x);
}

double f2(double x)
{
    return exp(-1000.0 * x);
}

double f3(double x)
{
    return 1.0 / sqrt(1.0 - x * x);
}

int main(int argc, char **argv)
{
    size_t N;
    double int_g, int_s;

    if (argc != 2)
        return 1;

    if (sscanf(argv[1], "%zu", &N) != 1)
        return 1;

    int_g = gauss_integral_n(0.0, M_PI, f1, N);
    int_s = simpson_integral_n(0.0, M_PI, f1, N);
    printf("cos(100x)\t%e\t%e\t%e\t%e\n", int_g, fabs(int_g - 0.0),
                                          int_s, fabs(int_s - 0.0));

    int_g = gauss_integral_n(0.0, 1.0, f2, N);
    int_s = simpson_integral_n(0.0, 1.0, f2, N);
    printf("exp(-1000x)\t%e\t%e\t%e\t%e\n", int_g, fabs(int_g - 1e-3),
                                            int_s, fabs(int_s - 1e-3));

    int_g = gauss_integral_n(-1.0, 1.0, f3, N);
    int_s = simpson_integral_n(-1.0, 1.0, f3, N);
    printf("1/sqrt(1-x^2)\t%e\t%e\t%e\t%e\n", int_g, fabs(int_g - M_PI),
                                              int_s, fabs(int_s - M_PI));

    return 0;
}
