#include <stdio.h>
#include <math.h>

#include "libintegral_2d.h"

double f(double x, double y);

double f(double x, double y)
{
    return exp(x) * sin(y);
}

int main(int argc, char **argv)
{
    size_t nx, ny;
    double ans;
    point_t top = POINT_ZERO;
    point_t bottom = POINT_ZERO;

    if (argc != 7)
        return 1;

    if (sscanf(argv[1], "%zu", &nx) != 1)
        return 2;

    if (sscanf(argv[2], "%zu", &ny) != 1)
        return 2;

    if (sscanf(argv[3], "%lf", &bottom.x) != 1)
        return 2;

    if (sscanf(argv[4], "%lf", &bottom.y) != 1)
        return 2;

    if (sscanf(argv[5], "%lf", &top.x) != 1)
        return 2;

    if (sscanf(argv[6], "%lf", &top.y) != 1)
        return 2;

    ans = integral_2d(f, bottom, top, nx, ny);
    printf("%e\n%e\n", ans, fabs((M_E - 1.0) * (1.0 - cos(1.0)) - ans));
    return 0;
}
