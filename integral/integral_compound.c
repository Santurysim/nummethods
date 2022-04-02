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

int main(void)
{
    return 0;
}
