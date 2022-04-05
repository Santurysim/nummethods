#include <stddef.h>
#include <math.h>

#include "libintegral.h"

// Simple integral formulas

double gauss_integral(double a, double b, double (*f)(double))
{
    double xm, x0, xp;
    xm = (a + b) * 0.5 - (b - a) * 0.5 * sqrt(0.6);
    x0 = (a + b) * 0.5;
    xp = (a + b) * 0.5 + (b - a) * 0.5 * sqrt(0.6);
    return (b - a) * (5.0 * f(xm) + 8.0 * f(x0) + 5.0 * f(xp)) / 18.0;
}

double simpson_integral(double a, double b, double (*f)(double))
{
    return (b - a) * (f(a) + 4.0 * f((a + b) * 0.5) + f(b)) / 6.0;
}

// Compound integral formulas

double gauss_integral_n(double a, double b, double (*f)(double), size_t N)
{
    double result, xm, x0, xp, right, left;
    result = 0.0;
    for (size_t i = 0; i < N; i++) {
        left  = a + (b - a) * ((double)i / (double)N);
        right = a + (b - a) * ((double)(i + 1) / (double)N);
        xm = (left + right) * 0.5 - (right - left) * 0.5 * sqrt(0.6);
        x0 = (left + right) * 0.5;
        xp = (left + right) * 0.5 + (right - left) * 0.5 * sqrt(0.6);
        result += ((right - left) / 18.0) * (5.0 * f(xm) + 8.0 * f(x0)
                                             + 5.0 * f(xp));
    }
    return result;
}

double simpson_integral_n(double a, double b, double (*f)(double), size_t N)
{
    double result, right, left;
    result = 0.0;
    for (size_t i = 0; i < N; i++) {
        left  = a + (b - a) * ((double)i / (double)N);
        right = a + (b - a) * ((double)(i + 1) / (double)N);
        result +=  ((right - left) / 6.0) * (f(left)
                                             + 4.0 * f((left + right) * 0.5)
                                             + f(right));
    }
    return result;
}

