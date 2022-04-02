#ifndef LIBINTEGRAL_H
#define LIBINTEGRAL_H

typedef double (*function_t)(double);

double gauss_integral(double a, double b, double (*f)(double));
double simpson_integral(double a, double b, function_t f);

double gauss_integral_n(double a, double b, double (*f)(double), size_t N);
double simpson_integral_n(double a, double b, double (*f)(double), size_t N);

#endif /* LIBINTEGRAL_H */

