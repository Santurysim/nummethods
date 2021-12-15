#include <math.h>

#include <assert.h>

double f1(double, double*, int);
double f2(double, double*, int);
double f3(double, double*, int);

double ref1(double);
double ref2(double);
double ref3(double);

double f1(double x, double *y, int n) {
	assert(n >= 3);
	(void)x;
	return 3 * y[0] - y[1] + y[2];
}

double f2(double x, double *y, int n) {
	assert(n >= 3);
	(void)x;
	return y[0] + y[1] + y[2];
}

double f3(double x, double *y, int n) {
	assert(n >= 3);
	(void)x;
	return 4 * y[0] - y[1] + 4 * y[2];
}

double ref1(double x) {
	return exp(x) + exp(2 * x) + exp(5 * x);
}

double ref2(double x) {
	return exp(x) - 2 * exp(2 * x) + exp(5 * x);
}

double ref3(double x) {
	return -exp(x) - 3 * exp(2 * x) + 3 * exp(5 * x);
}
