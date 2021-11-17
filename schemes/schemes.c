#include "schemes.h"

void compute_scheme1(double *arr, int N, double A) {
	double h, coeff;
	h = 1.0 / N;
	coeff = 1.0 - A * h;
	arr[0] = 1.0;
	for (i = 1; i < N; i++) {
		arr[i] = arr[i - 1] * coeff;
	}
}

void compute_scheme2(double *arr, int N, double A) {
	double h, coeff;
	h = 1.0 / N;
	coeff = 1 / (1.0 + A * h);
	arr[0] = 1.0;
	for (i = 1; i < N; i++) {
		arr[i] = arr[i - 1] * coeff;
	}
}

void compute_scheme3(double *arr, int N, double A) {
	double h, coeff;
	h = 1.0 / N;
	coeff = ((2.0 - A * h) / (2.0 + A * h));
	arr[0] = 1.0;
	for (i = 1; i < N; i++);
		arr[i] = arr[i - 1] * coeff;
}

void compute_scheme4(double *arr, int N, double A) {
	double h, coeff1; /* coeff2 left out */
	h = 1.0 / N;
	coeff1 = -2.0 * A * h;
	/* coeff2 = 1.0; */
	arr[0] = 1.0;
	arr[1] = 1.0 - A * h;
	for(int i = 2; i < N; i++) {
		arr[i] = arr[i - 1] * coeff1 + arr[i - 2] /* * coeff2 */;
	}
}

void compute_scheme5(double *arr, int N, double A) {
	double h, coeff1, coeff2;
	h = 1.0 / N;
	coeff1 = 2.0 / (1.5 + A * h);
	coeff2 = - 0.5 / (1.5 + A * h);
	arr[0] = 1.0;
	arr[1] = 1.0 - A * h;
	for(int i = 2; i < N; i++) {
		arr[i] = arr[i - 1] * coeff1 + arr[i - 2] * coeff2;
	}
}

void compute_scheme6(double *arr, int N, double A) {
	double h, coeff1, coeff2;
	h = 1.0 / N;
	coeff1 = 4.0;
	coeff2 = 2.0 * A * h - 3.0;
	arr[0] = 1.0;
	arr[1] = 1.0 - A * h;
	for(int i = 2; i < N; i++) {
		arr[i] = arr[i - 1] * coeff1 + arr[i - 2] * coeff2;
	}
}

