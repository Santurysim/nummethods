#include <math.h>

#include "schemes.h"

void compute_scheme1(double *arr, int N, double A) {
	double h, coeff;
	h = 1.0 / N;
	coeff = 1.0 - A * h;
	arr[0] = 1.0;
	for (int i = 1; i <= N; i++) {
		arr[i] = arr[i - 1] * coeff;
		if (fabs(arr[i]) <= HAZARDOUS_LOWER_BOUND) {
			arr[i] = 0.0;
		}
	}
}

void compute_scheme2(double *arr, int N, double A) {
	double h, coeff;
	h = 1.0 / N;
	coeff = 1 / (1.0 + A * h);
	arr[0] = 1.0;
	for (int i = 1; i <= N; i++) {
		arr[i] = arr[i - 1] * coeff;
		if (fabs(arr[i]) <= HAZARDOUS_LOWER_BOUND) {
			arr[i] = 0.0;
		}
	}
}

void compute_scheme3(double *arr, int N, double A) {
	double h, coeff;
	h = 1.0 / N;
	coeff = ((2.0 - A * h) / (2.0 + A * h));
	arr[0] = 1.0;
	for (int i = 1; i <= N; i++) {
		arr[i] = arr[i - 1] * coeff;
		if (fabs(arr[i]) <= HAZARDOUS_LOWER_BOUND) {
			arr[i] = 0.0;
		}
	}
}

void compute_scheme4(double *arr, int N, double A) {
	double h, coeff1; /* coeff2 left out */
	int is_going_to_overflow = 0;
	h = 1.0 / N;
	coeff1 = -2.0 * A * h;
	/* coeff2 = 1.0; */
	arr[0] = 1.0;
	arr[1] = 1.0 - A * h;
	for (int i = 2; i <= N; i++) {
		if (is_going_to_overflow) {
			arr[i] = INFINITY;
		} else {
			arr[i] = arr[i - 1] * coeff1 + arr[i - 2] /* * coeff2 */;
			if (fabs(arr[i]) <= HAZARDOUS_LOWER_BOUND) {
				arr[i] = 0.0;
			}
			is_going_to_overflow = (fabs(arr[i]) >= HAZARDOUS_UPPER_BOUND);
		}
	}
}

void compute_scheme5(double *arr, int N, double A) {
	double h, coeff1, coeff2;
	int is_going_to_overflow = 0;
	h = 1.0 / N;
	coeff1 = 2.0 / (1.5 + A * h);
	coeff2 = - 0.5 / (1.5 + A * h);
	arr[0] = 1.0;
	arr[1] = 1.0 - A * h;
	for (int i = 2; i <= N; i++) {
		if (is_going_to_overflow) {
			arr[i] = INFINITY;
		} else {
			arr[i] = arr[i - 1] * coeff1 + arr[i - 2] * coeff2;
			if (fabs(arr[i]) <= HAZARDOUS_LOWER_BOUND) {
				arr[i] = 0.0;
			}
			is_going_to_overflow = (fabs(arr[i]) >= HAZARDOUS_UPPER_BOUND);
		}
	}
}

void compute_scheme6(double *arr, int N, double A) {
	double h, coeff1, coeff2;
	int is_going_to_overflow = 0;
	h = 1.0 / N;
	coeff1 = 4.0;
	coeff2 = 2.0 * A * h - 3.0;
	arr[0] = 1.0;
	arr[1] = 1.0 - A * h;
	for (int i = 2; i <= N; i++) {
		if (is_going_to_overflow) {
			arr[i] = INFINITY;
		} else {
			arr[i] = arr[i - 1] * coeff1 + arr[i - 2] * coeff2;
			if (fabs(arr[i]) <= HAZARDOUS_LOWER_BOUND) {
				arr[i] = 0.0;
			}
			is_going_to_overflow = (fabs(arr[i]) >= HAZARDOUS_UPPER_BOUND);
		}
	}
}

void compute_reference(double *arr, int N, double A) {
	for (int i = 0; i <= N; i++) {
		if ((-A * ((double)i / (double)N)) <= HAZARDOUS_LOWER_BOUND_LOG) {
			arr[i] = 0.0;
		} else {
			arr[i] = exp(-A * ((double)i / (double)N));
		}
	}
}
