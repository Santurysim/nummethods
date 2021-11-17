#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "schemes.h"

#define MESH_COUNT_1 10
#define MESH_COUNT_2 100
#define MESH_COUNT_3 1000
#define MESH_COUNT_6 1000000

#define MESH_COUNT_MAX MESH_COUNT_6

double get_approximation_error(double *approximation, double *reference, int N);

void estimate_scheme(scheme_t compute_scheme, int A,
					 double *approximation_cache, double *reference_cache,
					 double *errors_buffer);

int main(void) {
	double *approximation_cache, reference_cache;
	double errors[4];
	puts(".TS");

}

void estimate_scheme(scheme_t compute_scheme, double A,
					 double *approximation_cache, double *reference_cache,
					 double *errors_buffer) {
}

double get_approximation_error(double *approximation, double *reference,
							   int N) {
	double t;
	double result = 0.0;
	for (int i = 0; i <= N; i++) {
		t = fabs(reference[i * (MESH_COUNT_MAX / N)] - approximation[i]);
		if (t > result) {
			result = t;
		}
	}
	return result;
}

