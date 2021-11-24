#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "schemes.h"

#define MESH_COUNT_1 10
#define MESH_COUNT_2 100
#define MESH_COUNT_3 1000
#define MESH_COUNT_6 1000000

#define MESH_COUNT_MAX MESH_COUNT_6

// FIXME

double get_approximation_error(double *approximation, double *reference, int N);

void estimate_schemes(double A, double *approximation_cache,
	double *reference_cache, double errors_table[4][6]);

void estimate_schemes_with_mesh_count(double A, int N,
	double *approximation_cache, double *reference_cache,
	double *errors_buffer);

const double PARAMETERS[] = {1.0, 10.0, 1000.0};
const int APPROXIMATION_ORDERS[] = {1, 1, 2, 2, 2, 2};

int main(void) {
	double *approximation_cache, *reference_cache;
	double errors[4][6];

	approximation_cache = (double*)malloc((MESH_COUNT_MAX + 1) *
		sizeof(double));
	if (!approximation_cache) {
		return 1;
	}

	reference_cache = (double*)malloc((MESH_COUNT_MAX + 1) * sizeof(double));
	if (!reference_cache) {
		free(approximation_cache);
		return 1;
	}

	/* tbl prologue */
	puts(".TS");
	puts("center, box, allbox;");
	puts("c c c c c c c.");

	puts("No\tE1\tE2\tE3\tE6\tm\tA");
	
	/* tbl data */
	for (int i = 0; i < 3; i++) {
		estimate_schemes(PARAMETERS[i], approximation_cache, reference_cache,
			errors);
		for (int j = 0; j < 6; j++) {
			printf("%d\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%d\t%.2lf\n", j + 1,
				errors[0][j], errors[1][j], errors[2][j], errors[3][j],
				APPROXIMATION_ORDERS[j], PARAMETERS[i]);
		}
	}
	
	/* tbl epilgue */
	puts(".TE");

	free(approximation_cache);
	free(reference_cache);
	return 0;
}

void estimate_schemes_with_mesh_count(double A, int N,
	double *approximation_cache, double *reference_cache,
	double *errors_buffer) {
	compute_reference(reference_cache, N, A);

	compute_scheme1(approximation_cache, N, A);
	errors_buffer[0] = get_approximation_error(approximation_cache,
							reference_cache, N);

	compute_scheme2(approximation_cache, N, A);
	errors_buffer[1] = get_approximation_error(approximation_cache,
							reference_cache, N);

	compute_scheme3(approximation_cache, N, A);
	errors_buffer[2] = get_approximation_error(approximation_cache,
							reference_cache, N);

	compute_scheme4(approximation_cache, N, A);
	errors_buffer[3] = get_approximation_error(approximation_cache,
							reference_cache, N);

	compute_scheme5(approximation_cache, N, A);
	errors_buffer[4] = get_approximation_error(approximation_cache,
							reference_cache, N);

	compute_scheme6(approximation_cache, N, A);
	errors_buffer[5] = get_approximation_error(approximation_cache,
							reference_cache, N);
}

void estimate_schemes(double A, double *approximation_cache,
	double *reference_cache, double errors_table[4][6]) {

	estimate_schemes_with_mesh_count(A, MESH_COUNT_1, approximation_cache,
		reference_cache, errors_table[0]);
	estimate_schemes_with_mesh_count(A, MESH_COUNT_2, approximation_cache,
		reference_cache, errors_table[1]);
	estimate_schemes_with_mesh_count(A, MESH_COUNT_3, approximation_cache,
		reference_cache, errors_table[2]);
	estimate_schemes_with_mesh_count(A, MESH_COUNT_6, approximation_cache,
		reference_cache, errors_table[3]);
}

double get_approximation_error(double *approximation, double *reference,
							   int N) {
	double t;
	double result = 0.0;
	for (int i = 0; i <= N; i++) {
		t = fabs(reference[i] - approximation[i]);
		if (t > result) {
			result = t;
		}
	}
	return result;
}

