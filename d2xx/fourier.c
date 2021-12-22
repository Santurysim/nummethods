#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SQUARE(x) ((x)*(x))

typedef double (*function_t)(double);

double f(double);

double reference(double);

const double b = -1.0;

double eigenvector(int n, int k, double h);
double eigenvalue(int n, double h, int N); // N to reduce error

double dot_product_with_eigenvector(function_t func, int n, int N);

int main(int argc, char **argv) {
	int N;
	int result;
	double h, error;
	double *solution;
	double *solution_coords;

	if(argc != 2) {
		return 1;
	}

	result = sscanf(argv[1], "%d", &N);
	if(!result)
		return 1;
	h = 2.0 / (2 * N - 1);

	solution = (double*)malloc((unsigned long)(N - 1) * sizeof(double));
	if(!solution)
		return 1;

	solution_coords = (double*)malloc((unsigned long)(N - 1) * sizeof(double));
	if(!solution_coords) {
		free(solution);
		return 1;
	}

	for(int i = 0; i < N - 1; i++)
		solution_coords[i] = dot_product_with_eigenvector(f, i + 1, N)
			/ eigenvalue(i + 1, h, N);

	for(int i = 0; i < N - 1; i++) {
		solution[i] = 0.0;
		for(int j = 0; j < N - 1; j++)
			solution[i] += solution_coords[j] * eigenvector(j + 1, i, h);
	}

	error = 0.0;
	for(int i = 0; i < N - 1; i++)
		error += SQUARE(solution[i] - reference((i + 1) * h));
	error *= h;
	error = sqrt(error);

	printf("%e\n", error);

	free(solution);
	free(solution_coords);
	return 0;
}
