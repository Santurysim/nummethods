#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double f(double);

double b(double);

double reference(double);

int main(int argc, char **argv) {
	int N;
	int result;
	double h, error;
	double *solution;
	double *alpha, *beta;

	if(argc != 2) {
		return 1;
	}

	result = sscanf(argv[1], "%d", &N);
	if(!result)
		return 1;
	h = 2.0 / (2 * N - 1);

	alpha = (double*)malloc((long unsigned)(N - 1) * sizeof(double*));
	if(!alpha)
		return 1;

	beta = (double*)malloc((long unsigned)(N - 1) * sizeof(double*));
	if(!beta) {
		free(alpha);
		return 1;
	}

	solution = (double*)malloc((long unsigned)(N - 1) * sizeof(double*));
	if(!solution) {
		free(alpha);
		free(beta);
		return 1;
	}

	alpha[0] = 2.0 - b(h) * h * h;
	beta[0] = (f(h) * h * h) / (b(h) * h * h - 2.0);

	for(int i = 1; i < N - 1; i++) {
		double denominator = (2.0 - b((i + 1) * h) * h * h - alpha[i - 1]);
		alpha[i] = 1.0 / denominator;
		beta[i] = (-f((i + 1) * h) * h * h + beta[i - 1]) / denominator;
	}

	// Obtain solution
	solution[N - 2] = (-f((N - 1) * h) * h * h + beta[N - 2])
		/ (3 - b((N - 1) * h) * h * h - alpha[N - 2]);

	for(int i = 0; i < N - 2; i++) {
		solution[N - 3 - i] = alpha[N - 2 - i] * solution[N - 2 - i]
			+ beta[N - 2 - i];
	}

	error = 0.0;
	for(int i = 0; i < N - 1; i++)
		error += solution[i] * reference((i + 1) * h);
	error *= h;
	error = sqrt(error);

	printf("%e\n", error);

	free(alpha);
	free(beta);
	free(solution);
	return 0;
}
