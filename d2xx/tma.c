#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SQUARE(x) ((x)*(x))

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

	solution = (double*)malloc((long unsigned)(N - 2) * sizeof(double*));
	if(!solution) {
		free(alpha);
		free(beta);
		return 1;
	}

	alpha[0] = 1.0 / (2.0 + b(h) * h * h);
	beta[0] = (f(h) * h * h) / (2.0 + b(h) * h * h);

	for(int i = 0; i < N - 3; i++) {
		double denominator = (2.0 + b((i + 2) * h) * h * h - alpha[i]);
		alpha[i + 1] = 1.0 / denominator;
		beta[i + 1] = (f((i + 2) * h) * h * h + beta[i]) / denominator;
	}

	// Obtain solution
	solution[N - 2] = (f((N - 1) * h) * h * h + beta[N - 3])
		/ (3 + b((N - 1) * h) * h * h - alpha[N - 3]);

	for(int i = 0; i < N - 2; i++) {
		solution[N - i - 3] = alpha[N - i - 3] * solution[N - i - 2]
			+ beta[N - i - 3];
	}

	error = 0.0;
	for(int i = 0; i < N - 1; i++)
		error += SQUARE(solution[i] - reference((i + 1) * h));
	error *= h;
	error = sqrt(error);

	printf("%e\n", error);

//	for(int i = 0; i < N - 1; i++) {
//		printf("%e\n", solution[i]);
//	}

	free(alpha);
	free(beta);
	free(solution);
	return 0;
}

double f(double x) {
	return 2.0 * M_PI * M_PI * exp(M_PI * x) * cos(M_PI * x)
		* (sin(M_PI * x) - 1);
}

double b(double x) {
	return 2.0 * M_PI * M_PI * cos(M_PI * x);
}

double reference(double x) {
	return exp(M_PI * x) * sin(M_PI * x);
}
