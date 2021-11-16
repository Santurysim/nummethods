#include <stdio.h>
#include <math.h>
#include <stdio.h>

#define SQUARE(x) ((x) * (x))

double eigenvector(int n, int k, double h);
double eigenvalue(int n, double h, int N); // N to reduce error

int main(int argc, char **argv) {
	int result, N;
	double h, residual1, residual2, temp;

	if (argc < 2) {
		fprintf(stderr, "Specify mesh size\n");
		return 1;
	}

	result = sscanf(argv[1], "%d", &N);
	if (result != 1) {
		fprintf(stderr, "Bad mesh size\n");
		return 1;
	}
	if (N < 3) {
		fprintf(stderr, "Small mesh size\n");
		return 1;
	}

	h = 2 / (double)(2 * N - 1);


	// Compute residual 1
	residual1 = -INFINITY;
	for (int i = 1; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			temp = .0;
			for (int k = 1; k < N; k++)
				temp += eigenvector(i, k, h) * eigenvector(j, k, h);
			if (temp > residual1)
				residual1 = temp;
		}
	}

	// Compute residual 2
	residual2 = .0;
	for (int i = 1; i < N; i++) {
		temp = SQUARE((-2 * eigenvector(i, 1, h) +
				eigenvector(i, 2, h)) / (h * h) +
				eigenvalue(i, h, N) * eigenvector(i, 1, h));
		for (int k = 2; k < N - 1; k++) {
			temp += SQUARE((eigenvector(i, k - 1, h) -
				2 * eigenvector(i, k, h) +
				eigenvector(i, k + 1, h)) / (h * h)  +
				eigenvalue(i, h, N) * eigenvector(i, k, h));
		}
		temp += SQUARE((eigenvector(i, N - 2, h) -
			3 * eigenvector(i, N - 1, h)) / (h * h) +
			eigenvalue(i, h, N) * eigenvector(i, N - 1, h));
		temp = sqrt(temp) / eigenvalue(i, h, N);
		if (temp > residual2)
			residual2 = temp;
	}

	printf("Residual 1: %e\nResidual 2: %e\n", residual1, residual2);
	return 0;
}

double eigenvector(int n, int k, double h) {
	return sin(M_PI * n * k * h) * sqrt(2);
}

double eigenvalue(int n, double h, int N) {
	return SQUARE(sin(M_PI * n * h / 2)) * SQUARE(2 * N - 1);
}
