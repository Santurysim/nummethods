#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

typedef double (*function_t)(double, double*, int);

#define COORD(i, j, order) ((i) * (order) + (j))

#define EPS 1e-16

extern double f1(double x, double *y, int n);
extern double f2(double x, double *y, int n);
extern double f3(double x, double *y, int n);

extern double ref1(double x);
extern double ref2(double x);
extern double ref3(double x);

double dist(double*, double*, int N);

void adams_moulton(function_t *func, double *solution, double *new_val,
   	int order, int N, int k);

void adams_moulton_step(function_t *func, double *solution, double *tmp,
	int order, int N, int k);

int main(int argc, char **argv) {
	int order, N, result;
	double *approximation, *approximation2;
	double *reference;
	order = 3;

	if(argc != 2) {
		return 1;
	}

	result = sscanf(argv[1], "%d", &N);
	if(result != 1) {
		return 1;
	}

	approximation = malloc((long long)N * order * sizeof(double));
	if(!approximation) {
		return 1;
	}

	approximation2 = malloc((long long) 2 * N * order * sizeof(double));
	if(!approximation2) {
		free(approximation1);
		return 1;
	}

	// Do things;
	
	free(approximation);
	free(approximation2);
	return 0;
}

void adams_moulton(function_t *func, double *solution, double *tmp,
	int order, int N) {
	// TODO: k <= 2: fill with reference solution
	double eps = 10.0;
	for(int k = 3; k <= N; k++) {
		while(eps >= EPS) {
			adams_moulton_step(func, solution, tmp, order, N, k);
			diff = dist(solution + k * order, tmp, order);
			memcpy(solution + k * order, tmp, order * sizeof(double));
		}
	}
}

void adams_moulton_step(function_t *func, double *solution, double *new_val,
		int order, int N, int k) {
	assert((k > 2) && (k <= N));

	for (int i = 0; i < order, i++) {
		new_val[i] = 9.0 * func[i]((double)k / N, solution + k * order, order)
			+ 19.0 * func[i]((double)(k - 1) / N, solution + (k - 1) * order,
					order)
			- 5.0 * func[i]((double)(k - 2) / N, solution + k - 2, order)
			+ func[i]((double)(k - 3) / N, solution + k - 3, order);
		new_val[i] /= (24.0 * N);
		new_val[i] += solution[COORD(k, i, order)];
	}
}

double dist(double* x, double *y, int n) { // Euclidean norm for now
	double result_square = 0.0;
	for(int i = 0; i < n; i++) {
		result_square += (x[i] - y[i]) * (x[i] - y[i]);
	}
	return sqrt(result_square);
}
