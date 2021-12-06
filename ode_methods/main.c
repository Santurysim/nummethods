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

void adams_moulton_step(function_t *func, double *solution, double *tmp,
	int order, int N, int k);

void shift(double *array, int n, int m);

int main(int argc, char **argv) {
	int order, N, result;
	double *approximation_last, *approximation2_last, *tmp;
	double h;
	double *reference;
	function_t FUNCTIONS[3] = {f1, f2, f3};
	order = 3;

	if(argc != 2) {
		return 1;
	}

	result = sscanf(argv[1], "%d", &N);
	if(result != 1) {
		return 1;
	}
	h = 1.0 / N;

	approximation_last = (double*)malloc(4ll * order * sizeof(double));
	if(!approximation_last) {
		return 1;
	}

	approximation2_last = (double*)malloc(4ll * order * sizeof(double));
	if(!approximation2_last) {
		free(approximation_last);
		return 1;
	}

	// TODO Initial conditions in 3 points
	for(int i = 0; i < 3; i++) {
		approximation_last[i * order] = ref1((double)i * h);
		approximation_last[i * order + 1] = ref2((double)i * h);
		approximation_last[i * order + 2] = ref3((double)i * h);

		approximation2_last[
	}
	
	free(approximation_last);
	free(approximation2_last);
	return 0;
}

void adams_moulton_step(function_t *func, double *solution, double *new_val,
		int order, int N, int k) {
	double difference;
	assert((k > 2) && (k <= N));
	// TODO Initial values

	for(;;) {
		// Obtain new vector for y_k
		for (int i = 0; i < order; i++) {
			new_val[i] = 9.0 * func[i]((double)k / N, solution + k * order,
					order)
				+ 19.0 * func[i]((double)(k - 1) / N,
					solution + (k - 1) * order,order)
				- 5.0 * func[i]((double)(k - 2) / N, solution + k - 2, order)
				+ func[i]((double)(k - 3) / N, solution + k - 3, order);
			new_val[i] /= (24.0 * N);
			new_val[i] += solution[COORD(k, i, order)];
		}
		difference = dist(solution + k * order, new_val, order);
		memcpy(solution + k * order, new_val, order * sizeof(double));
		if(difference < EPS) {
			break;
		}
	}
}

void shift(double *array, int n, int m) {
	for(int i = 0; i < n - 1; i++) {
		memcpy(array + i * m, array + (i + 1) * m, m * sizeof(double)); 
	}
}

double dist(double* x, double *y, int n) { // Euclidean norm for now
	double result_square = 0.0;
	for(int i = 0; i < n; i++) {
		result_square += (x[i] - y[i]) * (x[i] - y[i]);
	}
	return sqrt(result_square);
}
