#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

typedef double (*function_t)(double, double*, int);

#define COORD(i, j, order) ((i) * (order) + (j))

#define EPS 1e-16

double dist(double*, double*, int N);

void adams_moulton(function_t *func, double *solution, double *new_val,
   	int system_order, int N, int k);

void adams_moulton_step(function_t *func, double *solution, double *tmp,
	int system_order, int N, int k);

int main(void) {
	int system_order, N;
	double *solution;
//	double *reference
	return 0;
}

void adams_moulton(function_t *func, double *solution, double *tmp,
	int system_order, int N) {
	// TODO: k <= 2: fill with reference solution
	double eps = 10.0;
	for(int k = 3; k <= N; k++) {
		while(eps >= EPS) {
			adams_moulton_step(func, solution, tmp, system_order, N, k);
			eps = dist(solution + k, tmp, system_order);
			memcpy(solution + k, tmp, system_order * sizeof(double));
		}
	}
}

void adams_moulton_step(function_t *func, double *solution, double *new_val,
		int system_order, int N, int k) {
	assert((k > 2) && (k <= N));

	for (int i = 0; i < system_order, i++) {
		new_val[i] = func[i](((double)k) / ((double)N), solution + k,
			system_order);
	}
}
