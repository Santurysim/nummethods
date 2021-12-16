#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>

typedef double (*function_t)(double, double*, int);

#define COORD(i, j, order) ((i) * (order) + (j))

#define EPS 1e-14

extern double f1(double x, double *y, int n);
extern double f2(double x, double *y, int n);
extern double f3(double x, double *y, int n);

extern double ref1(double x);
extern double ref2(double x);
extern double ref3(double x);

double dist(double *x, double *y, int N);

void adams_moulton_step(function_t *func, double *solution, double *tmp,
	int order, int N, int k);

void shift(double *array, int n, int m);

void print_line(double x, double *y, int order, double error,
	double runge_error);

int main(int argc, char **argv) {
	int N, result, order;
	double *approximation_last, *approximation2_last, *tmp;
	double h, runge_coeff, error;
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

	approximation_last = (double*)malloc(4 * (unsigned int)order
		* sizeof(double));
	if(!approximation_last) {
		return 1;
	}

	approximation2_last = (double*)malloc(4 * (unsigned int)order
		* sizeof(double));
	if(!approximation2_last) {
		free(approximation_last);
		return 1;
	}

	tmp = (double*)malloc((unsigned int)order * sizeof(double));
	if(!tmp) {
		free(approximation_last);
		free(approximation2_last);
		return 1;
	}

	for(int i = 0; i < 3; i++) {
		approximation_last[i * order] = ref1(i * h);
		approximation_last[i * order + 1] = ref2(i * h);
		approximation_last[i * order + 2] = ref3(i * h);
	}

	print_line(0.0, approximation_last, order, 0.0, 0.0);
	print_line(h, approximation_last + order, order, 0.0, 0.0);
	print_line(2 * h, approximation_last + 2 * order, order, 0.0, 0.0);
	
	for(int i = 3; i <= N; i++) {
		adams_moulton_step(FUNCTIONS, approximation_last, tmp, order, N, i);

		tmp[0] = ref1(i * h);
		tmp[1] = ref2(i * h);
		tmp[2] = ref3(i * h);
		error = dist(approximation_last + 3 * order, tmp, order);

		//two steps for approximation2
		
		memcpy(approximation2_last, approximation_last,
			3 * (unsigned long)order * sizeof(double));

		adams_moulton_step(FUNCTIONS, approximation2_last, tmp, order, 2 * N,
			2 * i - 1);
		shift(approximation2_last, 4, order);
		adams_moulton_step(FUNCTIONS, approximation2_last, tmp, order, 2 * N,
			2 * i);

		// runge coefficient
		runge_coeff = dist(approximation_last + 3 * order,
			approximation2_last + 3 * order, order) / 15.0;

		shift(approximation_last, 4, order);

		print_line(i * h, approximation_last + 3 * order, order, error,
			runge_coeff * h * h * h * h);
	}

	free(tmp);
	free(approximation_last);
	free(approximation2_last);
	return 0;
}

void print_line(double x, double *y, int order, double error,
	double runge_error) {
	printf("%lf ", x);
	for(int i = 0; i < order; i++) {
		printf("%lf ", y[i]);
	}
	printf("%e %e\n", error, runge_error);
}

void adams_moulton_step(function_t *func, double *solution, double *new_val,
		int order, int N, int k) {
	double difference, h;
	assert((k > 2) && (k <= N));
	h = 1.0 / N;

	// Initial value - value at previous point
	memcpy(solution + 3 * order, solution + 2 * order, (unsigned long)order
		* sizeof(double));

	for(;;) {
		// Obtain new vector for y_k
		for (int i = 0; i < order; i++) {
			new_val[i] = 9.0 * func[i](k * h, solution + 3 * order, order)
				+ 19.0 * func[i]((k - 1) * h, solution + 2 * order, order)
				- 5.0 * func[i]((k - 2) * h, solution + order, order)
				+ func[i]((k - 3) * h, solution, order);
			new_val[i] *= (h / 24.0);
			new_val[i] += solution[2 * order + i];
		}
		difference = dist(solution + 3 * order, new_val, order);
		memcpy(solution + 3 * order, new_val, (unsigned long)order
			* sizeof(double));
		if(difference < EPS) break;
//		printf("%e\n", difference);
	}
}

void shift(double *array, int n, int m) {
	for(int i = 0; i < n - 1; i++) {
		memcpy(array + i * m, array + (i + 1) * m,
			(unsigned long)m * sizeof(double)); 
	}
}

double dist(double* x, double *y, int n) { // Euclidean norm for now
	double result = 0.0;
	for(int i = 0; i < n; i++) {
		result += fabs(x[i] - y[i]);
	}
	return result;
}
