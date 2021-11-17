#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "schemes.h"

double scheme_error(double *solution, double *approximation, int N);

int main(void) {
	return 0;
}

double scheme_error(double *solution, double *approximation, int N) {
	double t;
	double result = 0.0;
	for (int i = 0; i <= N; i++) {
		t = fabs(solution[i] - approximation[i]);
		if (t > result) {
			result = t;
		}
	}
	return result;
}

