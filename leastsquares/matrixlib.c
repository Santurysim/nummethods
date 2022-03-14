/*
 *  Copyright 2020 Peter Shkenev
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <math.h>
#include <string.h>

#include "matrixlib.h"
#include "common.h"

int qr_decompose(double *matrix, double *y, size_t *map, size_t n, size_t m)
{
    double s, norm1, norm2_square, tmp;
    size_t j1;

	for(size_t i = 0; i < m; i++)
		map[i] = i;

    for (size_t i = 0; i < MIN(n, m); i++) {
        j1 = i;
        norm1 = 0.0;
        for (size_t j = i; j < m; j++) {
            s = 0.0;
            for (size_t k = i; k < n; k++)
                s += SQUARE(matrix[COORD(j, k, n)]);
            if (norm1 > s) {
                j1 = j;
                norm1 = s;
            }
        }
        map[j1] = i;
        map[i] = j1;

        // Swap
        for (j = 0; j < n; j++) {
            tmp = matrix[COORD(i, j, n)];
            matrix[COORD(i, j, n)] = matrix[COORD(j1, j, n)];
            matrix[COORD(j1, j, n)] = tmp;
        }

        // norm1 contains squared norm of column i
        s = norm1 - SQUARE(matrix[COORD(i, i, n)]);
        norm1 = sqrt(norm1);

        if (norm1 < EPS) {
            return -1; // non-invertible matrix TODO
        }

        if (s < EPS) {
            continue; // nothing to do there
        }

        matrix[COORD(i, i, n)] -= norm1;
        norm2_square = SQUARE(matrix[COORD(i, i, n)]) + s;

        norm2_square = 2.0 / norm2_square;

        // Vector of reflection is ready, now we need to operate on matrices

        for (size_t j = i + 1; j < m; j++) {
            s = 0.0;
            for (size_t k = i; k < n; k++) {
                s += matrix[COORD(i, k, n)] * matrix[COORD(j, k, n)];
            }

            s *= norm2_square;
            for (size_t k = i; k < n; k++) {
                matrix[COORD(j, k, n)] -= s * matrix[COORD(i, k, n)];
            }
        }

        s = 0.0;
        for (size_t k = i; k < n; k++) {
            s += matrix[COORD(i, k, n)] * y[k];
        }

        s *= norm2_square;
        for (size_t k = i; k < n; k++) {
            y[k] -= s * matrix[COORD(i, k, n)];
        }

        // Finalize: set the i-th subcolumn of matrix
        matrix[COORD(i, i, n)] = norm1;
    }

    return 0;
}

void gauss_back_substitute(double *matrix, double *x, double *y, size_t n,
						   size_t m)
{
    double s, tmp;
    for (size_t i = 0; i < n; i++) {
        // Divide i-th row of result by matrix[i, i]

        s = matrix[COORD(n - 1 - i, n - 1 - i, m)];
        if (fabs(s) < EPS) continue;
        tmp = (result[n - 1 - i] /= s);

        // Substract (order - 1 - i)-th element of result multiplied by
        // matrix[k, order - 1 - i] from k-th element of result
        // for k = 0, ..., order - 2 - i
        for (size_t k = 0; k < order - 1 - i; k++) {
            result[k] -= tmp * matrix[COORD(order - 1 - i, k, order)];
        }
    }
}

int solve_system(double *matrix, double *result, size_t *map, size_t n,
                 size_t m)
{

    // Back substitution of Gaussian method
    // We know that the matrix is inversible at the moment
    // Note: no action is required on matrix


    // And... here we go
    return 0;
}
