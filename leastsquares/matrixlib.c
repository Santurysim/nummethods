/*
 *  Copyright 2020-2022 Peter Shkenev
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
#include <sys/types.h>
#include <assert.h>

#include "matrixlib.h"
#include "common.h"

ssize_t qr_decompose(double *matrix, double *y, size_t *map,
                     size_t n, size_t m)
{
    double s, norm1, norm2_square, tmp;
    size_t j1, tmpi;

	for(size_t i = 0; i < m; i++)
		map[i] = i;

    for (size_t i = 0; i < MIN(n, m); i++) {
        j1 = i;
        norm1 = 0.0;
        for (size_t j = i; j < m; j++) {
            s = 0.0;
            for (size_t k = i; k < n; k++)
                s += SQUARE(matrix[COORD(j, k, n)]);
            if (s > norm1) {
                j1 = j;
                norm1 = s;
            }
        }
        tmpi = map[j1];
        map[j1] = map[i];
        map[i] = tmpi;

        // Swap
        for (size_t j = 0; j < n; j++) {
            tmp = matrix[COORD(i, j, n)];
            matrix[COORD(i, j, n)] = matrix[COORD(j1, j, n)];
            matrix[COORD(j1, j, n)] = tmp;
        }

        // norm1 contains squared norm of column i
        if (norm1 < EPS) {
            return (ssize_t)i; // non-invertible matrix
        }

        s = norm1 - SQUARE(matrix[COORD(i, i, n)]);
        norm1 = sqrt(norm1);

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

        assert(matrix[COORD(i, i, n)] >= 0.0); 
    }

    return MIN(n, m);
}

void gauss_back_substitute(double *matrix, double *x, double *y, size_t n,
						   size_t m)
{
    size_t rank;
    double s, tmp;
    for (rank = 0;
         (rank < MIN(n, m)) && (matrix[COORD(rank, rank, n)] > EPS);
         rank++);
    for(size_t i = 0; i < rank; i++) {

        s = matrix[COORD(rank - 1 - i, rank - 1 - i, n)];
        y[rank - 1 - i] /= s;

        tmp = y[rank - 1 - i];
        for(size_t k = 0; k < rank - 1 - i; k++) {
            y[k] -= tmp * matrix[COORD(rank - 1 - i, k, n)];
        }
        x[rank - 1 - i] = y[rank - 1 - i];
    }
}

ssize_t solve_system(double *matrix, double *x, double *y, size_t *map,
                     size_t n, size_t m)
{
    ssize_t temp;
    memset(x, 0, m * sizeof(double));
    temp = qr_decompose(matrix, y, map, n, m);

    gauss_back_substitute(matrix, x, y, n, m);

    // Remap variables
    for(size_t i = 0; i < m; i++) {
        double tmp;
        tmp = x[i];
        x[i] = x[map[i]];
        x[map[i]] = tmp;

        map[map[i]] = map[i];
        map[i] = i;
    } 

    return 0;
}
