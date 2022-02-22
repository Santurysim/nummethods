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

int solve_system(double *matrix, double *result, size_t order) {
    double s, norm1, norm2_square, tmp;

    // Cast the matrix to upper triangular type
    for(size_t i = 0; i < order; i++) {
        s = 0.0;
        for(size_t j = i + 1; j < order; j++){
            s += SQUARE(matrix[COORD(i, j, order)]);
        }

        norm1 = sqrt(SQUARE(matrix[COORD(i, i, order)]) + s);

        if(norm1 < EPS) {
            return 1; // non-invertible matrix
        }

        if(s < EPS) {
            continue; // nothing to do there
        }

        matrix[COORD(i, i, order)] -= norm1;
        norm2_square = SQUARE(matrix[COORD(i, i, order)]) + s;

        norm2_square = 2.0 / norm2_square;

        // Vector of reflection is ready, now we need to operate on matrices

        for(size_t j = i + 1; j < order; j++) {
            s = 0.0;
            for(size_t k = i; k < order; k++) {
                s += matrix[COORD(i, k, order)] * matrix[COORD(j, k, order)];
            }

            s *= norm2_square;
            for(size_t k = i; k < order; k++) {
                matrix[COORD(j, k, order)] -= s * matrix[COORD(i, k, order)];
            }
        }

        s = 0.0;
        for(size_t k = i; k < order; k++) {
            s += matrix[COORD(i, k, order)] * result[k];
        }

        s *= norm2_square;
        for(size_t k = i; k < order; k++) {
            result[k] -= s * matrix[COORD(i, k, order)];
        }

        // Finalize: set the i-th subcolumn of matrix
        matrix[COORD(i, i, order)] = norm1;
    }

    // Back substitution of Gaussian method
    // We know that the matrix is inversible at the moment
    // Note: no action is required on matrix

    for(size_t i = 0/*order - 1*/; i < order/*>= 0*/; i++) {
        // Divide i-th row of result by matrix[i, i]

        s = matrix[COORD(order - 1 - i, order - 1 - i, order)];
        result[i] /= s;

        // Substract (order - 1 - i)-th element of result multiplied by
        // matrix[k, order - 1 - i] from k-th element of result
        // for k = 0, ..., order - 2 - i
        tmp = result[i];
        for(size_t k = 0; k < order - 1 - i; k++) {
            result[k] -= tmp * matrix[COORD(order - 1 - i, k, order)];
        }
    }

    // And... here we go
    return 0;
}
