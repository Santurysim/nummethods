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

#pragma once

int solve_system(double *matrix, double *x, double *y, size_t *map, size_t n,
                 size_t m);

int qr_decompose(double *matrix, double *y, size_t *map, size_t n, size_t m);

void gauss_back_substitute(double *matrix, double *x, double *y, size_t n,
						   size_t m);
