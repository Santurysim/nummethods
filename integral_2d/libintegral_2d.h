#ifndef LIBINTEGRAL_2D_H
#define LIBINTEGRAL_2D_H

typedef struct point {
    double x;
    double y;
} point_t;

#define POINT_ZERO { .x = 0.0, .y = 0.0 }

typedef double (*function_t)(double, double);

double triangle_integral(function_t func, point_t *a, point_t *b, point_t *c);

double integral_2d(function_t func, point_t bottom, point_t top, size_t nx,
                   size_t ny);

#endif /* LIBINTEGRAL_2D_H */
