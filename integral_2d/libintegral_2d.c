#include <math.h>

#include "libintegral_2d.h"

double triangle_integral(function_t func, point_t a, point_t b, point_t c)
{
    double ab_x, ab_y, ac_x, ac_y, s, t;
    ab_x = b.x - a.x;
    ab_y = b.y - a.y;
    ac_x = c.x - a.x;
    ac_y = c.y - a.y;
    s = fabs(ab_x * ac_y - ab_y * ac_x) / 2.0;
    t = func((a.x + b.x) / 2.0, (a.y + b.y) / 2.0)
        + func((a.x + c.x) / 2.0, (a.y + c.y) / 2.0)
        + func((b.x + c.x) / 2.0, (b.y + c.y) / 2.0);
    return (s * t) / 3.0;
}

double integral_2d(function_t func, point_t bottom, point_t top, size_t nx,
                   size_t ny)
{
    point_t a = POINT_ZERO;
    point_t b = POINT_ZERO;
    point_t c = POINT_ZERO;
    double res = 0.0;
    for (size_t i = 0; i < nx, i++) {
        a.x = bottom.x + (top.x - bottom.x) * ((double)i / (double) nx);
        a.y = bottom.y;
        b.x = a.x;
        b.y = bottom.y + (top.y - bottom.y) / (double)ny;
        c.x = bottom.x + (top.x - bottom.x) * (double)(i + 1) / (double) nx;
        c.y = b.y;
        for (size_t j = 0; j < ny; j++) {
            res += triangle_integral(func, a, b, c);
            a.y = b.y;
            b.y += (top.y - bottom.y) / (double)ny;
            c.y = b.y;
        }
    }
}
