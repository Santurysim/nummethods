#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "geom.h"

int point_cmp(void const *a, void const *b);

int main(int argc, char **argv)
{
    size_t n;
    point_t *points;
    double t1, t2;
    FILE *in;

    if (argc != 2)
        return 1;

    in = fopen(argv[1], "r");
    if (!in)
        return 1;

    n = 0;
    while (fscanf(in, "%d %d", &t1, &t2) == 2)
        n++;
    fclose(in);

    if (n < 3)
        return 1;

    points = (point_t*)malloc(n * sizeof(point_t));
    if (!points)
        return 1;

    lists = (adj_list_t*)malloc(n * sizeof(adj_list_t));
    if (!lists)
        return 1;

    in = fopen(argv[1], "r");
    if (!in)
        return 1;

    for (size_t i = 0; i < n; i++)
        if (fscanf("%d %d", &points[i].x, &points[i].y) != 2)
            return 1;

    fclose(in);

    qsort(points, n, sizeof(point_t), point_cmp);

    return 0;
}

int point_cmp(void const *a, void const *b)
{
    point_t *pa = (point_t*)a;
    point_t *pb = (point_t*)b;

    if (fabs(pa.x - pb.x) < EPS) {
        if (fabs(pa.y - pb.y) < EPS)
            return 0;
        else
            return pa.y < pb.y ? -1 : 1;
    }
    return pa.x < pb.x ? -1 : 1;
}
