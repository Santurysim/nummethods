#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vector>
#include <algorithm>

#include "geom.h"

using namespace std;

int point_cmp(void const *a, void const *b);

int main(int argc, char **argv)
{
    vector<Point> points;
    double t1, t2;
    FILE *in;

    if (argc != 2)
        return 1;

    in = fopen(argv[1], "r");
    if (!in)
        return 1;

    while (fscanf(in, "%d %d", &t1, &t2) == 2) {
        Point p;
        p.x = t1; p.y = t2;
        points.add(point);
    }
    fclose(in);

    sort(points.begin(), points.end, point_cmp);

    return 0;
}

bool point_cmp(const Point &pa, const Point &pb)
{
    if (fabs(pa.x - pb.x) < EPS) {
        if (fabs(pa.y - pb.y) < EPS)
            return 0;
        else
            return pa.y < pb.y ? -1 : 1;
    }
    return pa.x < pb.x ? -1 : 1;
}
