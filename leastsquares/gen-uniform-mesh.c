#include <stdio.h>

int main(int argc, char **argv)
{
    size_t n;
    double a, b;
    if (argc != 4)
        return 1;

    if (sscanf(argv[1], "%lf", &a) != 1)
        return 2;

    if (sscanf(argv[2], "%lf", &b) != 1)
        return 2;

    if (sscanf(argv[3], "%zu", &n) != 1)
        return 2;

    for (size_t i = 0; i < n; i++) {
        printf("%lf ", a + (b - a) * (double)i / (double)(n - 1));
    }
    return 0;
}
