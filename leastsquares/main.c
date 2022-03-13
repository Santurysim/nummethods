#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <matrixlib.h>

typedef double (*function_t)(double);

// TODO give a system of functions



int main(int argc, char **argv)
{
    int result;
    size_t n, N;
    double *x, *y, *alpha, *matrix, *solution, t1, t2, t3;
    FILE *in;

    if(argc != 2) {
        return 1;
    }

    in = fopen(argv[1], "r");
    if(!in) {
        perror("fatal: fopen");
        return 1;
    }

    while((result = fscanf(in, "%lf %lf %lf", &t1, &t2, &t3) != EOF)) {
        if(result != 3) {
            fprintf(stderr, "fatal: malformed data\n");
            fclose(in);
            return 1;
        }
        N++;
    }

    fclose(in);

    x = (double*)malloc(N * sizeof(double));
    if(!x) {
        perror("fatal: malloc");
        return 2;
    }

    y = (double*)malloc(N * sizeof(double));
    if(!y) {
        perror("fatal: malloc");
        free(x);
        return 2;
    }

    alpha = (double*)malloc(N * sizeof(double));
    if(!alpha) {
        perror("fatal: malloc");
        free(x);
        free(y);
        return 2;
    }

    in = fopen(argv[1], "r");
    if(!in) {
        perror("fatal: fopen");
        free(x);
        free(y);
        free(alpha);
        return 1;
    }

    for(size_t i = 0; i < N; i++) {
        result = fscanf("%lf %lf %lf", x + i, y + i, alpha + i);
        if(result != 3) {
            fprintf(stderr, "fatal: malformed data\n");
            free(x);
            free(y);
            free(alpha);
            fclose(in);
            return 1;
        }
    }

    fclose(in);

    // 1. Reduce prbolem to case of alpha = (1.0, ..., 1.0)
    // This can be achieved by multiplying y[i] by sqrt(alpha[i])
    // and dividing obtained solution[i] by sqrt(alpha[i])

    for(size_t i = 0; i < N; i++)
        y[i] *= sqrt(alpha[i]);

    solve_system(...);
    
    free(x);
    free(y);
    free(alpha);
    return 0;
}
