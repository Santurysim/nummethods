#define HAZARDOUS_UPPER_BOUND 1e307
#define HAZARDOUS_LOWER_BOUND 1e-303

// log(1e-307) rounded up
#define HAZARDOUS_LOWER_BOUND_LOG -706.0

void compute_scheme1(double *arr, int N, double A);
void compute_scheme2(double *arr, int N, double A);
void compute_scheme3(double *arr, int N, double A);
void compute_scheme4(double *arr, int N, double A);
void compute_scheme5(double *arr, int N, double A);
void compute_scheme6(double *arr, int N, double A);

void compute_reference(double *arr, int N, double A);
