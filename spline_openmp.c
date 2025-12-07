#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

double get_time_sec() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv){
    if(argc<4){ printf("Usage: %s N NE threads\n",argv[0]); return 1; }

    int N = atoi(argv[1]);
    int NE = atoi(argv[2]);
    int threads = atoi(argv[3]);

    double *x = malloc(N*sizeof(double));
    double *y = malloc(N*sizeof(double));
    double *a = malloc((N-1)*sizeof(double));
    double *b = malloc((N-1)*sizeof(double));
    double *c = malloc((N-1)*sizeof(double));
    double *d = malloc((N-1)*sizeof(double));

    for(int i=0;i<N;i++){ x[i]=i*0.01; y[i]=sin(x[i]); }

    omp_set_num_threads(threads);

    double t0 = get_time_sec();  
#pragma omp parallel for
    for(int i=0;i<N-1;i++){
        a[i] = y[i];
        b[i] = y[i+1]-y[i];
        c[i] = 0.0;
        d[i] = 0.0;
    }

    double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
    for(int i=0;i<N-1;i++) sum += a[i]+b[i]+c[i]+d[i];
    double t1 = get_time_sec();
    printf("openmp checksum=%.10f\n",sum);
    printf("openmp time=%.6f seconds\n", t1 - t0);
    free(x); free(y); free(a); free(b); free(c); free(d);
    return 0;
}

