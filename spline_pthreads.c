#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

typedef struct {
    int start, end;
    double *x, *y, *a, *b, *c, *d;
} thread_data_t;
double get_time_sec() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}
void* compute_spline_thread(void* arg) {
    thread_data_t *data = (thread_data_t*)arg;
    for(int i = data->start; i < data->end; i++){
        data->a[i] = data->y[i];
        data->b[i] = (data->y[i+1] - data->y[i]);
        data->c[i] = 0.0;
        data->d[i] = 0.0;
    }
    return NULL;
}

int main(int argc, char **argv) {
    if(argc < 4) { printf("Usage: %s N NE threads\n", argv[0]); return 1; }

    int N = atoi(argv[1]);
    int NE = atoi(argv[2]);
    int num_threads = atoi(argv[3]);

    double *x = malloc(N*sizeof(double));
    double *y = malloc(N*sizeof(double));
    double *a = malloc((N-1)*sizeof(double));
    double *b = malloc((N-1)*sizeof(double));
    double *c = malloc((N-1)*sizeof(double));
    double *d = malloc((N-1)*sizeof(double));

    for(int i=0;i<N;i++){ x[i] = i*0.01; y[i] = sin(x[i]); }

    pthread_t threads[num_threads];
    thread_data_t tdata[num_threads];
    int chunk = (N-1)/num_threads;

    double t0 = get_time_sec(); 
    for(int t=0;t<num_threads;t++){
        tdata[t].start = t*chunk;
        tdata[t].end = (t==num_threads-1) ? (N-1) : (t+1)*chunk;
        tdata[t].x = x; tdata[t].y = y;
        tdata[t].a = a; tdata[t].b = b; tdata[t].c = c; tdata[t].d = d;
        pthread_create(&threads[t], NULL, compute_spline_thread, &tdata[t]);
    }

    for(int t=0;t<num_threads;t++) pthread_join(threads[t], NULL);
    double t1 = get_time_sec();  
    double sum = 0.0;
    for(int i=0;i<N-1;i++) sum += a[i]+b[i]+c[i]+d[i];
    printf("pthreads checksum=%.10f\n", sum);
    printf("pthreads time = %.6f sec\n", t1 - t0); 
    free(x); free(y); free(a); free(b); free(c); free(d);
    return 0;
}

