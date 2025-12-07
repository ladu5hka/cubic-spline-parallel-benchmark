#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

double get_time_sec() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv){
    MPI_Init(&argc,&argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if(argc<3){
        if(rank==0) printf("Usage: %s N NE\n",argv[0]);
        MPI_Finalize(); return 1;
    }

    int N = atoi(argv[1]);
    int NE = atoi(argv[2]);

    double *x = malloc(N*sizeof(double));
    double *y = malloc(N*sizeof(double));
    double *a = malloc((N-1)*sizeof(double));
    double *b = malloc((N-1)*sizeof(double));
    double *c = malloc((N-1)*sizeof(double));
    double *d = malloc((N-1)*sizeof(double));

    for(int i=0;i<N;i++){ x[i]=i*0.01; y[i]=sin(x[i]); }

    int chunk = (N-1)/size;
    int start = rank*chunk;
    int end = (rank==size-1)?(N-1):(rank+1)*chunk;
    MPI_Barrier(MPI_COMM_WORLD);    
    double t0 = get_time_sec();

    for(int i=start;i<end;i++){
        a[i] = y[i];
        b[i] = y[i+1]-y[i];
        c[i] = 0.0;
        d[i] = 0.0;
    }

    double local_sum = 0.0;
    for(int i=start;i<end;i++) local_sum += a[i]+b[i]+c[i]+d[i];

    double global_sum = 0.0;
    MPI_Reduce(&local_sum,&global_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = get_time_sec();
    if(rank == 0){
        printf("mpi C checksum=%.10f\n", global_sum);
        printf("mpi C time=%.6f seconds\n", t1 - t0);
    }
    free(x); free(y); free(a); free(b); free(c); free(d);
    MPI_Finalize();
    return 0;
}

