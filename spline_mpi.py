#!/opt/software/python/3.11.8/bin/python3
from mpi4py import MPI
import numpy as np
import math
import sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if len(sys.argv) != 3:
    if rank == 0:
        print("Usage: python spline_mpi.py N NE")
    sys.exit(1)

N = int(sys.argv[1])
NE = int(sys.argv[2])  

x = np.arange(N, dtype=np.float64) * 0.01
y = np.sin(x)

a = np.copy(y)
b = y[1:] - y[:-1]
c = np.zeros(N-1, dtype=np.float64)
d = np.zeros(N-1, dtype=np.float64)

chunk = (N - 1) // size
start = rank * chunk
end = (rank + 1) * chunk if rank != size - 1 else (N - 1)

t0 = MPI.Wtime()

local_sum = np.sum(a[start:end] + b[start:end] + c[start:end] + d[start:end])

t1 = MPI.Wtime()

checksum = comm.reduce(local_sum, op=MPI.SUM, root=0)

if rank == 0:
    print(f"mpi Python checksum={checksum:.10f}")
    print(f"mpi Python time = {t1 - t0:.6f} sec")

