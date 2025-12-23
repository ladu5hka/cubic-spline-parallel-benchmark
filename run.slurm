#!/bin/bash
#SBATCH --job-name=spline_all
#SBATCH --output=report.out
#SBATCH --nodes=4
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --partition=cascade

module purge
module load mpi/openmpi/4.1.6/gcc/11
PYTHON_BIN=/opt/software/python/3.11.8/bin/python3

N=50000000
NE=50000000
THREADS=(2 4 8 16 48)
PROCS=(2 4 16 28 112)

echo "===== Job Info ====="
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo "Nodes Allocated   = $SLURM_JOB_NUM_NODES"
echo "Tasks Allocated   = $SLURM_NTASKS"
echo "CPUs per Task     = $SLURM_CPUS_PER_TASK"
echo "==================="

# ----------------------------
# Compile all programs
# ----------------------------
echo "Compiling..."
gcc -O2 -pthread -o spline_pthreads spline_pthreads.c -lm
gcc -O2 -fopenmp -o spline_openmp spline_openmp.c -lm
mpicc -O2 -o spline_mpi spline_mpi.c -lm

echo
echo "=== RUN pthreads ==="
for t in "${THREADS[@]}"; do
  echo "Threads = $t"
  ./spline_pthreads $N $NE $t
done

# ----------------------------
# Run OpenMP (single node)
# ----------------------------
echo
echo "=== RUN OpenMP ==="
for t in "${THREADS[@]}"; do
  echo "Threads = $t"
  OMP_NUM_THREADS=$t ./spline_openmp $N $NE $t
done

# ----------------------------
# Run MPI C (multi-node)
# ----------------------------
echo
echo "=== RUN MPI C ==="
for p in "${PROCS[@]}"; do
  echo "Procs = $p"
  mpirun --mca btl tcp,self -np $p ./spline_mpi $N $NE
done

# ----------------------------
# Run MPI Python (multi-node)
# ----------------------------
echo
echo "=== RUN MPI Python ==="
for p in "${PROCS[@]}"; do
  echo "Procs = $p"
  mpirun --mca btl tcp,self -np $p $PYTHON_BIN spline_mpi.py $N $NE
done
