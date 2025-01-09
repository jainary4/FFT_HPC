#!/bin/bash
#SBATCH --job-name=fft_scaling_mpiexec
#SBATCH --output=fft_timing_%j.txt
#SBATCH --error=fft_scaling_mpiexec_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00

# Load necessary modules (adjust as needed)
module load gcc/12.2.0
module load openmpi/4.1.4

# Set the number of OpenMP threads if needed
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Parameters
EXEC="./fft_hybrid"
K=25
MODE=3
N=$((1<<$K))  # 2^25

echo "Running FFT scaling test for N=$N (2^$K) mode=$MODE using mpirun"
echo "Number of OpenMP threads: $OMP_NUM_THREADS"

# We have allocated 2 nodes with 4 tasks each, total 8 tasks.
# We'll test a few different MPI process counts:
for p in 1 2 4 6 8
do
    echo "-----------------------------------------"
    echo "Running with $p MPI tasks..."
    mpirun -n $p $EXEC $K $MODE
    echo "Finished run with $p MPI tasks."
done

echo "All runs completed."

