#!/bin/bash
#SBATCH --job-name=gustafson_scaling          # Job name
#SBATCH --output=gustafson_scaling_%j.out     # Slurm job stdout
#SBATCH --error=gustafson_scaling_%j.err      # Slurm job stderr
#SBATCH --nodes=2 
#SBATCH --ntasks-per-node=1                   # Single MPI task (using OpenMP)
#SBATCH --cpus-per-task=16                    # Up to 32 CPUs for max scale
#SBATCH --time=01:00:00                       # Max runtime

# The output results file
results_file="gustafson_scaling_results.txt"

# Write a header to the results file
echo "N (2^k)   Threads   Runtime Output" > "$results_file"

# Define the problem sizes and corresponding threads
# k values: 19 to 24
k_values=(19 20 21 22 23 24)
thread_values=(1 2 4 8 16 32)

# Loop over the indices of these arrays
for i in "${!k_values[@]}"; do
    k=${k_values[$i]}
    threads=${thread_values[$i]}
    
    # Check that we have enough CPUs allocated
    if [ $threads -le $SLURM_CPUS_PER_TASK ]; then
        export OMP_NUM_THREADS=$threads
        N=$((1<<$k))
        
        echo "-----------------------------------------" >> "$results_file"
        echo "Running for N=2^$k (N=$N) with $threads thread(s)..." >> "$results_file"
        
        # Run the code in mode 3 for timing results
        ./fft $k 3 >> "$results_file" 2>&1
        
        echo "Finished run for N=2^$k with $threads thread(s)." >> "$results_file"
    else
        echo "Skipping N=2^$k with $threads threads as it exceeds allocated CPUs ($SLURM_CPUS_PER_TASK)" >> "$results_file"
    fi
done

echo "Gustafson's scaling analysis completed. Results saved to $results_file." >> "$results_file"

