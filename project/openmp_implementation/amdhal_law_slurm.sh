#!/bin/bash
#SBATCH --job-name=fft_scaling               # Job name
#SBATCH --output=fft_scaling_results_%j.out  # Output file for Slurm job stdout
#SBATCH --error=fft_scaling_results_%j.err   # Error file for Slurm job stderr
#SBATCH --ntasks=1                           # Single task (using OpenMP only)
#SBATCH --cpus-per-task=16                   # Allocate 16 CPU cores
#SBATCH --time=01:00:00                      # Maximum job runtime

# Define the output results file
results_file="fft_scaling_results_$SLURM_JOB_ID.txt"

# Write header to results file
echo "Threads   Code Output" > "$results_file"

threads_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)

for threads in "${threads_list[@]}"; do
    if [ $threads -le $SLURM_CPUS_PER_TASK ]; then
        export OMP_NUM_THREADS=$threads
        echo "-----------------------------------------" >> "$results_file"
        echo "Running with $threads thread(s)..." >> "$results_file"
        
        # Run the code and append its output (including timing) to the results file
        ./fft 23 3 >> "$results_file" 2>&1

        echo "Finished run with $threads thread(s)." >> "$results_file"
    else
        echo "Skipping $threads threads as it exceeds allocated CPUs ($SLURM_CPUS_PER_TASK)" >> "$results_file"
    fi
done

echo "Scaling analysis completed. Results saved to $results_file." >> "$results_file"

