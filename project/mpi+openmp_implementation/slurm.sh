#!/bin/bash
#SBATCH --job-name=fft_scaling               # Job name
#SBATCH --output=fft_scaling_results_%j.out  # Output file for job results, %j is replaced with job ID
#SBATCH --error=fft_scaling_results_%j.err   # Error file for job errors
#SBATCH --ntasks=1                           # Single task (using OpenMP, not MPI)
#SBATCH --cpus-per-task=16

#SBATCH --time=01:00:00                      # Maximum job runtime

# Define the output results file
results_file="fft_scaling_results_$SLURM_JOB_ID.txt"

echo "Threads Runtime" > $results_file

threads_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)

for threads in "${threads_list[@]}"; do
    if [ $threads -le $SLURM_CPUS_PER_TASK ]; then
        export OMP_NUM_THREADS=$threads
        echo "Running with $threads thread(s)..."
        
        # Measure the runtime using `time` command or internal timing in code
        start_time=$(date +%s.%N)
        ./fft 22 4
        end_time=$(date +%s.%N)
        
      
        runtime=$(echo "$end_time - $start_time" | bc -l)
        
        # Save the results to the file
        echo "$threads $runtime" >> $results_file
    else
        echo "Skipping $threads threads as it exceeds allocated CPUs ($SLURM_CPUS_PER_TASK)"
    fi
done

echo "Scaling analysis completed. Results saved to $results_file."

