#!/bin/bash
#SBATCH --job-name=test_helper_function         # Job name
#SBATCH --output=output_%j.log      # Standard output log (%j = job ID)
#SBATCH --error=error_%j.log        # Error log
#SBATCH --ntasks=1                  # Number of tasks (1 for serial jobs)
#SBATCH --cpus-per-task=1           # Number of CPU cores
#SBATCH --mem=96G                    # Memory per node
#SBATCH --time=00:10:00             # Time limit (hh:mm:ss)
#SBATCH --partition=debug        # Partition name (adjust as needed)

# Load R module if needed
module load R

# Run the R script
Rscript R/0_helper_fct.R

