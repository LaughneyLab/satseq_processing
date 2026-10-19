#!/bin/bash
#SBATCH --job-name=MAST_run
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --error=logs/array_%A_%a.err
#SBATCH --array=1-382                # <--- Defines the range of i
#SBATCH --ntasks=1                   # <--- 1 core per task
#SBATCH --cpus-per-task=4            # <--- Change if MAST uses multi-threading
#SBATCH --mem=64G                     # <--- Adjust memory per task
#SBATCH --time=02:00:00              # <--- Adjust time per task

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Load R module
module load R

# Run the script with the array index as the argument
Rscript --vanilla run_MAST_interaction_variant_cGAMP.R $SLURM_ARRAY_TASK_ID
