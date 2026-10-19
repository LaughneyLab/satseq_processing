#!/bin/bash
# ============================================================
#  submit_MAST_celltype.sh
#
#  Submits a SLURM array job to run MAST for every combination
#  of cell type (11) x variant (3) = 33 tasks, all vs REF.
#
#  Each task runs:
#    Rscript --vanilla run_MAST_celltype_variant_vs_REF.R <task_id>
#
#  After all tasks complete, run the collation step:
#    Rscript --vanilla collate_MAST_to_excel.R
# ============================================================

#SBATCH --job-name=MAST_ct_var
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --error=logs/array_%A_%a.err
#SBATCH --array=1-36                 # 3 variants x 12 cell types = 36 tasks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00              # MAST on a cell-type subset is usually <2h;
                                     # 4h gives comfortable headroom

# Prevent BLAS/LAPACK from spawning extra threads
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Create log and output directories if they don't exist
mkdir -p logs
mkdir -p MAST_output/raw

# Load R (adjust module name to match your HPC environment)
module load R/4.5.2

# Run the analysis for this array task
Rscript --vanilla run_MAST_celltype_variant_vs_REF.R $SLURM_ARRAY_TASK_ID
