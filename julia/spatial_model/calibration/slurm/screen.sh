#!/bin/bash
# Sobol screening of the ϑ box (spatial_screen.jl) with SLURM_CPUS_PER_TASK − 1 workers.
# Submit from the repository root, e.g.
#   sbatch julia/spatial_model/calibration/slurm/screen.sh --run-dir julia/spatial_model/calibration/runs/screen --first 513 --n 1024
# A requeued job resumes: indices already in <run-dir>/screen.jsonl are skipped.
#SBATCH -J sp-screen
#SBATCH -p econ-grad
#SBATCH -t 4:00:00
#SBATCH -c 64
#SBATCH --mem-per-cpu 2G
#SBATCH --export=ALL
#SBATCH -o julia/spatial_model/calibration/slurm/logs/%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
exec julia --project=julia/spatial_model/calibration/estimation julia/spatial_model/calibration/spatial_screen.jl "$@"
