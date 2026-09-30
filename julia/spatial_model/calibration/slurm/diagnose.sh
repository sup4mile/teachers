#!/bin/bash
# Diagnostics at a point ϑ (spatial_diagnose.jl) with SLURM_CPUS_PER_TASK − 1 workers:
# 17 CPUs cover a forward-difference Jacobian of the 8 parameters twice over, or a
# central-difference Jacobian in one round. Submit from the repository root, e.g.
#   sbatch julia/spatial_model/calibration/slurm/diagnose.sh --theta theta0 --consistency --jacobian
#   sbatch julia/spatial_model/calibration/slurm/diagnose.sh --theta julia/spatial_model/calibration/runs/<id> \
#          --polish --jacobian --consistency --out julia/spatial_model/calibration/runs/<id>/diagnose
#SBATCH -J sp-diagnose
#SBATCH -p econ-grad
#SBATCH -t 12:00:00
#SBATCH -c 17
#SBATCH --mem-per-cpu 2G
#SBATCH --export=ALL
#SBATCH -o julia/spatial_model/calibration/slurm/logs/%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
exec julia --project=julia/spatial_model/calibration/estimation julia/spatial_model/calibration/spatial_diagnose.jl "$@"
