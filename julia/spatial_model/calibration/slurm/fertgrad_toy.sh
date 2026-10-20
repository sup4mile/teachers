#!/bin/bash
# Toy-economy checks of the ability-graded fertility shock (runs/fertgrad-20261018/tests/toy_check.jl), serial.
# Submit from the repository root:  sbatch julia/spatial_model/calibration/slurm/fertgrad_toy.sh
#SBATCH -J sp-fertgrad-toy
#SBATCH -p econ-grad
#SBATCH -t 01:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH --export=ALL
#SBATCH -o /home/j/jpryan7/slurm_logs/%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
exec julia --project=julia/spatial_model/calibration/estimation \
    julia/spatial_model/calibration/runs/fertgrad-20261018/tests/toy_check.jl "$@"
