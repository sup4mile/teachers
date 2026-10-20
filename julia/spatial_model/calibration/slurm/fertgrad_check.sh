#!/bin/bash
# Ability-graded fertility shock code checks (runs/fertgrad-20261018/tests/calib_check.jl), one case per job:
# bundled | split, and the smoke slope. Submit from the repository root:
#   sbatch julia/spatial_model/calibration/slurm/fertgrad_check.sh bundled 0.3
#SBATCH -J sp-fertgrad-check
#SBATCH -p econ-grad
#SBATCH -t 04:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH --export=ALL
#SBATCH -o /home/j/jpryan7/slurm_logs/%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
exec julia --project=julia/spatial_model/calibration/estimation \
    julia/spatial_model/calibration/runs/fertgrad-20261018/tests/calib_check.jl "$@"
