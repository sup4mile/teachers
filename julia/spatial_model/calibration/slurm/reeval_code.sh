#!/bin/bash
# One evaluation at a point ϑ (spatial_diagnose.jl without --polish: report.txt, theta.toml and moments.toml in --out)
# on a frozen code tree runs/_code/<name> (CODE), or on the main tree when CODE is unset; serial, one worker-free
# process. Used for bit-for-bit inertness checks (2026-10-26, constant class size). Submit from the repository root:
#   CODE=julia/spatial_model/calibration/runs/_code/pre-ccs-20261026 sbatch -o /home/j/jpryan7/slurm_logs/%x-%j.out \
#       julia/spatial_model/calibration/slurm/reeval_code.sh --theta julia/spatial_model/calibration/runs/exact-polish-base \
#       --out julia/spatial_model/calibration/runs/ccs-20261026/inert/pre
#SBATCH -J sp-reeval
#SBATCH -p econ-grad
#SBATCH -t 03:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu 4G
#SBATCH --export=ALL
#SBATCH -o /home/j/jpryan7/slurm_logs/%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
root="${CODE:-.}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
echo "code root: $root"
exec julia --project=julia/spatial_model/calibration/estimation \
    "$root/julia/spatial_model/calibration/spatial_diagnose.jl" --workers 0 "$@"
