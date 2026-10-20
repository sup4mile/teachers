#!/bin/bash
# NLSY79 probability of ever having a child (data/spatial/nlsy_fertility.py): cuts the fertility
# variables out of the 6.7 GB full-cohort CSV in the zip (one core, ~10 GB RAM), then runs the
# estimates and the household-cluster bootstrap. Submit from the repository root, e.g.
#   rdmount   # on the login node first (refresh the Kerberos ticket for ResearchDrive)
#   sbatch julia/spatial_model/calibration/slurm/nlsy_fertility.sh --extract
# Extra arguments go to the script (--extract, --B N, --quick). Log is on NFS.
#SBATCH -J nlsy-fertility
#SBATCH -p econ-grad
#SBATCH -t 02:00:00
#SBATCH -c 1
#SBATCH --mem 16G
#SBATCH --export=ALL
#SBATCH -o /home/j/jpryan7/slurm_logs/%x-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export PATH="/home/j/jpryan7/.local/bin:$PATH"
exec uv run --project data/spatial python data/spatial/nlsy_fertility.py "$@"
