#!/bin/bash
# Table 5 sensitivity panel: one Levenberg–Marquardt refit (spatial_diagnose.jl --polish)
# per case of panel_cases.tsv, started from the baseline fit, as a throttled array job.
# Submit from the repository root with the baseline ϑ (a theta.toml, or a directory
# holding one) as the first argument; extra arguments pass through:
#   n=$(grep -vc '^#' julia/spatial_model/calibration/slurm/panel_cases.tsv)
#   sbatch --array=1-$n%8 julia/spatial_model/calibration/slurm/panel.sh julia/spatial_model/calibration/runs/<fit>
# Results go to julia/spatial_model/calibration/runs/panel/<case>/ (report.txt, theta.toml,
# moments.toml). Other case lists and output roots (the 2026-10 experiments): set
# PANEL_CASES and PANEL_OUT, e.g.
#   PANEL_CASES=julia/spatial_model/calibration/slurm/cases_exact_panel.tsv \
#   PANEL_OUT=julia/spatial_model/calibration/runs/panel-exact sbatch --array=1-$n%16 ... panel.sh <fit>
#SBATCH -J sp-panel
#SBATCH -p econ-grad
#SBATCH -t 08:00:00
#SBATCH -c 17
#SBATCH --mem-per-cpu 2G
#SBATCH --export=ALL
#SBATCH -o julia/spatial_model/calibration/slurm/logs/%x-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
base="$1"; shift
cases="${PANEL_CASES:-julia/spatial_model/calibration/slurm/panel_cases.tsv}"
outroot="${PANEL_OUT:-julia/spatial_model/calibration/runs/panel}"
line=$(grep -v '^#' "$cases" | sed -n "${SLURM_ARRAY_TASK_ID}p")
name=$(cut -f1 <<< "$line")
read -r -a case_args <<< "$(cut -f2 <<< "$line")"
echo "case $SLURM_ARRAY_TASK_ID: $name: ${case_args[*]}"
exec julia --project=julia/spatial_model/calibration/estimation julia/spatial_model/calibration/spatial_diagnose.jl \
    --theta "$base" --polish --max-iter "${PANEL_MAX_ITER:-30}" --consistency --out "$outroot/$name" \
    "${case_args[@]}" "$@"
